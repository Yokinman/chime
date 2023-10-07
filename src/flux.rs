//! ...

use std::cmp::Ordering;
use std::marker::PhantomData;
use std::vec::IntoIter;
use time::{Time, TimeUnit};
use crate::change::{LinearIso, LinearValue, Scalar};
use crate::degree::{IsBelowDeg, IsDeg, MaxDeg};
use crate::polynomial::{Poly, Roots};

 // Convenient implementations (Vec<T>, RefCell<T>, etc.):
mod flux_impl;
pub use self::flux_impl::*;

/// A value that can change over time.
pub trait FluxValue: Sized {
	type Iso: LinearIso;
	type Degree: IsDeg;
	
	/// Initial value.
	fn value(&self) -> <Self::Iso as LinearIso>::Linear; // ??? Unsure if this should return a reference or not
	
	/// Accumulate change over time.
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>);
	
	/// Apply change over time.
	fn update(&mut self, time: Time); 
	
	fn calc(&self, time: Time) -> <Self::Iso as LinearIso>::Linear {
		let mut changes = ChangeAccum::new(
			ChangeAccumArgs::Sum((time >> TimeUnit::Nanosecs) as f64));
		self.change(&mut changes);
		self.value() + changes.sum
	}
	
	fn at(&self, time: Time) -> Self::Iso {
		Self::Iso::map(self.calc(time))
	}
	
	fn poly(&self) -> Poly<<Self::Iso as LinearIso>::Linear, Self::Degree> {
		let constant = coeff_calc(self, 0.0, 0.0);
		let mut coeff_list = Self::Degree::array_from(
			<Self::Iso as LinearIso>::Linear::ZERO);
		
		for degree in 1..=Self::Degree::USIZE {
			coeff_list[degree - 1] = coeff_calc(self, 0.0, degree as f64);
		}
		
		Poly(constant, coeff_list)
	}
	
	fn when<'v, T>(&'v self, cmp_order: Ordering, cmp_value: &'v T)
		-> When<'v, Self, T>
	where
		T: FluxValue<Iso=Self::Iso>,
		Self::Degree: MaxDeg<T::Degree>,
		<Self::Iso as LinearIso>::Linear: Roots<<Self::Degree as MaxDeg<T::Degree>>::Max> + PartialOrd,
	{
		When {
			value: self,
			cmp_value,
			cmp_order,
		}
	}
	
	fn when_eq<'v, T>(&'v self, eq_value: &'v T)
		-> WhenEq<'v, Self, T>
	where
		T: FluxValue<Iso=Self::Iso>,
		Self::Degree: MaxDeg<T::Degree>,
		<Self::Iso as LinearIso>::Linear: Roots<<Self::Degree as MaxDeg<T::Degree>>::Max> + PartialEq,
	{
		WhenEq {
			value: self,
			eq_value,
		}
	}
}

fn coeff_calc<T: FluxValue>(value: &T, depth: f64, degree: f64)
	-> <T::Iso as LinearIso>::Linear
{
	if depth == 0.0 && degree == 0.0 {
		value.value()
	} else {
		let mut accum = ChangeAccum::new(ChangeAccumArgs::Coeff(degree));
		accum.depth = depth;
		value.change(&mut accum);
		let coeff_sum = accum.sum * Scalar(1.0 / (depth + 1.0));
		if depth >= degree {
			value.value() + coeff_sum
		} else {
			coeff_sum
		}
	}
	// https://www.desmos.com/calculator/tn0udn7swy
}

/// [`FluxValue`] predictive comparison.
#[derive(Copy, Clone, Debug)]
pub struct When<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialOrd
{
	value: &'v A,
	cmp_value: &'v B,
	cmp_order: Ordering,
}

impl<'v, A, B> When<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialOrd
{
	fn time_iter(&self) -> IntoIter<(Time, Time)> {
		let poly = self.value.poly() - self.cmp_value.poly();
		
		 // Find Initial Order:
		let mut degree = <<A::Degree as MaxDeg<B::Degree>>::Max as IsDeg>::USIZE;
		let initial_order = loop {
			let coeff = poly.term(degree).unwrap();
			let order = coeff.partial_cmp(&(coeff * Scalar(0.0)));
			if degree == 0 || order != Some(Ordering::Equal) {
				break if degree % 2 == 0 {
					order
				} else {
					order.map(|o| o.reverse())
				}
			}
			degree -= 1;
		};
		
		 // Convert Roots to Ranges:
		let range_list = match poly.real_roots() {
			Ok(roots) => {
				let mut list = Vec::with_capacity(
					if self.cmp_order == Ordering::Equal {
						roots.len()
					} else {
						1 + (roots.len() / 2)
					}
				);
				let mut prev_point = if initial_order == Some(self.cmp_order) {
					Some(f64::NEG_INFINITY)
				} else {
					None
				};
				for point in roots {
					if let Some(prev) = prev_point {
						if point != prev {
							list.push((prev, point));
						}
						prev_point = None;
					} else if self.cmp_order == Ordering::Equal {
						list.push((point, point));
					} else {
						prev_point = Some(point);
					}
				}
				if let Some(point) = prev_point {
					list.push((point, f64::INFINITY));
				}
				list
			},
			Err(_) => vec![],
		};
		
		 // Convert Ranges to Time:
		let vec: Vec<(Time, Time)> = range_list.into_iter()
			.filter_map(|(a, b)| {
				assert!(a <= b);
				if b < 0.0 || a > (u64::MAX as f64) {
					None
				} else {
					Some((
						(a as u64)*TimeUnit::Nanosecs,
						(b as u64)*TimeUnit::Nanosecs
					))
				}
			})
			.collect();
		
		vec.into_iter()
	}
}

impl<'v, A, B> IntoIterator for When<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialOrd
{
	type Item = (Time, Time);
	type IntoIter = std::vec::IntoIter<Self::Item>;
	
	fn into_iter(self) -> Self::IntoIter {
		self.time_iter()
	}
}

/// [`FluxValue`] predictive equality comparison.
#[derive(Copy, Clone, Debug)]
pub struct WhenEq<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialEq
{
	value: &'v A,
	eq_value: &'v B,
}

impl<'v, A, B> WhenEq<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialEq
{
	fn time_iter(&self) -> IntoIter<Time> {
		let poly = self.value.poly() - self.eq_value.poly();
		let mut real_roots = poly.real_roots().unwrap_or(vec![]);
		
		 // Constant Equality:
		if
			real_roots.is_empty()
			&& poly.constant() == <A::Iso as LinearIso>::Linear::ZERO
			&& poly.coeff_iter().all(|&term| term == <A::Iso as LinearIso>::Linear::ZERO)
		{
			real_roots.push(0.0);
		}
		
		 // Convert Roots to Times:
		let vec: Vec<Time> = real_roots.into_iter()
			.filter_map(|t| {
				// let t = t + offset;
				if t < 0.0 || t > (u64::MAX as f64) {
					None
				} else {
					Some((t as u64)*TimeUnit::Nanosecs)
				}
			})
			.collect();
		
		vec.into_iter()
	}
}

impl<'v, A, B> IntoIterator for WhenEq<'v, A, B>
where
	A: FluxValue,
	B: FluxValue<Iso=A::Iso>,
	A::Degree: MaxDeg<B::Degree>,
	<A::Iso as LinearIso>::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialEq
{
	type Item = Time;
	type IntoIter = std::vec::IntoIter<Self::Item>;
	
	fn into_iter(self) -> Self::IntoIter {
		self.time_iter()
	}
}

/// Change accumulator.
#[derive(Debug)]
pub struct ChangeAccum<I: LinearIso, D: IsDeg> {
	sum: I::Linear,
	args: ChangeAccumArgs,
	depth: f64,
	degree: PhantomData<D>,
}

impl<I: LinearIso, D: IsDeg> ChangeAccum<I, D> {
	fn new(args: ChangeAccumArgs) -> Self {
		Self {
			sum: <I::Linear as LinearValue>::ZERO,
			args,
			depth: 0.0,
			degree: PhantomData,
		}
	}
	
	pub fn add<C>(&mut self, value: &C, unit: TimeUnit)
	where
		C: FluxValue<Iso=I>,
		C::Degree: IsBelowDeg<D>,
	{
		self.accum(Scalar(1.0), value, unit);
	}
	
	pub fn sub<C>(&mut self, value: &C, unit: TimeUnit)
	where
		C: FluxValue<Iso=I>,
		C::Degree: IsBelowDeg<D>,
	{
		self.accum(Scalar(-1.0), value, unit);
	}
	
	fn accum<C>(&mut self, scalar: Scalar, change: &C, unit: TimeUnit)
	where
		C: FluxValue<Iso=I>,
		C::Degree: IsBelowDeg<D>,
	{
		self.sum = self.sum + match &mut self.args {
			ChangeAccumArgs::Sum(time) => {
				 // Fetch Value of Change:
				let mut changes = ChangeAccum::new(ChangeAccumArgs::Sum(*time));
				changes.depth = self.depth + 1.0;
				change.change(&mut changes);
				let value = change.value() + changes.sum;
				
				let time = *time * ((unit >> TimeUnit::Nanosecs) as f64).recip();
				value * Scalar((time + self.depth) / (self.depth + 1.0)) * scalar
			},
			ChangeAccumArgs::Coeff(degree) => {
				let coeff = coeff_calc(change, self.depth + 1.0, *degree);
				if self.depth < *degree {
					let coeff = coeff * Scalar(((unit >> TimeUnit::Nanosecs) as f64).recip());
					if self.depth == 0.0 {
						coeff * scalar
					} else {
						(coeff + (coeff_calc(change, self.depth + 1.0, *degree + 1.0) * Scalar(self.depth))) * scalar
					}
				} else {
					coeff * Scalar(self.depth) * scalar
				}
			}
		};
	}
}

#[derive(Debug)]
enum ChangeAccumArgs {
	Sum(f64),
	Coeff(f64),
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	use crate::change::*;
	use crate::degree::*;
	
	#[derive(Debug, Default)] struct Pos { value: f64, spd: Spd, misc: Vec<Spd> }
	#[derive(Debug, Default)] struct Spd { value: f64, fric: Fric, accel: Accel }
	#[derive(Debug, Default)] struct Fric { value: f64 }
	#[derive(Debug, Default)] struct Accel { value: f64, jerk: Jerk }
	#[derive(Debug, Default)] struct Jerk { value: f64, snap: Snap }
	#[derive(Debug, Default)] struct Snap { value: f64 }
	
	impl FluxValue for Pos {
		type Iso = Linear<i64>; // ??? Iso could be generic, allowing different types of the same isomorphism to be compared
		type Degree = Deg<4>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
			changes.add(&self.spd, TimeUnit::Secs);
			changes.add(&self.misc, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.spd.update(time);
		}
	}
	
	impl FluxValue for Spd {
		type Iso = Linear<i64>;
		type Degree = Deg<3>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
			changes.sub(&self.fric, TimeUnit::Secs);
			changes.add(&self.accel, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.fric.update(time);
			self.accel.update(time);
		}
	}
	
	impl FluxValue for Fric {
		type Iso = Linear<i64>;
		type Degree = Deg<0>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, _changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	impl FluxValue for Accel {
		type Iso = Linear<i64>;
		type Degree = Deg<2>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
			changes.add(&self.jerk, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.jerk.update(time);
		}
	}
	
	impl FluxValue for Jerk {
		type Iso = Linear<i64>;
		type Degree = Deg<1>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
			changes.add(&self.snap, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.snap.update(time);
		}
	}
	
	impl FluxValue for Snap {
		type Iso = Linear<i64>;
		type Degree = Deg<0>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, _changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	fn position() -> Pos {
		Pos {
			value: 32.0, 
			spd: Spd {
				value: 0.0,
				fric: Fric { value: 3.5, },
				accel: Accel {
					value: 0.3,
					jerk: Jerk {
						value: 0.4,
						snap: Snap { value: -0.01 },
					},
				},
			},
			misc: Vec::new(),
		}
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		
		 // Values:
		assert_eq!(pos.at(0*Secs), Linear::from(32));
		assert_eq!(pos.at(10*Secs), Linear::from(-63));
		assert_eq!(pos.at(20*Secs), Linear::from(-113));
		assert_eq!(pos.at(100*Secs), Linear::from(8339));
		assert_eq!(pos.at(200*Secs), Linear::from(-209779));
		
		 // Update:
		pos.update(20*Secs);
		assert_eq!(pos.at(0*Secs), Linear::from(-113));
		assert_eq!(pos.at(80*Secs), Linear::from(8339));
		assert_eq!(pos.at(180*Secs), Linear::from(-209779));
		pos.update(55*Secs);
		assert_eq!(pos.at(25*Secs), Linear::from(8339));
		assert_eq!(pos.at(125*Secs), Linear::from(-209779));
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_eq!(
			pos.poly(),
			Poly(32.0, [
				-1.469166666666667e-9,
				-1.4045833333333335e-18,
				0.06416666666666669e-27,
				-0.0004166666666666668e-36
			])
		);
		pos.update(20*Secs);
		assert_eq!(
			pos.poly(),
			Poly(-112.55000000000007, [
				6.0141666666666615e-9,
				1.4454166666666668e-18,
				0.030833333333333334e-27,
				-0.0004166666666666668e-36
			])
		);
	}
	
	#[test]
	fn when() {
		let mut pos = position();
		pos.update(20*Secs);
		
		let vec: Vec<(Time, Time)> = pos.when(Ordering::Greater, &position())
			.into_iter().collect();
		
		let vec_eq: Vec<Time> = pos.when_eq(&position())
			.into_iter().collect();
		
		assert_eq!(vec, [(6110872304*Nanosecs, 87499325334*Nanosecs)]);
		assert_eq!(vec_eq, [6110872304*Nanosecs, 87499325334*Nanosecs]);
	}
	
	#[test]
	fn when_refcell() {
		use std::cell::RefCell;
		
		let a_pos = RefCell::new(position());
		let b_pos = RefCell::new(position());
		let wh = a_pos.when(Ordering::Greater, &b_pos);
		let wh_eq = a_pos.when_eq(&b_pos);
		
		 // Check Before:
		let vec: Vec<(Time, Time)> = wh.time_iter().collect();
		let vec_eq: Vec<Time> = wh_eq.time_iter().collect();
		assert_eq!(vec, []);
		assert_eq!(vec_eq, [0*Nanosecs]);
		
		 // Apply Changes:
		a_pos.borrow_mut().update(20*Secs);
		b_pos.borrow_mut().misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.borrow_mut().misc.push(Spd {
			value: 12.25,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		
		 // Check After:
		let vec: Vec<(Time, Time)> = wh.time_iter().collect();
		let vec_eq: Vec<Time> = wh_eq.time_iter().collect();
		assert_eq!(vec, [(8517857837*Nanosecs, 90130683345*Nanosecs)]);
		assert_eq!(vec_eq, [8517857837*Nanosecs, 90130683345*Nanosecs]);
	}
}