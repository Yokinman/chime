//! Utilities for describing how a type changes over time.

use std::borrow::Borrow;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Deref, DerefMut, Shr, Sub};
use std::rc::Rc;

use time::{Time, TimeUnit};
use crate::linear::*;
use crate::poly::*;

mod impls;
mod kind;
pub use self::impls::*;
pub use self::kind::*;

/// A value that can change over time.
pub trait FluxValue: Sized {
	/// The produced value.
	type Value;
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// The output accumulator of [`FluxValue::change`].
	type OutAccum<'a>: FluxAccum<'a, Self::Kind>
	where <Self::Kind as FluxKind>::Linear: 'a;
	
	/// Initial value.
	fn value(&self) -> <Self::Kind as FluxKind>::Linear;
	fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear);
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> Self::OutAccum<'a>;
	
	/// Apply change over time.
	fn advance(&mut self, time: Time);
	
	fn calc_value(&self, time: Time) -> <Self::Kind as FluxKind>::Linear {
		let mut value = self.value();
		let value_accum = FluxAccumKind::Sum { sum: &mut value, depth: 0, time };
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(value_accum));
		value
	}
	
	fn at(&self, time: Time) -> Self::Value
	where
		<Self::Kind as FluxKind>::Linear: LinearIso<Self::Value>,
	{
		self.calc_value(time).map()
	}
	
	fn poly(&self) -> Poly<Self::Kind> {
		let mut poly = Poly::default();
		poly.0 = self.value();
		let poly_accum = FluxAccumKind::Poly { poly: &mut poly, depth: 0 };
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(poly_accum));
		poly
	}
}

/// ...
#[allow(type_alias_bounds)]
type PredPoly<A: FluxValue> = Rc<RefCell<(
	Time,
	Poly<A::Kind>
	// ??? It may be worth storing the original value so that the polynomial
	// can be easily offset for comparison with the other value.
)>>;

/// A value with a predictable change over time.
pub struct PredValue<A: FluxValue> {
	value: A,
	poly: PredPoly<A>,
	time: Time,
}

impl<A: FluxValue> PredValue<A> {
	pub fn new(value: A) -> Self {
		let time = Time::default();
		let poly = Rc::new(RefCell::new((time, value.poly())));
		Self { value, poly, time }
	}
	
	/// Compares to any [`FluxValue`] or [`PredValue`] reference. If compared to
	/// a [`PredValue`], the returned [`When`] will update its predictions after
	/// any modifications to its inner flux value. 
	pub fn when<'a, B>(&'a self, cmp_order: Ordering, other: impl Into<PredValue<&'a B>>) -> When<A, B>
	where
		B: FluxValue + 'a,
		When<A, B>: IntoIterator
	{
		When {
			a_poly: self.poly.clone(),
			b_poly: other.into().poly,
			cmp_order,
		}
	}
	
	/// Compares to any [`FluxValue`] or [`PredValue`] reference. If compared to
	/// a [`PredValue`], the returned [`WhenEq`] will update its predictions
	/// after any modifications to its inner flux value.
	pub fn when_eq<'a, B>(&'a self, other: impl Into<PredValue<&'a B>>) -> WhenEq<A, B>
	where
		B: FluxValue + 'a,
		WhenEq<A, B>: IntoIterator
	{
		WhenEq {
			a_poly: self.poly.clone(),
			b_poly: other.into().poly,
		}
	}
	
	pub fn borrow_mut(&mut self) -> PredValueMut<A> {
		PredValueMut(self)
	}
	
	pub fn advance(&mut self, time: Time) {
		self.time += time;
		self.value.advance(time);
		// let (last_time, _) = *RefCell::borrow(&self.poly).last().unwrap();
		// self.poly.borrow_mut().push((last_time + time, self.value.poly()));
	}
}

impl<A: FluxValue> From<A> for PredValue<A> {
	fn from(value: A) -> Self {
		Self::new(value)
	}
}

impl<'a, A: FluxValue + 'a> From<&'a PredValue<A>> for PredValue<&'a A> {
	fn from(value: &'a PredValue<A>) -> Self {
		PredValue {
			value: &value.value,
			poly: value.poly.clone(),
			time: value.time,
		}
	}
}

impl<A: FluxValue> Deref for PredValue<A> {
	type Target = A;
	fn deref(&self) -> &Self::Target {
		&self.value
	}
}

impl<A: FluxValue> Borrow<A> for PredValue<A> {
	fn borrow(&self) -> &A {
		&self.value
	}
}

impl<A: FluxValue> Debug for PredValue<A> where A: Debug {
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		self.value.fmt(f)
	}
}

impl<A: FluxValue> Display for PredValue<A> where A: Display {
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		self.value.fmt(f)
	}
}

/// Borrowed mutable access to the value contained by [`PredValue`], which
/// updates any `When` or `WhenEq` predictions on drop. 
pub struct PredValueMut<'a, A: FluxValue>(&'a mut PredValue<A>);

impl<A: FluxValue> Drop for PredValueMut<'_, A> {
	fn drop(&mut self) {
		let (ref mut time, ref mut poly) = *self.0.poly.borrow_mut();
		*time = self.0.time;
		*poly = self.0.value.poly();
	}
}

impl<'a, A: FluxValue + 'a> Deref for PredValueMut<'a, A> {
	type Target = A;
	fn deref(&self) -> &Self::Target {
		&self.0.value
	}
}

impl<A: FluxValue> DerefMut for PredValueMut<'_, A> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0.value
	}
}

impl<A: FluxValue> Debug for PredValueMut<'_, A> where A: Debug {
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		self.0.value.fmt(f)
	}
}

impl<A: FluxValue> Display for PredValueMut<'_, A> where A: Display {
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		self.0.value.fmt(f)
	}
}

/// Iterator of [`Time`] ranges.
#[must_use]
pub struct TimeRanges(std::vec::IntoIter<(Time, Time)>);

impl Iterator for TimeRanges {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		self.0.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.0.size_hint()
	}
	fn count(self) -> usize {
		self.0.count()
	}
}

/// Iterator of [`Time`] values.
#[must_use]
pub struct Times(std::vec::IntoIter<Time>);

impl Iterator for Times {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		self.0.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.0.size_hint()
	}
	fn count(self) -> usize {
		self.0.count()
	}
}

/// [`FluxValue`] predictive comparison.
#[derive(Debug)]
pub struct When<A: FluxValue, B: FluxValue> {
	a_poly: PredPoly<A>,
	b_poly: PredPoly<B>,
	cmp_order: Ordering,
}

impl<A: FluxValue, B: FluxValue> Clone for When<A, B> {
	fn clone(&self) -> Self {
		Self {
			a_poly: self.a_poly.clone(),
			b_poly: self.b_poly.clone(),
			cmp_order: self.cmp_order,
		}
	}
}

impl<A, B> IntoIterator for When<A, B>
where
	A: FluxValue,
	B: FluxValue,
	B::Kind: FluxKind<Linear = <A::Kind as FluxKind>::Linear>,
	A::Kind: Add<B::Kind>,
	<A::Kind as FluxKind>::Linear: PartialOrd,
	<A::Kind as Add<B::Kind>>::Output: FluxKind<Linear = <A::Kind as FluxKind>::Linear> + Roots + PartialOrd,
	Poly<A::Kind>: Sub<Poly<B::Kind>, Output=Poly<<A::Kind as Add<B::Kind>>::Output>>,
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = RefCell::borrow(&self.a_poly).1
			- RefCell::borrow(&self.b_poly).1;
		
		 // Find Initial Order:
		let mut degree = <<A::Kind as Add<B::Kind>>::Output as FluxKind>::DEGREE;
		let initial_order = loop {
			if degree == 0 {
				break poly.constant().partial_cmp(&<A::Kind as FluxKind>::Linear::zero());
			} else {
				let coeff = poly.coeff(degree - 1).unwrap();
				let order = coeff.partial_cmp(&FluxKind::zero());
				if order != Some(Ordering::Equal) {
					break if degree % 2 == 0 {
						order
					} else {
						order.map(|o| o.reverse())
					}
				}
				degree -= 1;
			}
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
				for &point in roots.into_iter() {
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
				list.into_iter()
			},
			Err(_) => Default::default(),
		};
		
		 // Convert Ranges to Times:
		let list: Vec<(Time, Time)> = range_list
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
		
		TimeRanges(list.into_iter())
	}
}

/// [`FluxValue`] predictive equality comparison.
#[derive(Debug)]
pub struct WhenEq<A: FluxValue, B: FluxValue> {
	a_poly: PredPoly<A>,
	b_poly: PredPoly<B>,
}

impl<A, B> Clone for WhenEq<A, B>
where
	A: FluxValue,
	B: FluxValue,
{
	fn clone(&self) -> Self {
		Self {
			a_poly: self.a_poly.clone(),
			b_poly: self.b_poly.clone(),
		}
	}
}

impl<A, B> IntoIterator for WhenEq<A, B>
where
	A: FluxValue,
	B: FluxValue,
	B::Kind: FluxKind<Linear = <A::Kind as FluxKind>::Linear>,
	A::Kind: Add<B::Kind>,
	<A::Kind as FluxKind>::Linear: PartialEq,
	<A::Kind as Add<B::Kind>>::Output: FluxKind<Linear = <A::Kind as FluxKind>::Linear> + Roots + PartialEq,
	Poly<A::Kind>: Sub<Poly<B::Kind>, Output=Poly<<A::Kind as Add<B::Kind>>::Output>>,
{
	type Item = Time;
	type IntoIter = Times;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = RefCell::borrow(&self.a_poly).1
			- RefCell::borrow(&self.b_poly).1;
		let mut real_roots = poly.real_roots().unwrap_or(Box::default())
			.into_vec();
		
		 // Constant Equality:
		if
			real_roots.is_empty()
			&& poly.constant().is_zero()
			&& poly.coeff_iter().all(FluxKind::is_zero)
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
		
		Times(vec.into_iter())
	}
}

#[allow(type_alias_bounds)]
pub type Changes<'a, T: FluxValue> = <T::Kind as FluxKind>::Accum<'a>;

/// Change accumulator.
pub trait FluxAccum<'a, K: FluxKind> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self;
}

/// Nested summation change accumulator.
pub struct SumAccum<'a, K: FluxKind>(FluxAccumKind<'a, K>);

impl<'a, K: FluxKind> FluxAccum<'a, K> for SumAccum<'a, K> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self {
		Self(kind)
	}
}

impl<K: FluxKind> SumAccum<'_, K>
where
	K: Add<Output=K>,
{
	pub fn add<C: FluxValue>(self, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		self.accum(Scalar(1.0), value, unit)
	}
	
	pub fn sub<C: FluxValue>(self, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		self.accum(Scalar(-1.0), value, unit)
	}
	
	fn accum<C: FluxValue>(mut self, scalar: Scalar, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		match &mut self.0 {
			FluxAccumKind::Sum { sum, depth, time } => {
				let mut sub_sum = value.value();
				value.change(<C::Kind as FluxKind>::Accum::from_kind(FluxAccumKind::Sum {
					sum: &mut sub_sum,
					depth: *depth + 1,
					time: *time,
				}));
				let depth = *depth as f64;
				let time_scale = (time.as_nanos() as f64) / ((unit >> TimeUnit::Nanosecs) as f64);
				**sum = **sum + (sub_sum * Scalar((time_scale + depth) / (depth + 1.0)) * scalar);
			},
			FluxAccumKind::Poly { poly, depth } => {
				let mut sub_poly = Poly::default();
				sub_poly.0 = value.value();
				value.change(<C::Kind as FluxKind>::Accum::from_kind(FluxAccumKind::Poly {
					poly: &mut sub_poly,
					depth: *depth + 1,
				}));
				let depth = *depth as f64;
				let unit_scale = ((unit >> TimeUnit::Nanosecs) as f64).recip();
				**poly = **poly
					+ (sub_poly * (Scalar(depth / (depth + 1.0)) * scalar))
					+ ((sub_poly >> DegShift) * (Scalar(unit_scale / (depth + 1.0)) * scalar));
				// https://www.desmos.com/calculator/mhlpjakz32
			}
		}
		self
	}
}

/// General accumulator arguments.
#[non_exhaustive]
pub enum FluxAccumKind<'a, K: FluxKind> {
	Sum {
		sum: &'a mut K::Linear,
		depth: usize,
		time: Time,
	},
	Poly {
		poly: &'a mut Poly<K>,
		depth: usize,
	},
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	
	#[derive(Debug, Default)] struct Pos { value: f64, spd: Spd, misc: Vec<Spd> }
	#[derive(Debug, Default)] struct Spd { value: f64, fric: Fric, accel: Accel }
	#[derive(Debug, Default)] struct Fric { value: f64 }
	#[derive(Debug, Default)] struct Accel { value: f64, jerk: Jerk }
	#[derive(Debug, Default)] struct Jerk { value: f64, snap: Snap }
	#[derive(Debug, Default)] struct Snap { value: f64 }
	
	impl FluxValue for Pos {
		type Value = i64;
		type Kind = Deg<f64, 4>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.add(&self.spd, TimeUnit::Secs)
				.add(&self.misc, TimeUnit::Secs)
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.spd.advance(time);
		}
	}
	
	impl FluxValue for Spd {
		type Value = i64;
		type Kind = Deg<f64, 3>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.sub(&self.fric, TimeUnit::Secs)
				.add(&self.accel, TimeUnit::Secs)
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.fric.advance(time);
			self.accel.advance(time);
		}
	}
	
	impl FluxValue for Fric {
		type Value = i64;
		type Kind = Deg<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
		}
	}
	
	impl FluxValue for Accel {
		type Value = i64;
		type Kind = Deg<f64, 2>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.jerk, TimeUnit::Secs)
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.jerk.advance(time);
		}
	}
	
	impl FluxValue for Jerk {
		type Value = i64;
		type Kind = Deg<f64, 1>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.snap, TimeUnit::Secs)
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.snap.advance(time);
		}
	}
	
	impl FluxValue for Snap {
		type Value = i64;
		type Kind = Deg<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
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
		assert_eq!(pos.at(0*Secs), 32);
		assert_eq!(pos.at(10*Secs), -63);
		assert_eq!(pos.at(20*Secs), -113);
		assert_eq!(pos.at(100*Secs), 8339);
		assert_eq!(pos.at(200*Secs), -209779);
		
		 // Update:
		pos.advance(20*Secs);
		assert_eq!(pos.at(0*Secs), -113);
		assert_eq!(pos.at(80*Secs), 8339);
		assert_eq!(pos.at(180*Secs), -209779);
		pos.advance(55*Secs);
		assert_eq!(pos.at(25*Secs), 8339);
		assert_eq!(pos.at(125*Secs), -209779);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_eq!(
			pos.poly(),
			Poly(32.0, [
				Deg(-1.469166666666667e-9),
				Deg(-1.4045833333333335e-18),
				Deg(0.06416666666666669e-27),
				Deg(-0.00041666666666666684e-36),
			])
		);
		pos.advance(20*Secs);
		assert_eq!(
			pos.poly(),
			Poly(-112.55000000000007, [
				Deg(6.0141666666666615e-9),
				Deg(1.445416666666667e-18),
				Deg(0.03083333333333335e-27),
				Deg(-0.00041666666666666684e-36),
			])
		);
	}
	
	#[test]
	fn when() {
		let pos = PredValue::new(position());
		let acc = position().spd.accel;
		
		let pos_when = pos.when(Ordering::Greater, &acc);
		let vec: Vec<(Time, Time)> = pos_when.clone().into_iter().collect();
		
		let vec_eq: Vec<Time> = pos.when_eq(&acc)
			.into_iter().collect();
		
		assert_eq!(vec, [
			(0*Nanosecs, 4560099744*Nanosecs),
			(26912076290*Nanosecs, 127394131312*Nanosecs)
		]);
		assert_eq!(vec_eq, [
			4560099744*Nanosecs,
			26912076290*Nanosecs,
			127394131312*Nanosecs
		]);
		
		// pos.advance(20*Secs);
		// pos.borrow_mut().spd.value = -20.0;
		// pos.borrow_mut().spd.accel.jerk.value = 0.3;
		// 
		// assert_eq!(pos_when.into_iter().collect::<Vec<(Time, Time)>>(), [
		// 	(0*Nanosecs, 4560099744*Nanosecs),
		// 	(33544693273*Nanosecs, 157824014330*Nanosecs)
		// ]);
		// assert_eq!(pos.when(Ordering::Greater, &Snap { value: -116.0 }).into_iter().collect::<Vec<(Time, Time)>>(), [
		// 	(0*Nanosecs, 16114479737*Nanosecs),
		// 	(19315363058*Nanosecs, 20188850711*Nanosecs),
		// 	(29523688931*Nanosecs, 157716291466*Nanosecs)
		// ]);
	}
	
	#[test]
	fn when_borrow() {
		let mut a_pos = PredValue::new(position());
		let mut b_pos = PredValue::new(position());
		let wh = a_pos.when(Ordering::Greater, &b_pos);
		let wh_eq = a_pos.when_eq(&b_pos);
		
		 // Check Before:
		let vec: Vec<(Time, Time)> = wh.clone().into_iter().collect();
		let vec_eq: Vec<Time> = wh_eq.clone().into_iter().collect();
		assert_eq!(vec, []);
		assert_eq!(vec_eq, [0*Nanosecs]);
		
		 // Apply Changes:
		a_pos.borrow_mut().advance(20*Secs);
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
		let vec: Vec<(Time, Time)> = wh.into_iter().collect();
		let vec_eq: Vec<Time> = wh_eq.into_iter().collect();
		assert_eq!(vec, [(8517857837*Nanosecs, 90130683345*Nanosecs)]);
		assert_eq!(vec_eq, [8517857837*Nanosecs, 90130683345*Nanosecs]);
	}
}