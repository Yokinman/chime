//! ...

use std::borrow::Borrow;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};
use std::rc::Rc;
use time::{Time, TimeUnit};
use crate::change::{LinearIso, LinearValue, Scalar};
use crate::degree::{IsBelowDeg, IsDeg, MaxDeg};
use crate::polynomial::{Poly, Roots};

 // Convenient implementations (Vec<T>, ??? tuples, etc.):
mod flux_impl;
pub use self::flux_impl::*;

/// A value that can change over time.
pub trait FluxValue: Sized {
	/// The produced value.
	type Value;
	
	/// The inner linear isomorphic value (linear (+), exponential (*), etc.).
	type Linear: LinearValue + LinearIso<Self::Value>;
	
	/// The kind of change over time.
	type Degree: IsDeg;
	
	/// Initial value.
	fn value(&self) -> Self::Linear;
	fn set_value(&mut self, value: Self::Linear);
	
	/// Accumulates change over time.
	fn change(&self, changes: &mut Changes<Self>);
	
	/// Apply change over time.
	fn advance(&mut self, time: Time);
	
	fn calc_value(&self, time: Time) -> Self::Linear {
		let mut changes = ChangeAccum::new(
			ChangeAccumArgs::Sum(time.as_nanos() as f64));
		self.change(&mut changes);
		self.value() + changes.sum
	}
	
	fn at(&self, time: Time) -> Self::Value {
		self.calc_value(time).map()
	}
	
	fn poly(&self) -> Poly<Self::Linear, Self::Degree> {
		let constant = coeff_calc(self, 0.0, 0.0);
		let mut coeff_list = Self::Degree::array_from(Self::Linear::zero());
		
		for degree in 1..=Self::Degree::USIZE {
			coeff_list[degree - 1] = coeff_calc(self, 0.0, degree as f64);
		}
		
		Poly(constant, coeff_list)
	}
}

fn coeff_calc<T: FluxValue>(value: &T, depth: f64, degree: f64) -> T::Linear {
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

/// ...
#[allow(type_alias_bounds)]
type PredPoly<A: FluxValue> = Rc<RefCell<(
	Time,
	Poly<A::Linear, A::Degree>
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
		B: FluxValue<Linear=A::Linear> + 'a,
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
		B: FluxValue<Linear=A::Linear> + 'a,
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
pub struct When<A, B>
where
	A: FluxValue,
	B: FluxValue<Linear=A::Linear>,
{
	a_poly: PredPoly<A>,
	b_poly: PredPoly<B>,
	cmp_order: Ordering,
}

impl<A, B> Clone for When<A, B>
where
	A: FluxValue,
	B: FluxValue<Linear=A::Linear>,
{
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
	B: FluxValue<Linear=A::Linear>,
	A::Degree: MaxDeg<B::Degree>,
	A::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialOrd
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = RefCell::borrow(&self.a_poly).1
			- RefCell::borrow(&self.b_poly).1;
		
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
	/*{
		let mut a = RefCell::borrow(&self.a_poly).clone().into_iter();
		let mut b = RefCell::borrow(&self.b_poly).clone().into_iter();
		let (mut a_time, mut a_poly) = a.next().unwrap();
		let (mut b_time, mut b_poly) = b.next().unwrap();
		let mut next_a = a.next();
		let mut next_b = b.next();
		
		let mut list: Vec<(Time, Time)> = Vec::new();
		
		let mut can_loop = true;
		
		while can_loop {
		
		 // Apply Time Offset: ??? Could just store the original value, which wouldn't need extra logic for time offsets.
		fn binom(n: i32, k: i32) -> i32 {
			if n == k || k == 0 {
				return 1
			}
			binom(n - 1, k - 1) + binom(n - 1, k)
		}
		if a_time > b_time {
			let offset = ((a_time - b_time) >> TimeUnit::Nanosecs) as f64;
			let mut n = 1;
			for &coeff in b_poly.clone().coeff_iter() {
				let constant = &mut b_poly.0;
				*constant = *constant + (coeff * Scalar(offset.powi(n)));
				let mut k = 1;
				for mut_coeff in b_poly.coeff_iter_mut() {
					*mut_coeff = *mut_coeff + (coeff * Scalar(offset.powi(n - k) * (binom(n + 1, k) as f64)));
					if n == k {
						break
					}
					k += 1;
				}
				n += 1;
			}
		} else if b_time > a_time {
			let offset = ((b_time - a_time) >> TimeUnit::Nanosecs) as f64;
			let mut n = 1;
			for &coeff in a_poly.clone().coeff_iter() {
				let constant = &mut a_poly.0;
				*constant = *constant + (coeff * Scalar(offset.powi(n)));
				let mut k = 1;
				for mut_coeff in a_poly.coeff_iter_mut() {
					*mut_coeff = *mut_coeff + (coeff * Scalar(offset.powi(n - k) * (binom(n + 1, k) as f64)));
					if n == k {
						break
					}
					k += 1;
				}
				n += 1;
			}
		}
		
		let poly = a_poly - b_poly;
		
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
		
		let min_time = (a_time.max(b_time) >> TimeUnit::Nanosecs) as f64;
		
		match (next_a, next_b) {
			(Some((next_a_time, next_a_poly)), Some((next_b_time, next_b_poly))) => {
				match next_a_time.cmp(&next_b_time) {
					Ordering::Equal => {
						a_time = next_a_time;
						a_poly = next_a_poly;
						b_time = next_b_time;
						b_poly = next_b_poly;
						next_a = a.next();
						next_b = b.next();
					},
					Ordering::Greater => {
						b_time = next_b_time;
						b_poly = next_b_poly;
						next_b = b.next();
					},
					Ordering::Less => {
						a_time = next_a_time;
						a_poly = next_a_poly;
						next_a = a.next();
					},
				}
			},
			(Some((next_a_time, next_a_poly)), None) => {
				a_time = next_a_time;
				a_poly = next_a_poly;
				next_a = a.next();
			},
			(None, Some((next_b_time, next_b_poly))) => {
				b_time = next_b_time;
				b_poly = next_b_poly;
				next_b = b.next();
			},
			(None, None) => can_loop = false,
		}
		
		 // Convert Ranges to Times:
		let max_time = if can_loop {
			(a_time.max(b_time) >> TimeUnit::Nanosecs) as f64
		} else {
			(u64::MAX as f64) + 1.0
		};
		let mut append_list: Vec<(Time, Time)> = range_list
			.filter_map(|(a, b)| {
				assert!(a <= b);
				let a = a + min_time;
				let b = b + min_time;
				if b < min_time || a >= max_time {
					None
				} else {
					Some((
						(a.max(min_time) as u64)*TimeUnit::Nanosecs,
						(b.min(max_time) as u64)*TimeUnit::Nanosecs
					))
				}
			})
			.collect();
			
		if let Some((_, last_time)) = list.last_mut() {
			if let Some(&(first_time, next_time)) = append_list.first() {
				 // Merge adjacent ranges:
				if *last_time == first_time {
					*last_time = next_time;
					append_list.remove(0);
				}
			}
		}
		
		list.append(&mut append_list);
		
		}
		
		TimeRanges(list.into_iter())
	}*/
}

/// [`FluxValue`] predictive equality comparison.
#[derive(Debug)]
pub struct WhenEq<A, B>
where
	A: FluxValue,
	B: FluxValue<Linear=A::Linear>,
{
	a_poly: PredPoly<A>,
	b_poly: PredPoly<B>,
}

impl<A, B> Clone for WhenEq<A, B>
where
	A: FluxValue,
	B: FluxValue<Linear=A::Linear>,
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
	B: FluxValue<Linear=A::Linear>,
	A::Degree: MaxDeg<B::Degree>,
	A::Linear: Roots<<A::Degree as MaxDeg<B::Degree>>::Max> + PartialEq
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
			&& poly.constant() == A::Linear::zero()
			&& poly.coeff_iter().all(|&term| term == A::Linear::zero())
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
pub type Changes<T: FluxValue> = ChangeAccum<T::Value, T::Linear, T::Degree>;

/// Change accumulator.
#[derive(Debug)]
pub struct ChangeAccum<M, I: LinearIso<M>, D: IsDeg> {
	sum: I,
	args: ChangeAccumArgs,
	depth: f64,
	degree: PhantomData<D>,
	mapped: PhantomData<M>,
}

impl<M, I: LinearIso<M>, D: IsDeg> ChangeAccum<M, I, D> {
	fn new(args: ChangeAccumArgs) -> Self {
		Self {
			sum: I::zero(),
			args,
			depth: 0.0,
			degree: PhantomData,
			mapped: PhantomData,
		}
	}
	
	pub fn add<C>(&mut self, value: &C, unit: TimeUnit)
	where
		C: FluxValue<Value=M, Linear=I>,
		C::Degree: IsBelowDeg<D>,
	{
		self.accum(Scalar(1.0), value, unit);
	}
	
	pub fn sub<C>(&mut self, value: &C, unit: TimeUnit)
	where
		C: FluxValue<Value=M, Linear=I>,
		C::Degree: IsBelowDeg<D>,
	{
		self.accum(Scalar(-1.0), value, unit);
	}
	
	fn accum<C>(&mut self, scalar: Scalar, change: &C, unit: TimeUnit)
	where
		C: FluxValue<Value=M, Linear=I>,
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
	use crate::degree::*;
	
	#[derive(Debug, Default)] struct Pos { value: f64, spd: Spd, misc: Vec<Spd> }
	#[derive(Debug, Default)] struct Spd { value: f64, fric: Fric, accel: Accel }
	#[derive(Debug, Default)] struct Fric { value: f64 }
	#[derive(Debug, Default)] struct Accel { value: f64, jerk: Jerk }
	#[derive(Debug, Default)] struct Jerk { value: f64, snap: Snap }
	#[derive(Debug, Default)] struct Snap { value: f64 }
	
	impl FluxValue for Pos {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<4>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, changes: &mut Changes<Self>) {
			changes.add(&self.spd, TimeUnit::Secs);
			changes.add(&self.misc, TimeUnit::Secs);
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.spd.advance(time);
		}
	}
	
	impl FluxValue for Spd {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<3>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, changes: &mut Changes<Self>) {
			changes.sub(&self.fric, TimeUnit::Secs);
			changes.add(&self.accel, TimeUnit::Secs);
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.fric.advance(time);
			self.accel.advance(time);
		}
	}
	
	impl FluxValue for Fric {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<0>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, _changes: &mut Changes<Self>) {}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
		}
	}
	
	impl FluxValue for Accel {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<2>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, changes: &mut Changes<Self>) {
			changes.add(&self.jerk, TimeUnit::Secs);
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.jerk.advance(time);
		}
	}
	
	impl FluxValue for Jerk {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<1>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, changes: &mut Changes<Self>) {
			changes.add(&self.snap, TimeUnit::Secs);
		}
		fn advance(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.snap.advance(time);
		}
	}
	
	impl FluxValue for Snap {
		type Value = i64;
		type Linear = f64;
		type Degree = Deg<0>;
		fn value(&self) -> Self::Linear {
			self.value
		}
		fn set_value(&mut self, value: Self::Linear) {
			self.value = value;
		}
		fn change(&self, _changes: &mut Changes<Self>) {}
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
				-1.469166666666667e-9,
				-1.4045833333333335e-18,
				0.06416666666666669e-27,
				-0.0004166666666666668e-36
			])
		);
		pos.advance(20*Secs);
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