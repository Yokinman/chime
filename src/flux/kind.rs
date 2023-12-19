//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Mul, Sub};

use crate::Constant;
use crate::linear::{Linear, Scalar};
use crate::time::{Time, Times, TimeRanges};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind:
	'static + Copy + Clone + Debug
	+ From<<Self as FluxKind>::Value>
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Accum<'a>: FluxAccum<'a, Self>;
	
	fn value(&self) -> Self::Value;
	
	fn at(&self, time: Scalar) -> Self::Value;
	
	/// The order at or immediately preceding the value at time=0.
	/// 
	/// This should be the first non-zero [`FluxKind::value`] of this kind or
	/// its derivatives; reversed for odd derivatives.
	fn initial_order(&self) -> Option<Ordering>
	where
		Self::Value: PartialOrd;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

/// Combining [`FluxKind`] types.
/// 
/// Primarily this serves as a way to put two kinds of change-over-time into
/// the same space, for combination or comparison purposes.
pub mod ops {
	use std::ops;
	use super::FluxKind;
	use crate::linear::Scalar;
	
	/// Adding two kinds of change.
	pub trait Add<K: FluxKind = Self>: FluxKind<Value=K::Value> {
		type Output: FluxKind<Value=K::Value>;
		fn add(self, kind: K) -> <Self as Add<K>>::Output;
	}
	
	impl<A, B> Add<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value=A::Value>,
		<A as ops::Add<B>>::Output: FluxKind<Value=A::Value>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn add(self, kind: B) -> <A as ops::Add<B>>::Output {
			self + kind
		}
	}
	
	/// Differentiating two kinds of change.
	pub trait Sub<K: FluxKind = Self>: FluxKind<Value=K::Value> {
		type Output: FluxKind<Value=K::Value>;
		fn sub(self, kind: K) -> <Self as Sub<K>>::Output;
	}
	
	impl<A, B> Sub<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value=A::Value>,
		<A as ops::Add<B>>::Output: FluxKind<Value=A::Value>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn sub(self, kind: B) -> <A as ops::Add<B>>::Output {
			self + (kind * Scalar(-1.))
		}
	}
	
	/// Squaring a kind of change.
	pub trait Sqr: FluxKind {
		type Output: FluxKind<Value=Self::Value>;
		fn sqr(self) -> <Self as Sqr>::Output;
	}
	
	impl<K: FluxKind> Sqr for K
	where
		K: ops::Mul,
		<K as ops::Mul>::Output: FluxKind<Value=K::Value>,
	{
		type Output = <K as ops::Mul>::Output;
		fn sqr(self) -> <Self as Sqr>::Output {
			self * self
		}
	}
}

/// Change accumulator.
/// 
/// Converts a discrete pattern of change into a desired form.
pub trait FluxAccum<'a, K: FluxKind> {
	fn from_kind(kind: &'a mut K, depth: usize, time: Time) -> Self;
}

impl<K: FluxKind> FluxAccum<'_, K> for () {
	fn from_kind(_: &'_ mut K, _: usize, _: Time) -> Self {}
}

/// A polynomial wrapper for [`FluxKind`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<K> {
	inner: K,
	time: Time,
}

impl<K: FluxKind> Poly<K> {
	pub fn new(inner: K, time: Time) -> Self {
		Self {
			inner,
			time,
		}
	}
	
	pub fn with_value(value: K::Value) -> Self {
		Self {
			inner: K::from(value),
			time: Time::ZERO,
		}
	}
	
	pub fn with_time(time: Time) -> Self {
		Self {
			inner: K::zero(),
			time,
		}
	}
	
	pub fn time(&self) -> Time {
		self.time
	}
	
	pub fn at(&self, time: Time) -> K::Value {
		self.inner.at(Scalar(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
	
	pub fn sqr(self) -> Poly<<K as ops::Sqr>::Output>
	where
		K: ops::Sqr
	{
		Poly {
			inner: (*self).sqr(),
			time: self.time,
		}
	}
	
	/// All real-valued roots of this polynomial.
	pub fn real_roots(self) -> impl Iterator<Item=f64> + Send + Sync + Clone
	where
		K: Roots
	{
		self.roots().into_iter()
			.filter(|r| r.is_finite())
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	fn when_sign(&self, order: Ordering) -> TimeRanges
	where
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let basis = self.time;
		if let Some(initial_order) = self.initial_order() {
			TimeRanges::new(
				self.real_roots()
					.filter_map(move |r| time_try_from_secs(r, basis).ok()),
				basis,
				order,
				initial_order,
			)
		} else {
			TimeRanges::new(
				std::iter::empty(),
				basis,
				Ordering::Equal,
				Ordering::Equal,
			)
		}
	}
	
	/// Times when the value is equal to zero.
	fn when_zero(&self) -> Times
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		let basis = self.time;
		Times::new(self.real_roots()
			.filter_map(move |r| time_try_from_secs(r, basis).ok()))
	}
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Self {
			inner: K::zero(),
			time: Time::ZERO,
		}
	}
}

impl<K: FluxKind> From<K> for Poly<K> {
	fn from(value: K) -> Self {
		Self {
			inner: value,
			time: Time::ZERO,
		}
	}
}

impl<K> Deref for Poly<K> {
	type Target = K;
	fn deref(&self) -> &Self::Target {
		&self.inner
	}
}

impl<K> DerefMut for Poly<K> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.inner
	}
}

impl<A: FluxKind, B: FluxKind> Add<Poly<B>> for Poly<A>
where
	A: ops::Add<B>
{
	type Output = Poly<<A as ops::Add<B>>::Output>;
	fn add(self, rhs: Poly<B>) -> Self::Output {
		Poly {
			inner: (*self).add(*rhs),
			time: self.time,
		}
	}
}

impl<A: FluxKind, B: FluxKind> Sub<Poly<B>> for Poly<A>
where
	A: ops::Sub<B>
{
	type Output = Poly<<A as ops::Sub<B>>::Output>;
	fn sub(self, rhs: Poly<B>) -> Self::Output {
		Poly {
			inner: (*self).sub(*rhs),
			time: self.time,
		}
	}
}

impl<K: FluxKind> Mul<Scalar> for Poly<K> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		*self = *self * rhs;
		self
	}
}

/// Roots of a [`Poly`]nomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: FluxKind {
	type Output: IntoIterator<Item=f64, IntoIter = <Self as Roots>::IntoIter>;
	type IntoIter: Iterator<Item=f64> + Send + Sync + Clone;
	// !!! Item=Self::Value (requires Self::Value to implement some kind of
	// time conversion method)
	fn roots(self) -> <Self as Roots>::Output;
}

/// `Duration::try_from_secs_f64`, but always rounds down.
fn time_try_from_secs(mut t: f64, basis: Time) -> Result<Time, Time> {
	let sign = t.signum();
	t *= sign;
	
	const MANT_MASK: u64 = (1 << 52) - 1;
	const EXP_MASK: u64 = (1 << 11) - 1;
	const REM_MASK: u128 = (1 << 44) - 1;
	
	let bits = t.to_bits();
	let mant = (bits & MANT_MASK) | (MANT_MASK + 1);
	let exp = (((bits >> 52) & EXP_MASK) as i16) - 1023;
	
	let time = if exp < -30 {
		// Too small; `1ns < 2s^-30`.
		if sign <= 0. {
			Time::new(0, 1)
		} else {
			Time::ZERO
		}
	} else if exp < 0 {
		// No integer part.
		let nanos_tmp = ((mant as u128) << (44 + exp)) * 1_000_000_000;
		let mut nanos = (nanos_tmp >> (44 + 52)) as u32;
		if sign <= 0. && nanos_tmp & REM_MASK != 0 {
			nanos += 1;
		}
		Time::new(0, nanos)
	} else if exp < 52 {
		let secs = mant >> (52 - exp);
		let nanos_tmp = (((mant << exp) & MANT_MASK) as u128) * 1_000_000_000;
		let mut nanos = (nanos_tmp >> 52) as u32;
		if sign <= 0. && nanos_tmp & REM_MASK != 0 {
			nanos += 1;
		}
		Time::new(secs, nanos)
	} else if exp < 64 {
		// No fractional part.
		Time::from_secs(mant << (exp - 52))
	} else {
		// Too big.
		return if sign < 0. {
			Err(Time::ZERO)
		} else {
			Err(Time::MAX)
		}
	};
	
	if sign < 0. {
		basis.checked_sub(time).ok_or(Time::ZERO)
	} else {
		basis.checked_add(time).ok_or(Time::MAX)
	}
}

mod private {
	/// Sealed trait, only applied to [`super::Poly`] types.
	pub trait PolyValue<const SIZE: usize> {}
}

impl<K> private::PolyValue<1> for Poly<K> {}
impl<K, const SIZE: usize> private::PolyValue<SIZE> for [Poly<K>; SIZE] {}

/// [`crate::Flux::when`] predictive comparison.
pub trait When<B>: private::PolyValue<1> {
	fn when(self, order: Ordering, poly: Poly<B>) -> TimeRanges;
}

impl<A: FluxKind, B: FluxKind> When<B> for Poly<A>
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: Roots + PartialOrd,
	A::Value: PartialOrd,
{
	fn when(self, order: Ordering, poly: Poly<B>) -> TimeRanges {
		(self - poly).when_sign(order)
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B>: private::PolyValue<1> {
	fn when_eq(self, poly: Poly<B>) -> Times;
}

impl<A: FluxKind, B: FluxKind> WhenEq<B> for Poly<A>
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	fn when_eq(self, poly: Poly<B>) -> Times {
		(self - poly).when_zero()
	}
}

/// [`crate::FluxVec::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B, D>: private::PolyValue<SIZE> {
	fn when_dis(&self, poly: &[Poly<B>; SIZE], order: Ordering, dis: &Poly<D>) -> TimeRanges;
}

impl<A: FluxKind, const SIZE: usize, B: FluxKind, D> WhenDis<SIZE, B, D> for [Poly<A>; SIZE]
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: ops::Sqr,
	<<A as ops::Sub<B>>::Output as ops::Sqr>::Output:
		Add<Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	A::Value: Mul<Output=A::Value> + PartialOrd,
	D: FluxKind<Value=A::Value> + ops::Sqr + Add<Constant<D::Value>, Output=D>,
{
	fn when_dis(&self, poly: &[Poly<B>; SIZE], order: Ordering, dis: &Poly<D>) -> TimeRanges {
		use ops::*;
		
		let time = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		let mut abs_dis = A::Value::zero();
		
		for i in 0..SIZE {
			sum = sum + (*self[i]).sub(*poly[i]).sqr();
			
			 // Roughly the distance from zero at the next point of collision:
			let x = (*self[i]).at(Scalar(2.)) + (*poly[i]).at(Scalar(2.));
			abs_dis = abs_dis + x*x;
		}
		abs_dis = Mul::<Scalar>::mul(abs_dis.sqrt(), Scalar(0.5));
		
		let curr_dis = abs_dis + sum.at(Scalar(0.)).sqrt();
		let goal_dis = abs_dis + (**dis).at(Scalar(0.));
		
		 // Buffer Distance (undershoot prediction to avoid rounding):
		let epsilon = match curr_dis.partial_cmp(&goal_dis) {
			Some(Ordering::Greater) => goal_dis.next_up()   + Mul::<Scalar>::mul(goal_dis, Scalar(-1.)),
			Some(Ordering::Less)    => goal_dis.next_down() + Mul::<Scalar>::mul(goal_dis, Scalar(-1.)),
			_ => A::Value::zero()
		};
		
		let mut dis_poly = sum.sub((**dis + Constant::from(epsilon)).sqr());
		
		 // Guaranteed Start at Zero:
		if curr_dis == goal_dis + epsilon {
			dis_poly = dis_poly.sub(<D as Sqr>::Output::from(dis_poly.at(Scalar(0.))));
		}
		
		Poly::new(dis_poly, time)
			.when_sign(order)
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D>: private::PolyValue<SIZE> {
	fn when_dis_eq(&self, poly: &[Poly<B>; SIZE], dis: &Poly<D>) -> Times;
}

impl<A: FluxKind, const SIZE: usize, B: FluxKind, D> WhenDisEq<SIZE, B, D> for [Poly<A>; SIZE]
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: ops::Sqr,
	<<A as ops::Sub<B>>::Output as ops::Sqr>::Output:
		Add<Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	A::Value: Mul<Output=A::Value> + PartialOrd,
	D: FluxKind<Value=A::Value> + ops::Sqr + Add<Constant<D::Value>, Output=D>,
{
	fn when_dis_eq(&self, poly: &[Poly<B>; SIZE], dis: &Poly<D>) -> Times {
		use ops::*;
		
		let time = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		let mut abs_dis = A::Value::zero();
		
		for i in 0..SIZE {
			sum = sum + (*self[i]).sub(*poly[i]).sqr();
			
			 // Roughly the distance from zero at the next point of collision:
			let x = (*self[i]).at(Scalar(2.)) + (*poly[i]).at(Scalar(2.));
			abs_dis = abs_dis + x*x;
		}
		abs_dis = Mul::<Scalar>::mul(abs_dis.sqrt(), Scalar(0.5));
		
		let curr_dis = abs_dis + sum.at(Scalar(0.)).sqrt();
		let goal_dis = abs_dis + (**dis).at(Scalar(0.));
		
		 // Buffer Distance (undershoot prediction to avoid rounding):
		let epsilon = match curr_dis.partial_cmp(&goal_dis) {
			Some(Ordering::Greater) => goal_dis.next_up()   + Mul::<Scalar>::mul(goal_dis, Scalar(-1.)),
			Some(Ordering::Less)    => goal_dis.next_down() + Mul::<Scalar>::mul(goal_dis, Scalar(-1.)),
			_ => A::Value::zero()
		};
		
		let mut dis_poly = sum.sub((**dis + Constant::from(epsilon)).sqr());
		
		 // Guaranteed Start at Zero:
		if curr_dis == goal_dis + epsilon {
			dis_poly = dis_poly.sub(<D as Sqr>::Output::from(dis_poly.at(Scalar(0.))));
		}
		
		Poly::new(dis_poly, time)
			.when_zero()
	}
}