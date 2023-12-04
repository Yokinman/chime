//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Mul, Sub};

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
	
	// ??? Could instead make the polynomial evaluable at `-infinity`.
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
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self;
}

impl<K: FluxKind> FluxAccum<'_, K> for () {
	fn from_kind(_kind: FluxAccumKind<'_, K>) -> Self {}
}

/// General accumulator arguments.
#[non_exhaustive]
pub enum FluxAccumKind<'a, K: FluxKind> {
	Value {
		value: &'a mut K::Value,
		depth: usize,
		time: Time,
		base_time: Time,
	},
	Poly {
		poly: &'a mut K,
		depth: usize,
		time: Time,
		base_time: Time,
	},
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
	pub fn real_roots(self) -> impl IntoIterator<Item=f64>
	where
		K: Roots
	{
		self.roots().into_iter()
			.filter(|r| !r.is_nan())
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	fn when_sign(&self, order: Ordering) -> TimeRanges
	where
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let initial_order = self.initial_order();
		
		 // Convert Roots to Ranges:
		let range_list = {
			let roots = self.real_roots().into_iter();
			
			 // Sort:
			let (capacity, _) = roots.size_hint();
			let mut roots = roots.fold(
				Vec::with_capacity(capacity),
				|mut list, r| {
					list.insert(list.partition_point(|&x| r > x), r);
					list
				}
			).into_iter();
			
			let (min_size, max_size) = roots.size_hint();
			let mut list = Vec::with_capacity(
				if order == Ordering::Equal {
					max_size.unwrap_or(min_size)
				} else {
					1 + (max_size.unwrap_or(min_size) / 2)
				}
			);
			let mut prev_point = if initial_order == Some(order) {
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
				} else if order == Ordering::Equal {
					list.push((point, point));
				} else {
					prev_point = Some(point);
				}
			}
			if let Some(point) = prev_point {
				list.push((point, f64::INFINITY));
			}
			list.into_iter()
		};
		
		 // Convert Ranges to Times:
		let offset = if order == Ordering::Equal {
			Time::ZERO
		} else {
			Time::new(0, 1)
		};
		range_list
			.filter_map(|(a, b)| {
				if let (
					Ok(a) | Err(a @ Time::ZERO),
					Ok(b) | Err(b @ Time::MAX)
				) = (
					time_try_from_secs(a, self.time)
						.map(|t| t.checked_add(offset).unwrap_or(Time::MAX)),
					time_try_from_secs(b, self.time)
						.map(|t| t.checked_sub(offset).unwrap_or(Time::ZERO))
				) {
					if a <= b {
						Some((a, b))
					} else {
						None
					}
				} else {
					None
				}
			})
			.fold(Vec::new(), |mut ranges, (a, b)| {
				 // Combine Adjacent Ranges:
				match ranges.last_mut() {
					Some((_, prev)) if *prev == a => {
						*prev = b;
					},
					_ => ranges.push((a, b))
				}
				ranges
			})
			.into_iter()
			.collect()
	}
	
	/// Times when the value is equal to zero.
	fn when_zero(&self) -> Times
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		self.real_roots()
			.into_iter()
			.filter_map(|t| time_try_from_secs(t, self.time).ok())
			.collect()
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
	type Output: IntoIterator<Item=f64>; // !!! Item=Self::Value
	fn roots(self) -> <Self as Roots>::Output;
}

/// `Duration::try_from_secs_f64`, but rounds toward the `basis`.
fn time_try_from_secs(mut t: f64, basis: Time) -> Result<Time, Time> {
	let sign = t.signum();
	t *= sign;
	
	const MANT_MASK: u64 = (1 << 52) - 1;
	const EXP_MASK: u64 = (1 << 11) - 1;
	
	let bits = t.to_bits();
	let mant = (bits & MANT_MASK) | (MANT_MASK + 1);
	let exp = (((bits >> 52) & EXP_MASK) as i16) - 1023;
	
	let time = if exp < -30 {
		// Too small; `1ns < 2s^-30`.
		Time::ZERO
	} else if exp < 0 {
		// No integer part.
		let nanos = ((((mant as u128) << (44 + exp)) * 1_000_000_000) >> (44 + 52)) as u32;
		Time::new(0, nanos)
	} else if exp < 52 {
		let secs = mant >> (52 - exp);
		let nanos = (((((mant << exp) & MANT_MASK) as u128) * 1_000_000_000) >> 52) as u32;
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

impl<A: FluxKind, const SIZE: usize, B: FluxKind, D: FluxKind> WhenDis<SIZE, B, D> for [Poly<A>; SIZE]
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: ops::Sqr,
	<<A as ops::Sub<B>>::Output as ops::Sqr>::Output:
		Add<Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ ops::Sub<
			<D as ops::Sqr>::Output,
			Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	A::Value: PartialOrd,
	D: FluxKind<Value=A::Value> + ops::Sqr,
{
	fn when_dis(&self, poly: &[Poly<B>; SIZE], order: Ordering, dis: &Poly<D>) -> TimeRanges {
		use ops::*;
		
		let time = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		
		for i in 0..SIZE {
			let x = (*self[i]).sub(*poly[i]);
			sum = sum + x.sqr();
		}
		
		Poly::new(sum.sub((**dis).sqr()), time)
			.when_sign(order)
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D>: private::PolyValue<SIZE> {
	fn when_dis_eq(&self, poly: &[Poly<B>; SIZE], dis: &Poly<D>) -> Times;
}

impl<A: FluxKind, const SIZE: usize, B: FluxKind, D: FluxKind> WhenDisEq<SIZE, B, D> for [Poly<A>; SIZE]
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: ops::Sqr,
	<<A as ops::Sub<B>>::Output as ops::Sqr>::Output:
		Add<Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ ops::Sub<
			<D as ops::Sqr>::Output,
			Output = <<A as ops::Sub<B>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	A::Value: PartialEq,
	D: FluxKind<Value=A::Value> + ops::Sqr,
{
	fn when_dis_eq(&self, poly: &[Poly<B>; SIZE], dis: &Poly<D>) -> Times {
		use ops::*;
		
		let time = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		
		for i in 0..SIZE {
			let x = (*self[i]).sub(*poly[i]);
			sum = sum + x.sqr();
		}
		
		Poly::new(sum.sub((**dis).sqr()), time)
			.when_zero()
	}
}