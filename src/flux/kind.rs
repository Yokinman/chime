//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Mul, Sub};

use crate::linear::{Linear, LinearVec, Scalar};
use crate::time;
use crate::time::{Time, /*Times,*/ TimeRanges};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind:
	'static + Copy + Clone + Debug + Send + Sync
	+ From<<Self as FluxKind>::Value>
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Accum<'a>: FluxAccum<'a, Self>;
	
	type OutAccum<'a>: FluxAccum<'a, Self>;
	
	fn at(&self, time: Scalar) -> Self::Value;
	
	fn rate_at(&self, time: Scalar) -> Self::Value;
	
	fn to_time(self, time: Scalar) -> Self;
	
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
	
	fn value(&self) -> Self::Value {
		self.at(Scalar(0.))
	}
}

/// Multidimensional kind of change.
pub trait FluxKindVec<const SIZE: usize>: FluxKind
where
	Self::Value: LinearVec<SIZE>,
{
	type Kind: FluxKind<Value = <Self::Value as LinearVec<SIZE>>::Value>;
	fn index_kind(&self, index: usize) -> Self::Kind;
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
			// ??? This could pass through to `ops::Sub` directly, but then
			// stuff like [`sum::Sum`] would need a whole extra set of macro
			// implementations. For now, this will just reuse `ops::Add`.
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
	
	pub fn rate_at(&self, time: Time) -> K::Value {
		self.inner.rate_at(Scalar(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
	
	pub fn to_time(mut self, time: Time) -> Self {
		if self.time != time {
			self.inner = self.inner.to_time(Scalar(if time > self.time {
				(time - self.time).as_secs_f64()
			} else {
				-(self.time - time).as_secs_f64()
			}));
			self.time = time;
		}
		self
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
	fn when_sign(&self, order: Ordering, mut f: impl RootFilterMap + 'static) -> TimeRanges
	where
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let basis = self.time;
		let initial_order = self.initial_order().unwrap_or(Ordering::Equal);
		TimeRanges::new(
			self.real_roots()
				.filter_map(move |x| time_try_from_secs(x, basis).ok())
				.enumerate()
				.filter_map(move |(index, time)|
					f(time, (index % 2 == 1) == (initial_order == order))
				),
			basis,
			initial_order,
			order,
		)
	}
	
	/// Times when the value is equal to zero.
	fn when_zero(&self, f: impl RootFilterMap + 'static) -> TimeRanges
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		let basis = self.time;
		let times = time::Times::new(self.real_roots()
			.filter_map(move |x| time_try_from_secs(x, basis).ok()));
		
		#[derive(Clone)]
		struct RangeBuilder<F> {
			times: time::Times,
			next: Option<[Time; 2]>,
			f: F
		}
		impl<F: RootFilterMap> Iterator for RangeBuilder<F> {
			type Item = [Time; 2];
			fn next(&mut self) -> Option<Self::Item> {
				if let Some(time) = self.times.next() {
					let [a, mut b] = self.next.unwrap_or_else(|| [
						(self.f)(time, false).unwrap_or(Time::ZERO)
							.checked_sub(time::NANOSEC).unwrap_or(Time::ZERO),
						(self.f)(time, true).unwrap_or(Time::MAX)
							.checked_add(time::NANOSEC).unwrap_or(Time::MAX),
					]);
					debug_assert!(a <= b);
					while let Some(t) = self.times.peek() {
						debug_assert!(t > a);
						let x = (self.f)(t, false).unwrap_or(Time::ZERO)
							.checked_sub(time::NANOSEC).unwrap_or(Time::ZERO);
						let y = (self.f)(t, true).unwrap_or(Time::MAX)
							.checked_add(time::NANOSEC).unwrap_or(Time::MAX);
						debug_assert!(a <= x);
						if x <= b {
							if b < y {
								b = y;
							}
							self.times.next();
						} else {
							self.next = Some([x, y]);
							break
						}
					}
					Some([a, b])
				} else {
					None
				}
			}
			fn size_hint(&self) -> (usize, Option<usize>) {
				self.times.size_hint()
			}
		}
		
		TimeRanges::new(
			RangeBuilder { times, next: None, f }.flatten(),
			Time::ZERO,
			if self.is_zero() {
				Ordering::Less
			} else {
				Ordering::Greater
			},
			Ordering::Less,
		)
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
			inner: (*self).add(*rhs.to_time(self.time)),
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
			inner: (*self).sub(*rhs.to_time(self.time)),
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
		if sign == -1. && t != 0. {
			Time::new(0, 1)
		} else {
			Time::ZERO
		}
	} else if exp < 0 {
		// No integer part.
		let nanos_tmp = ((mant as u128) << (44 + exp)) * 1_000_000_000;
		let mut nanos = (nanos_tmp >> (44 + 52)) as u32;
		if sign == -1. && (nanos_tmp & REM_MASK) != 0 {
			nanos += 1;
		}
		Time::new(0, nanos)
	} else if exp < 52 {
		let secs = mant >> (52 - exp);
		let nanos_tmp = (((mant << exp) & MANT_MASK) as u128) * 1_000_000_000;
		let mut nanos = (nanos_tmp >> 52) as u32;
		if sign == -1. && (nanos_tmp & REM_MASK) != 0 {
			nanos += 1;
		}
		Time::new(secs, nanos)
	} else if exp < 64 {
		// No fractional part.
		Time::from_secs(mant << (exp - 52))
	} else {
		// Too big.
		return if sign == -1. {
			Err(Time::ZERO)
		} else {
			Err(Time::MAX)
		}
	};
	
	if sign == -1. {
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

/// Function that converts a root value to a Time, or ignores it.
trait RootFilterMap: FnMut(Time, bool) -> Option<Time> + Clone + Send + Sync {}
impl<T: FnMut(Time, bool) -> Option<Time> + Clone + Send + Sync> RootFilterMap for T {}

const ROOT_FILTER_TRIES: usize = 100; // Arbitrary number

fn root_filter_map<T: FluxKind>(
	a_poly: Poly<T>,
	b_poly: Poly<impl FluxKind<Value=T::Value>>,
	diff_poly: Poly<impl FluxKind<Value=T::Value>>,
) -> impl RootFilterMap
where
	T::Value: PartialEq
{
	move |mut time, is_end| {
		// Covers the local range and the range immediately preceding it, but
		// stops where the trend reverses. Careful, this logic is precise.
		let mut next_time = if is_end {
			time.checked_add(time::NANOSEC)?
		} else {
			time.checked_sub(time::NANOSEC)?
		};
		let sign = diff_poly.rate_at(time).sign();
		for i in 0..2 {
			if i != 0 {
				if is_end {
					break
				}
				time = next_time;
				next_time = if is_end {
					time.checked_add(time::NANOSEC)?
				} else {
					time.checked_sub(time::NANOSEC)?
				};
			}
			let diff = if is_end {
				T::Value::zero()
			} else {
				a_poly.at(time) - b_poly.at(time)
			};
			for _ in 0..ROOT_FILTER_TRIES {
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					return Some(time)
				}
				if diff != a_poly.at(next_time) - b_poly.at(next_time) {
					break
				}
				time = next_time;
				next_time = if is_end {
					time.checked_add(time::NANOSEC)?
				} else {
					time.checked_sub(time::NANOSEC)?
				};
			}
		}
		Some(time)
	}
}

fn dis_root_filter_map<T: FluxKind, const SIZE: usize>(
	pos_poly: Poly<T>,
	dis_poly: Poly<impl FluxKind<Value=T::Value>>,
	diff_poly: Poly<impl FluxKind<Value=T::Value>>,
	a_pos: [Poly<impl FluxKind<Value=T::Value>>; SIZE],
	b_pos: [Poly<impl FluxKind<Value=T::Value>>; SIZE],
) -> impl RootFilterMap
where
	T::Value: PartialEq + Mul<Output=T::Value>,
{
	move |mut time, is_end| {
		// Covers the local range, but stops where the trend reverses and
		// undershoots to avoid rounding past. Careful, this logic is precise.
		
		// !!! Tricky issue - if the peak of a value occurs near 0, and its rate
		// of change at 0 is briefly moving away from the target value, only to
		// immediately start towards the target, an event may be delayed by a
		// full nanosecond. This compounds with the issue of rounding, leading
		// to a value moving towards its target more than it should. This can
		// occur repeatedly in sequence, letting a value round past its target.
		// - Including the point at 0 would be awkward, since the rate of change
		//   is moving away. Events shouldn't have to expect this.
		// - Rewriting all system internals to use a custom Time type, which
		//   pairs `std::time::Duration` with a floating-point type to measure
		//   sub-nanosecond offsets would work. One issue is that events would
		//   have to opt-in, since Bevy's `Time::elapsed` method uses Duration.
		// - Could refactor system internals to make it so that when an event
		//   occurs, its prediction values are automatically shifted around to
		//   ensure that the rounding doesn't happen. I don't really know if
		//   this would work, and I think it would cause its own issues.
		
		let mut next_time = if is_end {
			time.checked_add(time::NANOSEC)?
		} else {
			time.checked_sub(time::NANOSEC)?
		};
		let sign = diff_poly.rate_at(time).sign();
		
		 // Rounding Buffer:
		if !is_end {
			let round_factor = Scalar(0.5 / (SIZE as f64).sqrt());
			for _ in 0..ROOT_FILTER_TRIES {
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					return Some(time)
				}
				let dis = dis_poly.at(next_time);
				let mut a_dis = T::Value::zero();
				let mut b_dis = T::Value::zero();
				let mut a_diff = T::Value::zero();
				for i in 0..SIZE {
					let a = a_pos[i].at(next_time);
					let b = b_pos[i].at(next_time);
					a_dis = a_dis + a*a;
					b_dis = b_dis + b*b;
					let x = a - b;
					a_diff = a_diff + x*x;
				}
				a_dis = a_dis.sqrt();
				b_dis = b_dis.sqrt();
				a_diff = Mul::<Scalar>::mul(a_diff.sqrt() - dis, round_factor);
				// !!! This could probably be refined, but it works for now:
				if a_dis != a_dis + a_diff && b_dis != b_dis + a_diff && dis != dis + a_diff {
					let b_diff = Mul::<Scalar>::mul(pos_poly.at(next_time).sqrt() - dis, round_factor);
					if a_dis != a_dis + b_diff && b_dis != b_dis + b_diff && dis != dis + b_diff {
						break
					}
				}
				time = next_time;
				next_time = if is_end {
					time.checked_add(time::NANOSEC)?
				} else {
					time.checked_sub(time::NANOSEC)?
				};
			}
			time = next_time;
			next_time = if is_end {
				time.checked_add(time::NANOSEC)?
			} else {
				time.checked_sub(time::NANOSEC)?
			};
		}
		
		 // Fully Cover Local Range:
		let mut diff = T::Value::zero();
		if !is_end {
			for i in 0..SIZE {
				let x = a_pos[i].at(time) - b_pos[i].at(time);
				diff = diff + x*x;
			}
			let dis = dis_poly.at(time);
			diff = diff - dis*dis;
		}
		for _ in 0..ROOT_FILTER_TRIES {
			let rate = diff_poly.rate_at(next_time);
			if sign != rate.sign() && !rate.is_zero() {
				return Some(time)
			}
			let mut pos = T::Value::zero();
			for i in 0..SIZE {
				let x = a_pos[i].at(next_time) - b_pos[i].at(next_time);
				pos = pos + x*x;
			}
			let dis = dis_poly.at(next_time);
			if diff != pos - dis*dis {
				break
			}
			time = next_time;
			next_time = if is_end {
				time.checked_add(time::NANOSEC)?
			} else {
				time.checked_sub(time::NANOSEC)?
			};
		}
		
		Some(time)
	}
}

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
		let diff_poly = self - poly;
		diff_poly
			.when_sign(order, root_filter_map(self, poly, diff_poly))
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B>: private::PolyValue<1> {
	fn when_eq(self, poly: Poly<B>) -> TimeRanges;
}

impl<A: FluxKind, B: FluxKind> WhenEq<B> for Poly<A>
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	fn when_eq(self, poly: Poly<B>) -> TimeRanges {
		let diff_poly = self - poly;
		diff_poly
			.when_zero(root_filter_map(self, poly, diff_poly))
	}
}

/// [`crate::FluxVec::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B, D>: private::PolyValue<SIZE> {
	fn when_dis(self, poly: [Poly<B>; SIZE], order: Ordering, dis: Poly<D>) -> TimeRanges;
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
	D: FluxKind<Value=A::Value> + ops::Sqr,
{
	fn when_dis(self, poly: [Poly<B>; SIZE], order: Ordering, dis: Poly<D>) -> TimeRanges {
		use ops::*;
		
		let basis = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + (*self[i].to_time(basis))
				.sub(*poly[i].to_time(basis))
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_sign(order, dis_root_filter_map(sum, dis, diff_poly, self, poly))
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D>: private::PolyValue<SIZE> {
	fn when_dis_eq(self, poly: [Poly<B>; SIZE], dis: Poly<D>) -> TimeRanges;
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
		+ PartialEq,
	A::Value: Mul<Output=A::Value> + PartialEq,
	D: FluxKind<Value=A::Value> + ops::Sqr,
{
	fn when_dis_eq(self, poly: [Poly<B>; SIZE], dis: Poly<D>) -> TimeRanges {
		use ops::*;
		
		let basis = if SIZE == 0 {
			Time::ZERO
		} else {
			self[0].time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + (*self[i].to_time(basis))
				.sub(*poly[i].to_time(basis))
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_zero(dis_root_filter_map(sum, dis, diff_poly, self, poly))
	}
}