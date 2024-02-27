//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Mul, Sub};

use crate::linear::{Linear, LinearVec, Scalar};
use crate::time;
use crate::time::{Time, TimeIter, TimeRanges};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind:
	'static + Copy + Clone + Debug + Send + Sync
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Accum<'a>;
	
	type OutAccum<'a>;
	
	fn from_value(value: Self::Value) -> Self;
	
	fn as_accum(&mut self, depth: usize, time: Time) -> Self::Accum<'_>;
	
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
pub trait FluxKindVec<const SIZE: usize>: Clone + Send + Sync + 'static {
	type Kind: FluxKind;
	fn index_kind(&self, index: usize) -> Self::Kind;
}

impl<const SIZE: usize, T: FluxKind> FluxKindVec<SIZE> for [T; SIZE] {
	type Kind = T;
	fn index_kind(&self, index: usize) -> Self::Kind {
		self[index]
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
		Self::new(K::from_value(value), Time::ZERO)
	}
	
	pub fn with_time(time: Time) -> Self {
		Self::new(K::zero(), time)
	}
	
	pub fn into_inner(self) -> K {
		self.inner
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
		Poly::new(self.inner.sqr(), self.time)
	}
	
	/// All real-valued roots of this polynomial.
	pub fn real_roots(&self) -> impl Iterator<Item=f64> + Send + Sync + Clone
	where
		K: Roots
	{
		self.inner.roots().into_iter()
			.filter(|r| r.is_finite())
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	fn when_sign(&self, order: Ordering, f: impl RootFilterMap + 'static) -> TimeRanges<impl TimeIter>
	where
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let basis = self.time;
		let basis_order = self.inner.initial_order().unwrap_or(Ordering::Equal);
		let times = self.real_roots()
			.filter_map(move |x| time_try_from_secs(x, basis).ok());
		TimeRanges::new(times, basis, basis_order, order)
			.into_filtered(f)
	}
	
	/// Times when the value is equal to zero.
	fn when_zero(&self, f: impl RootFilterMap + 'static) -> TimeRanges<impl TimeIter>
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		let basis = self.time;
		let basis_order = if self.inner.is_zero() {
			Ordering::Equal
		} else {
			Ordering::Greater
		};
		let times = self.real_roots()
			.filter_map(move |x| time_try_from_secs(x, basis).ok());
		TimeRanges::new(times, basis, basis_order, Ordering::Equal)
			.into_filtered(f)
	}
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Self::new(K::zero(), Time::ZERO)
	}
}

impl<K: FluxKind> From<K> for Poly<K> {
	fn from(value: K) -> Self {
		Self::new(value, Time::ZERO)
	}
}

impl<A: FluxKind, B: FluxKind> Add<Poly<B>> for Poly<A>
where
	A: ops::Add<B>
{
	type Output = Poly<<A as ops::Add<B>>::Output>;
	fn add(self, rhs: Poly<B>) -> Self::Output {
		Poly::new(
			self.inner.add(rhs.to_time(self.time).inner),
			self.time,
		)
	}
}

impl<A: FluxKind, B: FluxKind> Sub<Poly<B>> for Poly<A>
where
	A: ops::Sub<B>
{
	type Output = Poly<<A as ops::Sub<B>>::Output>;
	fn sub(self, rhs: Poly<B>) -> Self::Output {
		Poly::new(
			self.inner.sub(rhs.to_time(self.time).inner),
			self.time,
		)
	}
}

impl<K: FluxKind> Mul<Scalar> for Poly<K> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.inner = self.inner * rhs;
		self
	}
}

/// Multidimensional polynomial.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PolyVec<K, const SIZE: usize> {
	inner: K,
	time: Time,
}

impl<K, const SIZE: usize> PolyVec<K, SIZE> {
	pub fn into_inner(self) -> K {
		self.inner
	}
	
	pub fn time(&self) -> Time {
		self.time
	}
}

impl<K: FluxKind, const SIZE: usize> PolyVec<K, SIZE> {	
	pub fn at(&self, time: Time) -> K::Value {
		self.inner.at(Scalar(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
}

impl<K: FluxKindVec<SIZE>, const SIZE: usize> PolyVec<K, SIZE> {
	pub fn new(inner: K, time: Time) -> Self {
		Self {
			inner,
			time,
		}
	}
	
	pub fn index_poly(&self, index: usize) -> Poly<K::Kind> {
		Poly::new(self.inner.index_kind(index), self.time)
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

/// Function that converts a root value to a Time, or ignores it.
pub(crate) trait RootFilterMap: FnMut(Time, bool) -> Option<Time> + Clone + Send + Sync {}
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

fn dis_root_filter_map<T: FluxKind, const SIZE: usize, A, B>(
	pos_poly: Poly<T>,
	dis_poly: Poly<impl FluxKind<Value=T::Value>>,
	diff_poly: Poly<impl FluxKind<Value=T::Value>>,
	a_pos: PolyVec<A, SIZE>,
	b_pos: PolyVec<B, SIZE>,
) -> impl RootFilterMap
where
	T::Value: PartialEq + Mul<Output=T::Value>,
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	A::Kind: FluxKind<Value=T::Value>,
	B::Kind: FluxKind<Value=T::Value>,
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
				
				 // Calculate Accurate Distances:
				let dis = dis_poly.at(next_time);
				let mut a_dis = T::Value::zero();
				let mut b_dis = T::Value::zero();
				let mut a_diff = T::Value::zero();
				for i in 0..SIZE {
					let a = a_pos.index_poly(i).at(next_time);
					let b = b_pos.index_poly(i).at(next_time);
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
				let x = a_pos.index_poly(i).at(time) - b_pos.index_poly(i).at(time);
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
			
			 // Calculate Accurate Distance:
			let mut pos = T::Value::zero();
			for i in 0..SIZE {
				let x = a_pos.index_poly(i).at(next_time) - b_pos.index_poly(i).at(next_time);
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
pub trait When<B> {
	fn when(self, order: Ordering, poly: Poly<B>) -> TimeRanges<impl TimeIter>;
}

impl<A: FluxKind, B: FluxKind> When<B> for Poly<A>
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: Roots + PartialOrd,
	A::Value: PartialOrd,
{
	fn when(self, order: Ordering, poly: Poly<B>) -> TimeRanges<impl TimeIter> {
		let diff_poly = self - poly;
		diff_poly
			.when_sign(order, root_filter_map(self, poly, diff_poly))
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B> {
	fn when_eq(self, poly: Poly<B>) -> TimeRanges<impl TimeIter>;
}

impl<A: FluxKind, B: FluxKind> WhenEq<B> for Poly<A>
where
	A: ops::Sub<B>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	fn when_eq(self, poly: Poly<B>) -> TimeRanges<impl TimeIter> {
		let diff_poly = self - poly;
		diff_poly
			.when_zero(root_filter_map(self, poly, diff_poly))
	}
}

/// [`crate::FluxVec::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B: FluxKind, D> {
	fn when_dis(
		self,
		poly: PolyVec<impl FluxKindVec<SIZE, Kind=B>, SIZE>,
		order: Ordering,
		dis: Poly<D>
	) -> TimeRanges<impl TimeIter>;
}

impl<T, A: FluxKind, const SIZE: usize, B: FluxKind, D> WhenDis<SIZE, B, D> for PolyVec<T, SIZE>
where
	T: FluxKindVec<SIZE, Kind=A>,
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
	fn when_dis(
		self,
		poly: PolyVec<impl FluxKindVec<SIZE, Kind=B>, SIZE>,
		order: Ordering,
		dis: Poly<D>
	) -> TimeRanges<impl TimeIter> {
		use ops::*;
		
		let basis = if SIZE == 0 {
			Time::ZERO
		} else {
			self.index_poly(0).time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).to_time(basis).inner
				.sub(poly.index_poly(i).to_time(basis).inner)
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_sign(order, dis_root_filter_map(sum, dis, diff_poly, self, poly))
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B: FluxKind, D> {
	fn when_dis_eq(
		self,
		poly: PolyVec<impl FluxKindVec<SIZE, Kind=B>, SIZE>,
		dis: Poly<D>
	) -> TimeRanges<impl TimeIter>;
}

impl<T, A: FluxKind, const SIZE: usize, B: FluxKind, D> WhenDisEq<SIZE, B, D> for PolyVec<T, SIZE>
where
	T: FluxKindVec<SIZE, Kind=A>,
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
	fn when_dis_eq(
		self,
		poly: PolyVec<impl FluxKindVec<SIZE, Kind=B>, SIZE>,
		dis: Poly<D>
	) -> TimeRanges<impl TimeIter> {
		use ops::*;
		
		let basis = if SIZE == 0 {
			Time::ZERO
		} else {
			self.index_poly(0).time
		};
		
		let mut sum = <<A as Sub<B>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).to_time(basis).inner
				.sub(poly.index_poly(i).to_time(basis).inner)
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_zero(dis_root_filter_map(sum, dis, diff_poly, self, poly))
	}
}

#[test]
fn consistent_sign_pred() {
	use crate::sum::Sum;
	fn toast(time: Time) -> Vec<(Time, Time)> {
		let poly = PolyVec::new([
			Sum::new(-2., [5., -2.]),
			Sum::zero()
		], time);
		poly
			.when_dis(
				PolyVec::new([Sum::<f64, 2>::zero(); 2], time),
				Ordering::Greater,
				Poly::new(crate::Constant::from(1.), time),
			)
			.map(|(a, b)| (
				a.saturating_sub(time),
				if b == Time::MAX {
					b
				} else {
					b.saturating_sub(time)
				})
			)
			.collect::<Vec<_>>()
	}
	for t in 0..6 {
		assert_eq!(
			toast(Time::from_secs_f64((t as f64) * 0.5)),
			toast(Time::from_secs_f64(((t+1) as f64) * 0.5))
		);
	}
}