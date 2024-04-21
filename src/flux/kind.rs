//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub};

use crate::linear::{Linear, LinearIso, LinearIsoVec, LinearVec, Scalar};
use crate::time::Time;

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
	
	/// The order at or immediately preceding the value at a time.
	/// 
	/// This should be the first non-zero [`FluxKind::at`] value of this kind
	/// or its derivatives; reversed for odd derivatives.
	fn initial_order(&self, time: Scalar) -> Option<Ordering>
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
	type Value: LinearVec<SIZE, Value = <Self::Kind as FluxKind>::Value>;
	fn index_kind(&self, index: usize) -> Self::Kind;
}

impl<const SIZE: usize, T: FluxKind> FluxKindVec<SIZE> for [T; SIZE] {
	type Kind = T;
	type Value = [T::Value; SIZE];
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
pub struct Poly<K, I> {
	inner: K,
	time: Time,
	iso: PhantomData<I>,
}

impl<K, I> Clone for Poly<K, I>
where
	K: Clone
{
	fn clone(&self) -> Self {
		Self {
			inner: self.inner.clone(),
			time: self.time,
			iso: self.iso,
		}
	}
}

impl<K, I> Copy for Poly<K, I>
where
	K: Copy
{}

impl<K, I> Debug for Poly<K, I>
where
	K: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		f.debug_tuple("Poly")
			.field(&self.inner)
			.field(&self.time)
			.finish()
	}
}

impl<K, I> PartialEq for Poly<K, I>
where
	K: PartialEq
{
	fn eq(&self, other: &Self) -> bool {
		self.inner == other.inner && self.time == other.time
	}
}

impl<K: FluxKind> Poly<K, K::Value> {
	pub fn new(inner: K, time: Time) -> Self {
		Self {
			inner,
			time,
			iso: PhantomData,
		}
	}
	
	pub fn with_value(value: K::Value) -> Self {
		Self::new(K::from_value(value), Time::ZERO)
	}
	
	pub fn with_time(time: Time) -> Self {
		Self::new(K::zero(), time)
	}
}

impl<K: FluxKind, I> Poly<K, I> {
	pub fn with_iso<T: LinearIso<K::Value>>(self) -> Poly<K, T> {
		Poly {
			inner: self.inner,
			time: self.time,
			iso: PhantomData,
		}
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
	
	pub fn initial_order(&self, time: Time) -> Option<Ordering>
	where
		K::Value: PartialOrd
	{
		self.inner.initial_order(Scalar(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
}

impl<K: FluxKind, I: LinearIso<K::Value>> Poly<K, I> {
	pub fn sqr(self) -> Poly<<K as ops::Sqr>::Output, I>
	where
		K: ops::Sqr
	{
		Poly::new(self.inner.sqr(), self.time)
			.with_iso()
	}
	
	/// All real-valued roots of this polynomial.
	#[allow(dead_code)]
	pub(crate) fn real_roots(&self) -> impl Iterator<Item=f64> + Send + Sync + Clone
	where
		K: Roots,
		<K as Roots>::Output: IntoIterator<Item=f64>,
		<<K as Roots>::Output as IntoIterator>::IntoIter: Send + Sync + Clone,
	{
		self.inner.roots().into_iter()
			.filter(|r| r.is_finite())
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	pub(crate) fn when_sign<F>(self, order: Ordering, filter: F) -> crate::pred::PredFilter<crate::pred::Pred<K, I>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let pred = crate::pred::Pred {
			poly: self,
			order
		};
		crate::pred::PredFilter { pred, filter }
	}
	
	/// Times when the value is equal to zero.
	pub(crate) fn when_zero<F>(self, filter: F) -> crate::pred::PredFilter<crate::pred::PredEq<K, I>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		crate::pred::PredFilter {
			pred: crate::pred::PredEq { poly: self },
			filter,
		}
	}
}

impl<K: FluxKind, I: LinearIso<K::Value>> Default for Poly<K, I> {
	fn default() -> Self {
		Poly::new(K::zero(), Time::ZERO)
			.with_iso()
	}
}

impl<K: FluxKind> From<K> for Poly<K, K::Value> {
	fn from(value: K) -> Self {
		Self::new(value, Time::ZERO)
	}
}

impl<A: FluxKind, B: FluxKind, I, J> Add<Poly<B, J>> for Poly<A, I>
where
	A: ops::Add<B>
{
	type Output = Poly<<A as ops::Add<B>>::Output, I>;
	fn add(self, rhs: Poly<B, J>) -> Self::Output {
		Poly {
			inner: self.inner.add(rhs.to_time( self.time).inner),
			time: self.time,
			iso: PhantomData,
		}
	}
}

impl<A: FluxKind, B: FluxKind, I, J> Sub<Poly<B, J>> for Poly<A, I>
where
	A: ops::Sub<B>
{
	type Output = Poly<<A as ops::Sub<B>>::Output, I>;
	fn sub(self, rhs: Poly<B, J>) -> Self::Output {
		Poly {
			inner: self.inner.sub(rhs.to_time(self.time).inner),
			time: self.time,
			iso: PhantomData,
		}
	}
}

impl<K, I> Mul<Scalar> for Poly<K, I>
where
	K: Mul<Scalar, Output=K>
{
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.inner = self.inner * rhs;
		self
	}
}

/// Multidimensional polynomial.
pub struct PolyVec<const SIZE: usize, K, I> {
	inner: K,
	time: Time,
	iso: PhantomData<I>,
}

impl<const SIZE: usize, K, I> Clone for PolyVec<SIZE, K, I>
where
	K: Clone
{
	fn clone(&self) -> Self {
		Self {
			inner: self.inner.clone(),
			time: self.time,
			iso: self.iso,
		}
	}
}

impl<const SIZE: usize, K, I> Copy for PolyVec<SIZE, K, I>
where
	K: Copy
{}

impl<const SIZE: usize, K, I> PolyVec<SIZE, K, I> {
	pub fn into_inner(self) -> K {
		self.inner
	}
	
	pub fn time(&self) -> Time {
		self.time
	}
}

impl<const SIZE: usize, K: FluxKind, I> PolyVec<SIZE, K, I> {
	pub fn at(&self, time: Time) -> K::Value {
		self.inner.at(Scalar(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
}

impl<const SIZE: usize, K: FluxKindVec<SIZE>> PolyVec<SIZE, K, K::Value> {
	pub fn new(inner: K, time: Time) -> Self {
		Self {
			inner,
			time,
			iso: PhantomData,
		}
	}
}

impl<const SIZE: usize, K: FluxKindVec<SIZE>, I> PolyVec<SIZE, K, I> {
	pub fn with_iso<T: LinearIsoVec<SIZE, K::Value>>(self) -> PolyVec<SIZE, K, T> {
		PolyVec {
			inner: self.inner,
			time: self.time,
			iso: PhantomData,
		}
	}
}

impl<const SIZE: usize, K: FluxKindVec<SIZE>, I: LinearIsoVec<SIZE, K::Value>> PolyVec<SIZE, K, I> {
	pub fn index_poly(&self, index: usize) -> Poly<K::Kind, I::Value> {
		Poly::new(self.inner.index_kind(index), self.time)
			.with_iso()
	}
}

/// Roots of a [`Poly`]nomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: FluxKind {
	type Output: IntoTimes;
	fn roots(self) -> <Self as Roots>::Output;
}

/// Conversion from some [`Roots::Output`] into an iterator of time.
pub trait IntoTimes {
	type TimeIter: Iterator<Item=LinearTime>;
	fn into_times(self) -> Self::TimeIter;
}

impl<const N: usize> IntoTimes for [f64; N] {
	type TimeIter = std::array::IntoIter<LinearTime, N>;
	fn into_times(self) -> Self::TimeIter {
		let mut times = self.map(LinearTime::from_secs_f64);
		times.sort_unstable();
		times.into_iter()
	}
}

/// Interface for creating [`Time`] values to override conversion.
#[derive(Copy, Clone, Debug, Default)]
pub struct LinearTime(f64);

impl LinearTime {
	pub fn from_secs_f64(secs: f64) -> Self {
		Self(secs)
	}
	
	/// Conversion into [`Time`], but always rounds down.
	pub(crate) fn try_into_time(self, basis: Time) -> Result<Time, Time> {
		let LinearTime(mut t) = self;
		let sign = t.signum();
		t *= sign;
		
		const MANT_MASK: u64 = (1 << 52) - 1;
		const EXP_MASK: u64 = (1 << 11) - 1;
		
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
			if sign == -1. && (t * 1e9).fract() != 0. {
				nanos += 1;
			}
			Time::new(0, nanos)
		} else if exp < 52 {
			let secs = mant >> (52 - exp);
			let nanos_tmp = (((mant << exp) & MANT_MASK) as u128) * 1_000_000_000;
			let mut nanos = (nanos_tmp >> 52) as u32;
			if sign == -1. && (t * 1e9).fract() != 0. {
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
}

impl Ord for LinearTime {
	fn cmp(&self, other: &Self) -> Ordering {
		self.0.total_cmp(&other.0)
	}
}

impl PartialOrd for LinearTime {
	fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		Some(self.cmp(other))
	}
}

impl Eq for LinearTime {}

impl PartialEq for LinearTime {
	fn eq(&self, other: &Self) -> bool {
		self.cmp(other) == Ordering::Equal
	}
}

/// ...
#[derive(Clone)]
pub struct RootFilterMap<T> {
	pub(crate) times: T,
	pub(crate) basis: Time,
}

impl<T> Iterator for RootFilterMap<T>
where
	T: Iterator<Item = LinearTime>,
{
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		while let Some(root) = self.times.next() {
			if let Ok(time) = root.try_into_time(self.basis) {
				return Some(time)
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.times.size_hint().1)
	}
}