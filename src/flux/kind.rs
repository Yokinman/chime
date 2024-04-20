//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub};

use crate::linear::{Linear, LinearIso, LinearIsoVec, LinearVec, Scalar};
use crate::time;
use crate::time::{Time, TimeRanges};

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
	fn when_sign<F>(self, order: Ordering, filter: F) -> PredFilter<Pred<K, I>, F>
	where
		F: TimeFilterMap + 'static,
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let pred = Pred {
			poly: self,
			order
		};
		PredFilter { pred, filter }
	}
	
	/// Times when the value is equal to zero.
	fn when_zero<F>(self, filter: F) -> PredFilter<PredEq<K, I>, F>
	where
		F: TimeFilterMap + 'static,
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		PredFilter {
			pred: PredEq { poly: self },
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
	type TimeIter: Iterator<Item=LinearTime> + Send + Sync + Clone;
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
	times: T,
	basis: Time,
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

/// Function that converts a root value to a Time, or ignores it.
pub(crate) trait TimeFilterMap: Clone + Send + Sync {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time>;
}

impl<T: Fn(Time, bool) -> Option<Time> + Clone + Send + Sync> TimeFilterMap for T {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time> {
		self(time, is_end)
	}
}

/// ...
#[derive(Clone)]
pub struct PreTimeFilterMap;

impl TimeFilterMap for PreTimeFilterMap {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time> {
		if is_end {
			Some(time)
		} else {
			time.checked_sub(time::NANOSEC)
		}
	}
}

/// ...
pub struct DiffTimeFilterMap<A, B, D, I, J, L> {
	a_poly: Poly<A, I>,
	b_poly: Poly<B, J>,
	diff_poly: Poly<D, L>,
}

impl<A, B, D, I, J, L> Clone for DiffTimeFilterMap<A, B, D, I, J, L>
where
	A: Clone,
	B: Clone,
	D: Clone,
{
	fn clone(&self) -> Self {
		Self {
			a_poly: self.a_poly.clone(),
			b_poly: self.b_poly.clone(),
			diff_poly: self.diff_poly.clone(),
		}
	}
}

impl<A, B, D, I, J, L> TimeFilterMap for DiffTimeFilterMap<A, B, D, I, J, L>
where
	A: FluxKind,
	B: FluxKind<Value=A::Value>,
	D: FluxKind<Value=A::Value>,
	A::Value: PartialEq,
	I: LinearIso<A::Value>,
	J: LinearIso<A::Value>,
	L: LinearIso<A::Value>,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		
		let Self { a_poly, b_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).sign();
		
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = if is_end {
				time.checked_add(inc_time)
			} else {
				time.checked_sub(inc_time)
			} {
				 // Stop Before Rate Reverses:
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				if I::linear_id(a_poly.at(next_time)) != J::linear_id(b_poly.at(next_time)) {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub struct DisTimeFilterMap<const SIZE: usize, A, B, D, E, F, I, J, L>
where
	E: FluxKind,
	F: FluxKind,
{
	a_pos: PolyVec<SIZE, A, I>,
	b_pos: PolyVec<SIZE, B, J>,
	dis_poly: Poly<D, L>,
	pos_poly: Poly<E, E::Value>,
	diff_poly: Poly<F, F::Value>,
}

impl<const SIZE: usize, A, B, D, E, F, I, J, L> Clone
	for DisTimeFilterMap<SIZE, A, B, D, E, F, I, J, L>
where
	A: Clone,
	B: Clone,
	D: Clone,
	E: FluxKind,
	F: FluxKind,
{
	fn clone(&self) -> Self {
		Self {
			a_pos: self.a_pos.clone(),
			b_pos: self.b_pos.clone(),
			dis_poly: self.dis_poly.clone(),
			pos_poly: self.pos_poly.clone(),
			diff_poly: self.diff_poly.clone(),
		}
	}
}

impl<const SIZE: usize, A, B, D, E, F, I, J, L> TimeFilterMap
	for DisTimeFilterMap<SIZE, A, B, D, E, F, I, J, L>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	D: FluxKind,
	D::Value: PartialEq + Mul<Output=D::Value>,
	A::Kind: FluxKind<Value=D::Value>,
	B::Kind: FluxKind<Value=D::Value>,
	E: FluxKind<Value=D::Value>,
	F: FluxKind<Value=D::Value>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	L: LinearIso<D::Value>,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		// To handle rounding, the lower bound of equality is undershot.
		// For example, a pair of IVec2 points can round towards each other up
		// to `0.5` along each axis, or `sqrt(n)` in n-dimensional distance. 
		
		let Self { a_pos, b_pos, dis_poly, pos_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).sign();
		
		 // Rounding Buffer:
		if !is_end {
			let round_factor = Scalar(0.5 / (SIZE as f64).sqrt());
			loop {
				let mut inc_time = time::NANOSEC;
				while let Some(next_time) = time.checked_sub(inc_time) {
					 // Stop Before Rate Reverses:
					let rate = diff_poly.rate_at(next_time);
					if sign != rate.sign() && !rate.is_zero() {
						if inc_time == time::NANOSEC {
							return Some(time)
						}
						break
					}
					
					 // Calculate Actual Distances:
					let dis = dis_poly.at(next_time);
					let mut a_dis = D::Value::zero();
					let mut b_dis = D::Value::zero();
					let mut real_diff = D::Value::zero();
					for i in 0..SIZE {
						let a = a_pos.index_poly(i).at(next_time);
						let b = b_pos.index_poly(i).at(next_time);
						a_dis = a_dis + a*a;
						b_dis = b_dis + b*b;
						let x = a - b;
						real_diff = real_diff + x*x;
					}
					a_dis = I::Value::linear_id(a_dis.sqrt());
					b_dis = J::Value::linear_id(b_dis.sqrt());
					real_diff = Mul::<Scalar>::mul(real_diff.sqrt() - dis, round_factor);
					let c_dis = L::linear_id(dis);
					
					 // Undershoot Actual Distances:
					if
						a_dis != I::Value::linear_id(a_dis + real_diff) &&
						b_dis != J::Value::linear_id(b_dis + real_diff) &&
						c_dis != L::linear_id(c_dis + real_diff)
					{
						 // Undershoot Predicted Distances:
						let pred_diff = Mul::<Scalar>::mul(
							pos_poly.at(next_time).sqrt() - dis,
							round_factor
						);
						if
							a_dis != I::Value::linear_id(a_dis + pred_diff) &&
							b_dis != J::Value::linear_id(b_dis + pred_diff) &&
							c_dis != L::linear_id(c_dis + pred_diff)
						{
							break
						}
					}
					
					time = next_time;
					inc_time += inc_time;
				}
				if inc_time == time::NANOSEC {
					break
				}
			}
			
			return Some(time)
		}
		
		 // Fully Cover Zero Range:
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = time.checked_add(inc_time) {
				 // Stop Before Rate Reverses:
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				let mut pos = D::Value::zero();
				for i in 0..SIZE {
					let x = I::Value::linear_id(a_pos.index_poly(i).at(next_time))
						- J::Value::linear_id(b_pos.index_poly(i).at(next_time));
					pos = pos + x*x;
				}
				let dis = L::linear_id(dis_poly.at(next_time));
				if pos != dis*dis {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub trait Prediction {
	type TimeRanges: Iterator<Item = (Time, Time)>;
	fn into_time_ranges(self) -> Self::TimeRanges;
	
	/// Decrements the lower bound of each range by 1 nanosecond.  
	fn pre(self) -> PredFilter<Self, PreTimeFilterMap>
	where
		Self: Sized
	{
		PredFilter {
			pred: self,
			filter: PreTimeFilterMap,
		}
	}
}

macro_rules! impl_prediction {
	(for<$($param:ident),*> $pred:ty) => {
		impl<$($param),*> Prediction for $pred
		where
			Self: IntoIterator<Item = (Time, Time)>
		{
			type TimeRanges = <Self as IntoIterator>::IntoIter;
			fn into_time_ranges(self) -> Self::TimeRanges {
				self.into_iter()
			}
		}
		
		impl<$($param,)* T> std::ops::BitAnd<T> for $pred {
			type Output = PredInter<Self, T>;
			fn bitand(self, b_pred: T) -> Self::Output {
				PredInter {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param,)* T> std::ops::BitOr<T> for $pred {
			type Output = PredUnion<Self, T>;
			fn bitor(self, b_pred: T) -> Self::Output {
				PredUnion {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param,)* T> std::ops::BitXor<T> for $pred {
			type Output = PredSymDiff<Self, T>;
			fn bitxor(self, b_pred: T) -> Self::Output {
				PredSymDiff {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param),*> std::ops::Not for $pred {
			type Output = PredInv<Self>;
			fn not(self) -> Self::Output {
				PredInv {
					pred: self,
				}
			}
		}
	};
	($(for<$($param:ident),*> $pred:ty;)+) => {
		$(impl_prediction!{for<$($param),*> $pred})+
	};
}

impl_prediction!{
	for<K, I> Pred<K, I>;
	for<K, I> PredEq<K, I>;
	for<P, F> PredFilter<P, F>;
	for<A, B> PredInter<A, B>;
	for<A, B> PredUnion<A, B>;
	for<A, B> PredSymDiff<A, B>;
	for<P> PredInv<P>;
	for<> DynPred;
}

impl Prediction for std::ops::Range<Time> {
	type TimeRanges = TimeRanges<std::iter::Once<time::TimeRange>>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		TimeRanges::from_range(self)
	}
}

impl Prediction for std::ops::RangeInclusive<Time> {
	type TimeRanges = TimeRanges<std::iter::Once<time::TimeRange>>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		TimeRanges::from_range(self)
	}
}

impl Prediction for std::ops::RangeTo<Time> {
	type TimeRanges = TimeRanges<std::iter::Once<time::TimeRange>>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		TimeRanges::from_range(self)
	}
}

impl Prediction for std::ops::RangeToInclusive<Time> {
	type TimeRanges = TimeRanges<std::iter::Once<time::TimeRange>>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		TimeRanges::from_range(self)
	}
}

impl Prediction for std::ops::RangeFull {
	type TimeRanges = TimeRanges<std::iter::Once<time::TimeRange>>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		TimeRanges::from_range(self)
	}
}

impl<P: Prediction> Prediction for Option<P> {
	type TimeRanges = OptionTimeRanges<P::TimeRanges>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		OptionTimeRanges {
			times: self.map(|x| x.into_time_ranges())
		}
	}
}

/// ...
pub struct OptionTimeRanges<T> {
	times: Option<T>,
}

impl<T> Iterator for OptionTimeRanges<T>
where
	T: Iterator<Item = (Time, Time)>
{
	type Item = T::Item;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(times) = self.times.as_mut() {
			times.next()
		} else {
			None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		if let Some(times) = self.times.as_ref() {
			times.size_hint()
		} else {
			(0, Some(0))
		}
	}
}

/// ...
pub struct DynTimeRanges {
	inner: Box<dyn Iterator<Item = (Time, Time)>>,
}

impl Iterator for DynTimeRanges {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		self.inner.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.inner.size_hint()
	}
}

/// ...
pub struct Pred<K, I> {
	poly: Poly<K, I>,
	order: Ordering,
}

impl<K, I> Clone for Pred<K, I>
where
	K: Clone,
{
	fn clone(&self) -> Self {
		Self {
			poly: self.poly.clone(),
			order: self.order,
		}
	}
}

impl<K, I> IntoIterator for Pred<K, I>
where
	K: Roots + PartialOrd,
	K::Value: PartialOrd,
	I: LinearIso<K::Value>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>>;
	fn into_iter(self) -> Self::IntoIter {
		let basis = self.poly.time;
		let basis_order = self.poly
			.initial_order(Time::ZERO)
			.unwrap_or(Ordering::Equal);
		let times = RootFilterMap {
			times: self.poly.inner.roots().into_times(),
			basis,
		};
		TimeRanges::new(times, basis_order, self.order)
	}
}

/// ...
pub struct PredEq<K, I> {
	poly: Poly<K, I>,
}

impl<K, I> Clone for PredEq<K, I>
where
	K: Clone,
{
	fn clone(&self) -> Self {
		Self {
			poly: self.poly.clone(),
		}
	}
}

impl<K, I> IntoIterator for PredEq<K, I>
where
	K: Roots + PartialEq,
	K::Value: PartialEq,
	I: LinearIso<K::Value>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>>;
	fn into_iter(self) -> Self::IntoIter {
		let basis = self.poly.time;
		let basis_order = if self.poly.inner.is_zero() {
			Ordering::Equal
		} else {
			Ordering::Greater
		};
		let times = RootFilterMap {
			times: self.poly.inner.roots().into_times(),
			basis,
		};
		TimeRanges::new(times, basis_order, Ordering::Equal)
	}
}

/// ...
trait AnyPrediction {
	fn as_time_ranges(&mut self) -> DynTimeRanges;
}

impl<P> AnyPrediction for P
where
	P: Prediction + Clone + 'static
{
	fn as_time_ranges(&mut self) -> DynTimeRanges {
		DynTimeRanges {
			inner: Box::new(self.clone().into_time_ranges())
		}
	}
}

/// ...
pub struct DynPred {
	inner: Box<dyn AnyPrediction>,
}

impl DynPred {
	pub fn new(pred: impl Prediction + Clone + 'static) -> Self {
		Self {
			inner: Box::new(pred)
		}
	}
}

impl IntoIterator for DynPred {
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = DynTimeRanges;
	fn into_iter(mut self) -> Self::IntoIter {
		self.inner.as_time_ranges()
	}
}

/// ...
#[derive(Clone)]
pub struct PredFilter<P, F> {
	pred: P,
	filter: F,
}

impl<I, P, F> IntoIterator for PredFilter<P, F>
where
	I: time::TimeRangeIter,
	P: Prediction<TimeRanges = TimeRanges<I>>,
	F: TimeFilterMap + 'static,
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges<time::TimeFilter<I, F>>;
	fn into_iter(self) -> Self::IntoIter {
		self.pred.into_time_ranges()
			.into_filtered(self.filter)
	}
}

/// ...
#[derive(Clone)]
pub struct PredInter<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B, I, J> IntoIterator for PredInter<A, B>
where
	I: time::TimeRangeIter,
	J: time::TimeRangeIter,
	A: Prediction<TimeRanges = TimeRanges<I>>,
	B: Prediction<TimeRanges = TimeRanges<J>>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::InterTimeRanges<I, J>>;
	fn into_iter(self) -> Self::IntoIter {
		self.a_pred.into_time_ranges()
			.inter(self.b_pred.into_time_ranges())
	}
}

/// ...
#[derive(Clone)]
pub struct PredUnion<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B, I, J> IntoIterator for PredUnion<A, B>
where
	I: time::TimeRangeIter,
	J: time::TimeRangeIter,
	A: Prediction<TimeRanges = TimeRanges<I>>,
	B: Prediction<TimeRanges = TimeRanges<J>>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::UnionTimeRanges<I, J>>;
	fn into_iter(self) -> Self::IntoIter {
		self.a_pred.into_time_ranges()
			.union(self.b_pred.into_time_ranges())
	}
}

/// ...
#[derive(Clone)]
pub struct PredSymDiff<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B, I, J> IntoIterator for PredSymDiff<A, B>
where
	I: time::TimeRangeIter,
	J: time::TimeRangeIter,
	A: Prediction<TimeRanges = TimeRanges<I>>,
	B: Prediction<TimeRanges = TimeRanges<J>>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::DiffTimeRanges<I, J>>;
	fn into_iter(self) -> Self::IntoIter {
		self.a_pred.into_time_ranges()
			.sym_diff(self.b_pred.into_time_ranges())
	}
}

/// ...
#[derive(Clone)]
pub struct PredInv<P> {
	pred: P,
}

impl<P, I> IntoIterator for PredInv<P>
where
	I: time::TimeRangeIter,
	P: Prediction<TimeRanges = TimeRanges<I>>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = TimeRanges<time::InvTimes<I>>;
	fn into_iter(self) -> Self::IntoIter {
		self.pred.into_time_ranges()
			.inv()
	}
}

/// [`crate::Flux::when`] predictive comparison.
pub trait When<B, J> {
	type Pred: Prediction;
	fn when(self, order: Ordering, poly: Poly<B, J>) -> Self::Pred;
}

impl<A, B, I, J> When<B, J> for Poly<A, I>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind,
	I: LinearIso<A::Value>,
	J: LinearIso<B::Value>,
	<A as ops::Sub<B>>::Output: Roots + PartialOrd,
	A::Value: PartialOrd,
{
	type Pred = PredFilter<
		Pred<<A as ops::Sub<B>>::Output, I>,
		DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output, I, J, I>,
	>;
	fn when(self, order: Ordering, poly: Poly<B, J>) -> Self::Pred {
		let diff_poly = self - poly;
		diff_poly
			.when_sign(order, DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B, J> {
	type Pred: Prediction;
	fn when_eq(self, poly: Poly<B, J>) -> Self::Pred;
}

impl<A, B, I, J> WhenEq<B, J> for Poly<A, I>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind,
	I: LinearIso<A::Value>,
	J: LinearIso<B::Value>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	type Pred = PredFilter<PredEq<<A as ops::Sub<B>>::Output, I>, DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output, I, J, I>>;
	fn when_eq(self, poly: Poly<B, J>) -> Self::Pred {
		let diff_poly = self - poly;
		diff_poly
			.when_zero(DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::FluxVec::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B, D, J, L> {
	type Pred: Prediction;
	fn when_dis(
		self,
		poly: PolyVec<SIZE, B, J>,
		order: Ordering,
		dis: Poly<D, L>,
	) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D, I, J, L> WhenDis<SIZE, B, D, J, L>
	for PolyVec<SIZE, A, I>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	A::Kind: ops::Sub<B::Kind>,
	<A::Kind as ops::Sub<B::Kind>>::Output: ops::Sqr,
	<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	<A::Kind as FluxKind>::Value:
		Mul<Output = <A::Kind as FluxKind>::Value> + PartialOrd,
	D: FluxKind<Value = <A::Kind as FluxKind>::Value> + ops::Sqr,
	L: LinearIso<D::Value>,
{
	type Pred = PredFilter<
		Pred<
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output as FluxKind>::Value,
		>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			I, J, L
		>
	>;
	fn when_dis(
		self,
		poly: PolyVec<SIZE, B, J>,
		order: Ordering,
		dis: Poly<D, L>,
	) -> Self::Pred {
		use ops::*;
		
		let basis = self.time();
		
		let mut sum = <<A::Kind as Sub<B::Kind>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).inner
				.sub(poly.index_poly(i).to_time(basis).inner)
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_sign(order, DisTimeFilterMap {
				a_pos: self,
				b_pos: poly,
				dis_poly: dis,
				pos_poly: sum,
				diff_poly,
			})
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D, J, L> {
	type Pred: Prediction;
	fn when_dis_eq(
		self,
		poly: PolyVec<SIZE, B, J>,
		dis: Poly<D, L>,
	) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D, I, J, L> WhenDisEq<SIZE, B, D, J, L> for PolyVec<SIZE, A, I>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	A::Kind: ops::Sub<B::Kind>,
	<A::Kind as ops::Sub<B::Kind>>::Output: ops::Sqr,
	<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	<A::Kind as FluxKind>::Value:
		Mul<Output = <A::Kind as FluxKind>::Value> + PartialEq,
	D: FluxKind<Value = <A::Kind as FluxKind>::Value> + ops::Sqr,
	L: LinearIso<D::Value>,
{
	type Pred = PredFilter<
		PredEq<
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output as FluxKind>::Value
		>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			I, J, L
		>
	>;
	fn when_dis_eq(
		self,
		poly: PolyVec<SIZE, B, J>,
		dis: Poly<D, L>,
	) -> Self::Pred {
		use ops::*;
		
		let basis = self.time();
		
		let mut sum = <<A::Kind as Sub<B::Kind>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).inner
				.sub(poly.index_poly(i).to_time(basis).inner)
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_zero(DisTimeFilterMap {
				a_pos: self,
				b_pos: poly,
				dis_poly: dis,
				pos_poly: sum,
				diff_poly,
			})
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
		poly.when_dis(
			PolyVec::new([Sum::<f64, 2>::zero(); 2], time),
			Ordering::Greater,
			Poly::new(crate::Constant::from(1.), time),
		).into_iter()
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