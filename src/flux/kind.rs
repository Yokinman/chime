//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::{Debug, Formatter};
use std::ops::{Add, Mul, Sub};

use crate::linear::{Linear, LinearPlus, LinearPlusArray, Scalar, Vector};
use crate::time::Time;
use crate::{Change, Constant, Flux, FluxValue};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind: Clone + Debug + 'static {
	type Value: LinearPlus;
	
	const DEGREE: usize;
	
	fn from_value(value: <Self::Value as LinearPlus>::Inner) -> Self;
	
	fn add_value(self, value: <Self::Value as LinearPlus>::Inner) -> Self;
	
	fn deriv(self) -> Self;
	
	fn eval(&self, time: Scalar) -> <Self::Value as LinearPlus>::Inner;
	
	fn to_time(self, time: Scalar) -> Self;
	
	/// The order at or immediately preceding the value at a time.
	/// 
	/// This should be the first non-zero [`FluxKind::eval`] value of this kind
	/// or its derivatives; reversed for odd derivatives.
	fn initial_order(&self, time: Scalar) -> Option<Ordering>
	where
		<Self::Value as LinearPlus>::Inner: PartialOrd
	{
		// !!! Alternative: Translate polynomial using `to_time` and then check
		// leading terms in order. Unknown which is more precise/faster.
		
		use std::borrow::Cow;
		
		let mut deriv = Cow::Borrowed(self);
		
		for degree in 0..=Self::DEGREE {
			let order = deriv.eval(time)
				.partial_cmp(&Linear::zero());
			
			if order != Some(Ordering::Equal) || degree == Self::DEGREE {
				return if degree % 2 == 0 {
					order
				} else {
					order.map(Ordering::reverse)
				}
			}
			
			deriv = match deriv {
				Cow::Borrowed(x) => Cow::Owned(x.clone().deriv()),
				Cow::Owned(x) => Cow::Owned(x.deriv()),
			};
		}
		
		None
	}
	
	fn zero() -> Self {
		Self::from_value(Linear::zero())
	}
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
	
	fn value(&self) -> <Self::Value as LinearPlus>::Inner {
		self.eval(Scalar::from(0.))
	}
}

/// ...
pub trait FluxIntegral: FluxKind {
	type Integ: FluxKind<Value: LinearPlus<Inner = KindLinear<Self>>>;
	fn integ(self) -> Self::Integ;
}

/// Shortcut for [`Flux::change`] parameter.
pub type EmptyFluxAccum<T> = FluxAccum<Constant<<T as FluxKind>::Value>>;

/// ... [`FluxKind::Accum`]
pub struct FluxAccum<K> {
	/// Accumulated polynomial.
	pub poly: Poly<K>,
	
	/// The basis time of the changes being accumulated.
	pub time: Time,
}

impl<K: FluxKind> FluxAccum<K> {
	pub fn new(poly: Poly<K>, time: Time) -> Self {
		Self { poly, time }
	}
}

impl<A, B> Add<Change<B>> for FluxAccum<A>
where
	A: FluxKind + ops::Add<<B::Kind as FluxIntegral>::Integ>,
	B: Flux<Kind: FluxIntegral>,
	<B::Kind as FluxIntegral>::Integ:
		Mul<Scalar, Output = <B::Kind as FluxIntegral>::Integ>,
	// ??? ^ Maybe I should apply the Scalar some other way.
{
	type Output = FluxAccum<<A as ops::Add<<B::Kind as FluxIntegral>::Integ>>::Output>;
	fn add(self, rhs: Change<B>) -> Self::Output {
		let FluxAccum { poly, time } = self;
		let Change { rate, unit } = rhs;
		let poly_time = poly.time;
		let sub_poly = FluxValue::new(rate, time).poly(poly_time);
		let time_scale = unit.as_secs_f64().recip();
		FluxAccum {
			poly: poly + (sub_poly.integ() * Scalar::from(time_scale)),
			time,
		}
	}
}

impl<A, B> Sub<Change<B>> for FluxAccum<A>
where
	A: FluxKind + ops::Add<<B::Kind as FluxIntegral>::Integ>,
	B: Flux<Kind: FluxIntegral>,
	<B::Kind as FluxIntegral>::Integ:
		Mul<Scalar, Output = <B::Kind as FluxIntegral>::Integ>,
	// ??? ^ Maybe I should apply the Scalar some other way.
{
	type Output = FluxAccum<<A as ops::Add<<B::Kind as FluxIntegral>::Integ>>::Output>;
	fn sub(self, rhs: Change<B>) -> Self::Output {
		let FluxAccum { poly, time } = self;
		let Change { rate, unit } = rhs;
		let sub_poly = FluxValue::new(rate, time).poly(poly.time);
		let time_scale = unit.as_secs_f64().recip();
		FluxAccum {
			poly: poly + (sub_poly.integ() * Scalar::from(time_scale * -1.)),
			time,
		}
	}
}

impl<T: FluxKind, const SIZE: usize> FluxKind for [T; SIZE] {
	type Value = LinearPlusArray<T::Value, SIZE>;
	const DEGREE: usize = T::DEGREE;
	fn from_value(value: <Self::Value as LinearPlus>::Inner) -> Self {
		value.map(|x| T::from_value(x))
	}
	fn add_value(self, value: <Self::Value as LinearPlus>::Inner) -> Self {
		let mut values = value.into_iter();
		self.map(|x| x.add_value(values.next().unwrap()))
	}
	fn deriv(self) -> Self {
		self.map(T::deriv)
	}
	fn eval(&self, time: Scalar) -> <Self::Value as LinearPlus>::Inner {
		self.each_ref().map(|x| T::eval(x, time))
	}
	fn to_time(self, time: Scalar) -> Self {
		self.map(|x| x.to_time(time))
	}
}

/// Shortcut for the inner [`Linear`] type of a [`FluxKind`].
pub(crate) type KindLinear<T> = <<T as FluxKind>::Value as LinearPlus>::Inner;

/// Combining [`FluxKind`] types.
/// 
/// Primarily this serves as a way to put two kinds of change-over-time into
/// the same space, for combination or comparison purposes.
pub mod ops {
	use std::ops;
	use super::{FluxKind, KindLinear};
	use crate::linear::{LinearPlus, Scalar};
	
	/// Adding two kinds of change.
	pub trait Add<K: FluxKind = Self>: FluxKind {
		type Output: FluxKind<Value: LinearPlus<Inner = KindLinear<K>>>;
		fn add(self, kind: K) -> <Self as Add<K>>::Output;
	}
	
	impl<A, B> Add<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value: LinearPlus<Inner = KindLinear<A>>>,
		<A as ops::Add<B>>::Output: FluxKind<Value: LinearPlus<Inner = KindLinear<A>>>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn add(self, kind: B) -> <A as ops::Add<B>>::Output {
			self + kind
		}
	}
	
	/// Differentiating two kinds of change.
	pub trait Sub<K: FluxKind = Self>: FluxKind {
		type Output: FluxKind<Value: LinearPlus<Inner = KindLinear<K>>>;
		fn sub(self, kind: K) -> <Self as Sub<K>>::Output;
	}
	
	impl<A, B> Sub<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value: LinearPlus<Inner = KindLinear<A>>>
			+ ops::Mul<Scalar, Output=B>,
		<A as ops::Add<B>>::Output: FluxKind<Value: LinearPlus<Inner = KindLinear<A>>>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn sub(self, kind: B) -> <A as ops::Add<B>>::Output {
			// ??? This could pass through to `ops::Sub` directly, but then
			// stuff like [`sum::Sum`] would need a whole extra set of macro
			// implementations. For now, this will just reuse `ops::Add`.
			self + (kind * Scalar::from(-1.))
		}
	}
	
	/// Squaring a kind of change.
	pub trait Sqr: FluxKind {
		type Output: FluxKind<Value: LinearPlus<Inner = KindLinear<Self>>>;
		fn sqr(self) -> <Self as Sqr>::Output;
	}
	
	impl<K: FluxKind> Sqr for K
	where
		K: ops::Mul,
		<K as ops::Mul>::Output: FluxKind<Value: LinearPlus<Inner = KindLinear<K>>>,
	{
		type Output = <K as ops::Mul>::Output;
		fn sqr(self) -> <Self as Sqr>::Output {
			self.clone() * self
		}
	}
}

/// A [`FluxKind`] paired with a basis time.
/// e.g. `Poly<1 + 2x>` => `1 + 2(x-time)`.
pub struct Poly<K> {
	pub kind: K,
	pub time: Time,
}

/// ... [`<Poly as IntoIterator>::IntoIter`]
pub struct PolyIter<T> {
	iter: T,
	time: Time,
}

impl<K> Clone for Poly<K>
where
	K: Clone
{
	fn clone(&self) -> Self {
		Self {
			kind: self.kind.clone(),
			time: self.time,
		}
	}
}

impl<K> Copy for Poly<K>
where
	K: Copy
{}

impl<K> Debug for Poly<K>
where
	K: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		f.debug_tuple("Poly")
			.field(&self.kind)
			.field(&self.time)
			.finish()
	}
}

impl<K> PartialEq for Poly<K>
where
	K: PartialEq
{
	fn eq(&self, other: &Self) -> bool {
		self.kind == other.kind && self.time == other.time
	}
}

impl<K: FluxKind> Poly<K> {
	pub fn new(kind: K, time: Time) -> Self {
		Self { kind, time }
	}
	
	pub fn with_value(value: K::Value) -> Self {
		Self::new(K::from_value(value.into_inner()), Time::ZERO)
	}
	
	pub fn with_time(time: Time) -> Self {
		Self::new(K::zero(), time)
	}
}

impl<K: FluxKind> Poly<K> {
	pub fn eval(&self, time: Time) -> <K::Value as LinearPlus>::Inner {
		self.kind.eval(Scalar::from(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
	
	pub fn at(&self, time: Time) -> K::Value {
		K::Value::from_inner(self.eval(time))
	}
	
	pub fn rate_at(&self, time: Time) -> K::Value {
		K::Value::from_inner(self.kind.clone().deriv().eval(Scalar::from(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		})))
	}
	
	pub fn to_time(mut self, time: Time) -> Self {
		if self.time != time {
			self.kind = self.kind.to_time(Scalar::from(if time > self.time {
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
		KindLinear<K>: PartialOrd
	{
		self.kind.initial_order(Scalar::from(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}))
	}
	
	pub fn integ(self) -> Poly<K::Integ>
	where
		K: FluxIntegral,
	{
		Poly::new(self.kind.integ(), self.time)
	}
}

impl<K: FluxKind> Poly<K> {
	pub fn sqr(self) -> Poly<<K as ops::Sqr>::Output>
	where
		K: ops::Sqr
	{
		Poly::new(self.kind.sqr(), self.time)
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	pub(crate) fn when_sign<F>(self, order: Ordering, filter: F) -> crate::pred::PredFilter<crate::pred::Pred<K>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialEq,
		KindLinear<K>: PartialOrd,
	{
		let pred = crate::pred::Pred {
			poly: self,
			order
		};
		crate::pred::PredFilter { pred, filter }
	}
	
	/// Times when the value is equal to zero.
	pub(crate) fn when_zero<F>(self, filter: F) -> crate::pred::PredFilter<crate::pred::PredEq<K>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialEq,
		KindLinear<K>: PartialEq,
	{
		crate::pred::PredFilter {
			pred: crate::pred::PredEq { poly: self },
			filter,
		}
	}
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Poly::new(K::zero(), Time::ZERO)
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
		Poly {
			kind: self.kind.add(rhs.to_time(self.time).kind),
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
			kind: self.kind.sub(rhs.to_time(self.time).kind),
			time: self.time,
		}
	}
}

impl<K> Mul<Scalar> for Poly<K>
where
	K: Mul<Scalar, Output=K>
{
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.kind = self.kind * rhs;
		self
	}
}

impl<K, const SIZE: usize> Vector<SIZE> for Poly<K>
where
	K: Vector<SIZE>,
{
	type Output = Poly<K::Output>;
	fn index(&self, index: usize) -> Self::Output {
		Poly {
			kind: self.kind.index(index),
			time: self.time,
		}
	}
}

impl<K> IntoIterator for Poly<K>
where
	K: IntoIterator<Item: FluxKind>,
{
	type Item = Poly<K::Item>;
	type IntoIter = PolyIter<K::IntoIter>;
	fn into_iter(self) -> Self::IntoIter {
		PolyIter {
			iter: self.kind.into_iter(),
			time: self.time,
		}
	}
}

impl<T> Iterator for PolyIter<T>
where
	T: Iterator<Item: FluxKind>,
{
	type Item = Poly<T::Item>;
	fn next(&mut self) -> Option<Self::Item> {
		self.iter.next()
			.map(|x| Poly::new(x, self.time))
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.iter.size_hint()
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
	/// 
	/// Consistently rounding down is important so that tiny ranges of time
	/// aren't ignored due to rounding. It's also better if predictions catch
	/// events before they happen rather than after.
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
	pub(crate) prev_time: Time,
}

impl<T> Iterator for RootFilterMap<T>
where
	T: Iterator<Item = LinearTime>,
{
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		for root in self.times.by_ref() {
			if let Ok(time) = root.try_into_time(self.basis) {
				let prev_time = self.prev_time;
				self.prev_time = time;
				return Some(time - prev_time)
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.times.size_hint().1)
	}
}