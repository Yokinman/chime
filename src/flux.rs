//! Utilities for describing how a type changes over time.

pub mod kind;
mod impls;

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::ops::{Add, Deref, DerefMut, Mul};

use crate::{
	linear::*,
	kind::{*, ops as kind_ops},
};

pub use self::impls::*;

pub use time::{Time, TimeUnit};

pub use flux_proc_macro::flux;

/// A discrete interface for a value that changes over time.
pub trait Moment {
	type Flux: Flux<Moment=Self>;
	
	/// Constructs the entirety of a [`Flux`] from a single moment.
	fn to_flux(&self, time: Time) -> Self::Flux;
}

/// The continuous interface for a value that changes over time.
pub trait Flux {
	// !!! Deriving PartialEq, Eq should count `f(t) = 1 + 2t` and
	// `g(t) = 3 + 2(t-base_time)` as the same Flux if `base_time = 1`.
	
	type Moment: Moment<Flux=Self>;
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// The output accumulator of [`Flux::change`].
	type OutAccum<'a>: FluxAccum<'a, Self::Kind>;
	
	/// An evaluation of this flux at some point in time.
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value;
	
	/// The time of [`Flux::base_value`].
	fn base_time(&self) -> Time;
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> Self::OutAccum<'a>;
	
	/// A moment in the timeline.
	fn to_moment(&self, time: Time) -> Self::Moment;
	
	/// Sets a moment in the timeline (affects all moments).
	fn set_moment(&mut self, time: Time, moment: Self::Moment)
	where
		Self: Sized
	{
		*self = moment.to_flux(time);
	}
	
	/// A reference to a moment in the timeline.
	fn at(&self, time: Time) -> MomentRef<Self::Moment> {
		let moment = self.to_moment(time);
		MomentRef {
			moment,
			borrow: PhantomData,
		}
	}
	
	/// A unique, mutable reference to a moment in the timeline.
	/// 
	/// ```text
	/// let mut moment = self.at_mut(time);
	/// // modifications
	/// ```
	/// 
	/// equivalent to:
	/// 
	/// ```text
	/// let mut moment = self.at(time);
	/// // modifications
	/// self.set_moment(time, moment);
	/// ```
	fn at_mut(&mut self, time: Time) -> MomentRefMut<Self::Moment> {
		let moment = Some(self.to_moment(time));
		MomentRefMut {
			moment,
			time,
			borrow: self,
		}
	}
	
	/// A point in the timeline.
	/// 
	/// `self.value(self.base_time()) == self.base_value()`
	fn value(&self, time: Time) -> <Self::Kind as FluxKind>::Value {
		let mut value = self.base_value();
		let base_time = self.base_time();
		if time != base_time {
			let accum = FluxAccumKind::Value {
				value: &mut value,
				depth: 0,
				time,
				base_time,
			};
			self.change(<Self::Kind as FluxKind>::Accum::from_kind(accum));
		}
		value
	}
	
	/// A polynomial description of this flux at the given time.
	fn poly(&self, time: Time) -> Poly<Self::Kind> {
		let mut poly = Poly::with_value(self.value(time));
		let accum = FluxAccumKind::Poly {
			poly: &mut poly,
			depth: 0,
			time,
			base_time: self.base_time(),
		};
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(accum));
		poly
	}
	
	/// Ranges when this is above/below/equal to another flux.
	fn when<T: Flux>(&self, cmp_order: Ordering, other: &T) -> TimeRanges
	where
		When<Self::Kind, T::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let time = self.base_time().max(other.base_time());
		When {
			a_poly: self.poly(time),
			b_poly: other.poly(time),
			cmp_order,
			time,
		}.into_iter()
	}
	
	/// Times when this is equal to another flux.
	fn when_eq<T: Flux>(&self, other: &T) -> Times
	where
		WhenEq<Self::Kind, T::Kind>: IntoIterator<IntoIter=Times>
	{
		let time = self.base_time().max(other.base_time());
		WhenEq {
			a_poly: self.poly(time),
			b_poly: other.poly(time),
			time,
		}.into_iter()
	}
}

/// Immutable moment-in-time interface for [`Flux::at`].
pub struct MomentRef<'b, M: Moment> {
	moment: M,
	borrow: PhantomData<&'b M::Flux>,
}

impl<M: Moment> Deref for MomentRef<'_, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		&self.moment
	}
}

impl<M: Moment> Debug for MomentRef<'_, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<M: Moment> Display for MomentRef<'_, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

/// Mutable moment-in-time interface for [`Flux::at_mut`].
pub struct MomentRefMut<'b, M: Moment> {
	moment: Option<M>,
	time: Time,
	borrow: &'b mut M::Flux,
}

impl<M: Moment> Drop for MomentRefMut<'_, M> {
	fn drop(&mut self) {
		if let Some(moment) = std::mem::take(&mut self.moment) {
			self.borrow.set_moment(self.time, moment);
		}
	}
}

impl<M: Moment> Deref for MomentRefMut<'_, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		if let Some(moment) = self.moment.as_ref() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<M: Moment> DerefMut for MomentRefMut<'_, M> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		if let Some(moment) = self.moment.as_mut() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<M: Moment> Debug for MomentRefMut<'_, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<M: Moment> Display for MomentRefMut<'_, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

/// Multidimensional change over time.
pub trait FluxVec {
	type Kind: FluxKind;
	fn times(&self) -> Times;
	fn polys(&self, time: Time) -> Polys<Self::Kind>;
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<T: FluxVec + ?Sized, D: Flux>(
		&self,
		other: &T,
		cmp_order: Ordering,
		dis: &D,
	) -> TimeRanges
	where
		WhenDis<Self::Kind, T::Kind, D::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let time = std::iter::once(dis.base_time())
			.chain(self.times())
			.chain(other.times())
			.max()
			.unwrap_or_default();
		
		WhenDis {
			a_poly: self.polys(time),
			b_poly: other.polys(time),
			cmp_order,
			dis: dis.poly(time),
			time,
		}.into_iter()
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis_eq<T: FluxVec + ?Sized, D: Flux>(
		&self,
		other: &T,
		dis: &D,
	) -> Times
	where
		WhenDisEq<Self::Kind, T::Kind, D::Kind>: IntoIterator<IntoIter=Times>
	{
		let time = std::iter::once(dis.base_time())
			.chain(self.times())
			.chain(other.times())
			.max()
			.unwrap_or_default();
		
		WhenDisEq {
			a_poly: self.polys(time),
			b_poly: other.polys(time),
			dis: dis.poly(time),
			time,
		}.into_iter()
	}
}

impl<T: Flux> FluxVec for [T] {
	type Kind = T::Kind;
	fn times(&self) -> Times {
		self.iter()
			.map(|x| x.base_time())
			.collect()
	}
	fn polys(&self, time: Time) -> Polys<Self::Kind> {
		self.iter()
			.map(|x| x.poly(time))
			.collect()
	}
}

impl<T: Flux> FluxVec for Vec<T> {
	type Kind = T::Kind;
	fn times(&self) -> Times {
		self.as_slice().times()
	}
	fn polys(&self, time: Time) -> Polys<Self::Kind> {
		self.as_slice().polys(time)
	}
}

impl<T: Flux, const S: usize> FluxVec for [T; S] {
	type Kind = T::Kind;
	fn times(&self) -> Times {
		self.as_slice().times()
	}
	fn polys(&self, time: Time) -> Polys<Self::Kind> {
		self.as_slice().polys(time)
	}
}

// impl<A: Flux, B: Flux> FluxVec for (A, B)
// where
// 	B::Kind: FluxKind<Value = <A::Kind as FluxKind>::Value>,
// 	A::Kind: Add<B::Kind>,
// 	<A::Kind as Add<B::Kind>>::Output: FluxKind<Value = <A::Kind as FluxKind>::Value>
// {
// 	type Kind = <A::Kind as Add<B::Kind>>::Output;
// 	fn times(&self) -> Times {
// 		[self.0.basis_time(), self.1.basis_time()]
// 			.into_iter().collect()
// 	}
// 	fn polys(&self, time: Time) -> Polys<Self::Kind> {
// 		[self.0.poly(time) + Poly::<B::Kind>::default(), Poly::<A::Kind>::default() + self.1.poly(time)]
// 			.into_iter().collect()
// 	}
// }

/// ...
#[allow(type_alias_bounds)]
type Polys<K: FluxKind> = Box<[Poly<K>]>;

/// Iterator of [`Time`] ranges.
#[derive(Clone)]
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

impl FromIterator<(Time, Time)> for TimeRanges {
	fn from_iter<T: IntoIterator<Item=(Time, Time)>>(iter: T) -> Self {
		Self(iter.into_iter().collect::<Vec<_>>().into_iter())
	}
}

/// Iterator of [`Time`] values.
#[derive(Clone)]
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

impl FromIterator<Time> for Times {
	fn from_iter<T: IntoIterator<Item=Time>>(iter: T) -> Self {
		Self(iter.into_iter().collect::<Vec<_>>().into_iter())
	}
}

/// [`Flux::when`] predictive comparison.
#[derive(Copy, Clone, Debug)]
pub struct When<A: FluxKind, B: FluxKind> {
	a_poly: Poly<A>,
	b_poly: Poly<B>,
	cmp_order: Ordering,
	time: Time,
}

impl<A: FluxKind, B: FluxKind> IntoIterator for When<A, B>
where
	A: kind_ops::Sub<B>,
	<A as kind_ops::Sub<B>>::Output: Roots + PartialOrd,
	A::Value: PartialOrd,
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = self.a_poly - self.b_poly;
		poly.when_sign(self.cmp_order, self.time)
	}
}

/// [`Flux::when_eq`] predictive equality comparison.
#[derive(Copy, Clone, Debug)]
pub struct WhenEq<A: FluxKind, B: FluxKind> {
	a_poly: Poly<A>,
	b_poly: Poly<B>,
	time: Time,
}

impl<A: FluxKind, B: FluxKind> IntoIterator for WhenEq<A, B>
where
	A: kind_ops::Sub<B>,
	<A as kind_ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	type Item = Time;
	type IntoIter = Times;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = self.a_poly - self.b_poly;
		poly.when_zero(self.time)
	}
}

/// [`FluxVec::when_dis`] predictive distance comparison.
pub struct WhenDis<A: FluxKind, B: FluxKind, D: FluxKind> {
	a_poly: Box<[Poly<A>]>,
	b_poly: Box<[Poly<B>]>,
	cmp_order: Ordering,
	dis: Poly<D>,
	time: Time,
}

impl<A: FluxKind, B: FluxKind, D> IntoIterator for WhenDis<A, B, D>
where
	A: kind_ops::Sub<B>,
	<A as kind_ops::Sub<B>>::Output: kind_ops::Sqr,
	<<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output:
		Add<Output = <<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
		+ kind_ops::Sub<
			<D as kind_ops::Sqr>::Output,
			Output = <<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	A::Value: PartialOrd,
	D: FluxKind<Value=A::Value> + kind_ops::Sqr,
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges;
	
	fn into_iter(self) -> Self::IntoIter {
		let count = self.a_poly.len().max(self.b_poly.len());
		
		let mut a_iter = self.a_poly.iter().copied();
		let mut b_iter = self.b_poly.iter().copied();
		
		let mut sum = Poly
			::<<<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
			::default();
		
		for _ in 0..count {
			let a = a_iter.next().unwrap_or_default();
			let b = b_iter.next().unwrap_or_default();
			sum = sum + (a - b).sqr();
		}
		sum = sum - self.dis.sqr();
		
		sum.when_sign(self.cmp_order, self.time)
	}
}

/// [`FluxVec::when_dis_eq`] predictive distance comparison.
pub struct WhenDisEq<A: FluxKind, B: FluxKind, D: FluxKind> {
	a_poly: Box<[Poly<A>]>,
	b_poly: Box<[Poly<B>]>,
	dis: Poly<D>,
	time: Time,
}

impl<A: FluxKind, B: FluxKind, D> IntoIterator for WhenDisEq<A, B, D>
where
	A: kind_ops::Sub<B>,
	<A as kind_ops::Sub<B>>::Output: kind_ops::Sqr,
	<<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output:
		Add<Output = <<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
		+ kind_ops::Sub<
			<D as kind_ops::Sqr>::Output,
			Output = <<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	A::Value: PartialEq,
	D: FluxKind<Value=A::Value> + kind_ops::Sqr,
{
	type Item = Time;
	type IntoIter = Times;
	
	fn into_iter(self) -> Self::IntoIter {
		let count = self.a_poly.len().max(self.b_poly.len());
		
		let mut a_iter = self.a_poly.iter().copied();
		let mut b_iter = self.b_poly.iter().copied();
		
		let mut sum = Poly
			::<<<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output>
			::default();
		
		for _ in 0..count {
			let a = a_iter.next().unwrap_or_default();
			let b = b_iter.next().unwrap_or_default();
			sum = sum + (a - b).sqr();
		}
		sum = sum - self.dis.sqr();
		
		sum.when_zero(self.time)
	}
}

/// Used to construct a [`Change`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(TimeUnit::Secs)` 
pub trait Per: Sized {
	fn per(&self, unit: TimeUnit) -> Change<Self> {
		// !!! Make TimeUnit just a Duration
		Change {
			rate: self,
			unit
		}
	}
}

impl<T: Flux> Per for T {}

/// A description of a change over time for use with arithmetic operators.
pub struct Change<'t, T> {
	pub rate: &'t T,
	pub unit: TimeUnit,
}

/// No change over time.
/// 
/// Equivalent "constant" flux kinds should implement both `Into<Constant<T>>`
/// and `From<Constant<T>>` (e.g. `Sum<T, 0>`).
#[derive(Copy, Clone, Debug, Default)]
pub struct Constant<T> {
	value: T,
	time: Time,
}

impl<T> Deref for Constant<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.value
	}
}

impl<T> DerefMut for Constant<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.value
	}
}

impl<T: Linear> From<T> for Constant<T> {
	fn from(value: T) -> Self {
		Constant {
			value,
			time: Time::ZERO,
		}
	}
}

impl<T: Linear> Moment for T {
	type Flux = Constant<T>;
	fn to_flux(&self, time: Time) -> Self::Flux {
		Constant {
			value: *self,
			time,
		}
	}
}

impl<T: Linear> Flux for Constant<T> {
	type Moment = T;
	type Kind = Constant<T>;
	type OutAccum<'a> = ();
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		self.value
	}
	fn base_time(&self) -> Time {
		self.time
	}
	fn change<'a>(&self, _accum: <Self::Kind as FluxKind>::Accum<'a>) -> Self::OutAccum<'a> {}
	fn to_moment(&self, _time: Time) -> Self::Moment {
		self.value
	}
}

impl<T: Linear> FluxKind for Constant<T> {
	type Value = T;
	type Accum<'a> = ();
	fn value(&self) -> Self::Value {
		self.value
	}
	fn initial_order(&self) -> Option<Ordering> where T: PartialOrd {
		self.value.partial_cmp(&T::zero())
	}
	fn zero() -> Self {
		Self::from(T::zero())
	}
}

impl<T: Linear> Mul<Scalar> for Constant<T> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.value = self.value * rhs;
		self
	}
}

impl<T: Linear, K: FluxKind<Value=T>> Add<K> for Constant<T> {
	type Output = K;
	fn add(self, rhs: K) -> K {
		K::from(rhs.value() + self.value)
	}
}

impl<T> Mul<Constant<T>> for Constant<T>
where
	T: Mul<Output = T>
{
	type Output = Self;
	fn mul(mut self, rhs: Self) -> Self {
		self.value = self.value * rhs.value;
		self
	}
}

#[cfg(test)]
mod tests {
	use impl_op::impl_op;
	use super::*;
	use TimeUnit::*;
	use crate::sum::Sum;
	
	#[flux(Sum<f64, 4> = {value} + spd.per(Secs) + misc.per(Secs), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Pos {
		value: f64,
		spd: Spd,
		misc: Vec<Spd>,
	}
	
	#[flux(Sum<f64, 3> = {value} - fric.per(Secs) + accel.per(Secs), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Spd {
		value: f64,
		fric: Fric,
		accel: Accel,
	}
	
	#[flux(Sum<f64, 0> = {value}, crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Fric {
		value: f64,
	}
	
	#[flux(Sum<f64, 2> = {value} + jerk.per(Secs), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Accel {
		value: f64,
		jerk: Jerk,
	}
	
	#[flux(Sum<f64, 1> = {value} + snap.per(Secs), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Jerk {
		value: f64,
		snap: Snap,
	}
	
	#[flux(Constant<f64> = {value}, crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Snap {
		value: f64,
	}
	
	impl_op!{ *a -> f64 {
		Pos | Spd | Fric | Accel | Jerk | Snap => a.value
	}}
	
	fn position() -> <Pos as Moment>::Flux {
		let pos = Pos {
			value: 32.0, 
			spd: Spd {
				value: -4.4075 / 3.0,
				fric: Fric { value: 3.5 },
				accel: Accel {
					value: 2.0725 / 3.0,
					jerk: Jerk {
						value: 0.385,
						snap: Snap { value: -0.01 },
					},
				},
			},
			misc: Vec::new(),
		};
		let mut pos = pos.to_flux(Time::ZERO);
		let value = pos.value(10*Secs);
		pos.value.set_moment(10*Secs, value);
		pos.spd.accel.at_mut(20*Secs);
		pos.spd.accel.jerk.at_mut(10*Secs);
		pos
	}
	
	macro_rules! assert_times {
		($times:expr, $cmp_times:expr) => {{
			let times: Times = $times;
			let cmp_times = Times::from_iter($cmp_times.into_iter());
			assert_eq!(
				times.clone().count(),
				cmp_times.clone().count(),
				"a: {:?}, b: {:?}",
				times.collect::<Box<[_]>>(),
				cmp_times.collect::<Box<[_]>>(),
			);
			for (a, b) in times.zip(cmp_times) {
				let time = ((a.max(b) - b.min(a)).as_secs_f64() * 60.).floor();
				assert_eq!(time, 0., "a: {:?}, b: {:?}", a, b);
			}
		}};
	}
	
	macro_rules! assert_time_ranges {
		($ranges:expr, $cmp_ranges:expr) => {{
			let ranges: TimeRanges = $ranges;
			let cmp_ranges = TimeRanges::from_iter($cmp_ranges.into_iter());
			assert_eq!(
				ranges.clone().count(),
				cmp_ranges.clone().count(),
				"a: {:?}, b: {:?}",
				ranges.collect::<Box<[_]>>(),
				cmp_ranges.collect::<Box<[_]>>(),
			);
			for ((a, b), (x, y)) in ranges.zip(cmp_ranges) {
				let a_time = ((a.max(x) - x.min(a)).as_secs_f64() * 60.).floor();
				let b_time = ((b.max(y) - y.min(b)).as_secs_f64() * 60.).floor();
				assert_eq!(a_time, 0., "a: {:?}, x: {:?}", a, x);
				assert_eq!(b_time, 0., "b: {:?}, y: {:?}", b, y);
			}
		}};
	}
	
	macro_rules! assert_poly {
		($poly:expr, $cmp_poly:expr) => {{
			let poly = $poly;
			let cmp_poly = $cmp_poly;
			let dif_poly = poly - cmp_poly;
			for (degree, coeff) in dif_poly.into_iter().enumerate() {
				assert_eq!((coeff * 2_f64.powi((1 + degree) as i32)).round(), 0.,
					"poly: {:?}, cmp_poly: {:?}", poly, cmp_poly);
			}
		}};
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		
		 // Times:
		assert_eq!(pos.base_time(), 10*Secs);
		assert_eq!(pos.spd.base_time(), 0*Secs);
		assert_eq!(pos.spd.accel.base_time(), 20*Secs);
		assert_eq!(pos.spd.accel.jerk.base_time(), 10*Secs);
		
		 // Values:
		assert_eq!(pos.at(0*Secs).round(), 32.);
		assert_eq!(pos.at(10*Secs).round(), -63.);
		assert_eq!(pos.at(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
		
		 // Update:
		assert_eq!(pos.at_mut(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209778.);
		assert_eq!(pos.at_mut(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209778.);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_poly!(
			pos.poly(10*Secs),
			Poly::from(Sum::new(-63.15, [
				-11.9775,
				0.270416666666666,
				0.0475,
				-0.000416666666666,
			]))
		);
		for _ in 0..2 {
			pos.at_mut(20*Secs);
			assert_poly!(
				pos.poly(pos.base_time()),
				Poly::from(Sum::new(-112.55, [
					 6.0141666666666666666,
					 1.4454166666666666666,
					 0.0308333333333333333,
					-0.0004166666666666666,
				]))
			);
		}
		pos.at_mut(0*Secs);
		assert_poly!(
			pos.poly(pos.base_time()),
			Poly::from(Sum::new(32., [
				-1.4691666666666666666,
				-1.4045833333333333333,
				 0.0641666666666666666,
				-0.0004166666666666666,
			]))
		);
	}
	
	#[test]
	fn when() {
		let pos = position();
		let acc = Accel {
			value: 0.3,
			jerk: Jerk {
				value: 0.4,
				snap: Snap { value: -0.01 },
			},
		}.to_flux(Time::ZERO);
		
		assert_time_ranges!(pos.when(Ordering::Greater, &acc), [
			(Time::ZERO, Time::from_secs_f64(4.56)),
			(Time::from_secs_f64(26.912), Time::from_secs_f64(127.394))
		]);
		assert_times!(pos.when_eq(&acc), [
			Time::from_secs_f64(4.56),
			Time::from_secs_f64(26.912),
			Time::from_secs_f64(127.394)
		]);
	}
	
	#[test]
	fn when_at_mut() {
		let mut a_pos = position();
		let mut b_pos = position();
		
		 // Check Before:
		assert_time_ranges!(a_pos.when(Ordering::Greater, &b_pos), []);
		assert_time_ranges!(a_pos.when(Ordering::Equal, &b_pos), [(Time::ZERO, Time::MAX)]);
		assert_times!(a_pos.when_eq(&b_pos), []);
		a_pos.at_mut(20*Secs);
		
		 // Apply Changes:
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 12.,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.at_mut(10*Secs).value -= 100.0;
		
		 // Check After:
		assert_time_ranges!(a_pos.when(Ordering::Greater, &b_pos), [
			(0*Secs, 8*Secs),
			(50*Secs, Time::MAX)
		]);
		assert_times!(a_pos.when_eq(&b_pos), [8*Secs, 50*Secs]);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq)]
		#[flux(Sum<f64, 2> = {value} + spd.per(Mins), crate = "crate")]
		#[derive(Debug)]
		struct Pos {
			value: i64,
			spd: Spd,
		}
		
		#[derive(PartialOrd, PartialEq)]
		#[flux(Sum<f64, 1> = {value} + acc.per(Secs), crate = "crate")]
		#[derive(Debug)]
		struct Spd {
			value: i64,
			acc: Acc,
		}
		
		#[derive(PartialOrd, PartialEq)]
		#[flux(Constant<f64> = {value}, crate = "crate")]
		#[derive(Debug)]
		struct Acc {
			value: i64,
		}
		
		let a_pos = [
			Pos { value: 3, spd: Spd { value: 300, acc: Acc { value: -240 } } },
			Pos { value: -4, spd: Spd { value: 120, acc: Acc { value: 1080 } } }
		].to_flux(Time::ZERO);
		let b_pos = [
			Pos { value: 8, spd: Spd { value: 330, acc: Acc { value: -300 } } },
			Pos { value: 4, spd: Spd { value: 600, acc: Acc { value: 720 } } }
		].to_flux(Time::ZERO);
		
		let dis = Spd { value: 10, acc: Acc { value: 0 } }
			.to_flux(Time::ZERO);
		
		assert_time_ranges!(
			a_pos.when_dis(&b_pos, Ordering::Less, &dis),
			[
				(Time::ZERO, Time::from_secs_f64(0.0823337)),
				(Time::from_secs_f64(2.46704544), Time::from_secs_f64(4.116193987))
			]
		);
		
		let b_pos = b_pos.to_moment(Time::ZERO).to_flux(1*Secs);
		assert_times!(
			a_pos.when_dis_eq(&b_pos, &Constant::from(2.)),
			[Time::from_secs_f64(0.414068993), Time::from_secs_f64(0.84545191)]
		);
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}