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

/// Moment-in-time interface for [`Flux::at`] / [`Flux::at_mut`].
pub trait Moment {
	type Flux: Flux<Moment=Self>;
	
	/// Constructs the entirety of a [`Flux`] from a single moment.
	fn to_flux(self, time: Time) -> Self::Flux;
}

/// A value that can change over time.
pub trait Flux: Sized {
	type Moment: Moment<Flux=Self>;
	// !!! Might be worth making Flux generic over the Moment, so something like
	// a `Position` type can represent its x & y parts separately if desired.
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// The output accumulator of [`Flux::change`].
	type OutAccum<'a>: FluxAccum<'a, Self::Kind>;
	
	/// Initial value.
	fn value(&self) -> <Self::Kind as FluxKind>::Value;
	
	/// ...
	fn time(&self) -> Time;
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> Self::OutAccum<'a>;
	
	/// The moment of this value at the given time.
	fn at(&self, time: Time) -> Self::Moment;
	
	/// Sets the moment of this value at the given time (affects all moments).
	fn set_at(&mut self, time: Time, moment: Self::Moment) {
		*self = moment.to_flux(time);
	}
	
	/// Returns a mutable moment in time.
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
	fn at_mut(&mut self, time: Time) -> FluxMut<Self> {
		let moment = Some(self.at(time));
		FluxMut {
			inner: self,
			time,
			moment,
		}
	}
	
	/// The evaluation of this value at the given time.
	fn value_at(&self, time: Time) -> <Self::Kind as FluxKind>::Value {
		let mut value = self.value();
		let value_accum = FluxAccumKind::Value {
			value: &mut value,
			depth: 0,
			time,
			offset: self.time(),
		};
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(value_accum));
		value
	}
	
	/// A polynomial description of this flux value relative to `self.time()`.
	fn poly(&self, time: Time) -> Poly<Self::Kind> {
		let mut poly = Poly::default();
		if self.time() == time {
			poly.0 = self.value();
			let poly_accum = FluxAccumKind::Poly { poly: &mut poly, depth: 0 };
			self.change(<Self::Kind as FluxKind>::Accum::from_kind(poly_accum));
		} else {
			let new = self.at(time).to_flux(time); // !!! This may be lossy
			poly.0 = new.value();
			let poly_accum = FluxAccumKind::Poly { poly: &mut poly, depth: 0 };
			new.change(<Self::Kind as FluxKind>::Accum::from_kind(poly_accum));
		}
		poly
	}
	
	/// Ranges of time when this value will be above/below/equal to the given value.
	fn when<O: Flux>(&self, cmp_order: Ordering, other: &O) -> TimeRanges
	where
		When<Self::Kind, O::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let time = self.time().max(other.time());
		When {
			a_poly: self.poly(time),
			b_poly: other.poly(time),
			cmp_order,
			time,
		}.into_iter()
	}
	
	/// Times when this value will be equal to the given value.
	fn when_eq<O: Flux>(&self, other: &O) -> Times
	where
		WhenEq<Self::Kind, O::Kind>: IntoIterator<IntoIter=Times>
	{
		let time = self.time().max(other.time());
		WhenEq {
			a_poly: self.poly(time),
			b_poly: other.poly(time),
			time,
		}.into_iter()
	}
}

impl Moment for f64 {
	type Flux = Self;
	fn to_flux(self, _time: Time) -> Self::Flux {
		self
	}
}

impl Flux for f64 {
	type Moment = Self;
	type Kind = Constant<Self>;
	type OutAccum<'a> = ();
	fn value(&self) -> <Self::Kind as FluxKind>::Value {
		*self
	}
	fn time(&self) -> Time {
		Time::ZERO
	}
	fn change<'a>(&self, _accum: ()) -> Self::OutAccum<'a> {}
	fn at(&self, _time: Time) -> Self::Moment {
		*self
	}
}

/// Mutable moment-in-time interface for [`Flux::at_mut`].
pub struct FluxMut<'v, V: Flux> {
	inner: &'v mut V,
	time: Time,
	moment: Option<V::Moment>,
}

impl<V: Flux> Drop for FluxMut<'_, V> {
	fn drop(&mut self) {
		if let Some(moment) = std::mem::take(&mut self.moment) {
			self.inner.set_at(self.time, moment);
		}
	}
}

impl<V: Flux> Deref for FluxMut<'_, V> {
	type Target = V::Moment;
	fn deref(&self) -> &Self::Target {
		if let Some(moment) = self.moment.as_ref() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<V: Flux> DerefMut for FluxMut<'_, V> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		if let Some(moment) = self.moment.as_mut() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<V: Flux> Debug for FluxMut<'_, V>
where
	V::Moment: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<V::Moment as Debug>::fmt(self, f)
	}
}

impl<V: Flux> Display for FluxMut<'_, V>
where
	V::Moment: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<V::Moment as Display>::fmt(self, f)
	}
}

/// Multidimensional change over time.
pub trait FluxVec {
	type Kind: FluxKind;
	fn times(&self) -> Times;
	fn polys(&self, time: Time) -> Polys<Self::Kind>;
	
	/// Ranges of time when a vector is within a distance from another vector.
	fn when_dis<B: FluxVec + ?Sized, D: Flux>(
		&self,
		other: &B,
		cmp_order: Ordering,
		dis: &D,
	) -> TimeRanges
	where
		WhenDis<Self::Kind, B::Kind, D::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let time = std::iter::once(dis.time())
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
}

impl<T: Flux> FluxVec for [T] {
	type Kind = T::Kind;
	fn times(&self) -> Times {
		self.iter()
			.map(|x| x.time())
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
		(**self).times()
	}
	fn polys(&self, time: Time) -> Polys<Self::Kind> {
		(**self).polys(time)
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
// 		[self.0.time(), self.1.time()]
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
		poly.when(self.cmp_order, self.time)
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
		poly.when_eq(self.time)
	}
}

/// [`Flux::when_dis`] predictive distance comparison.
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
			Output = <<A as kind_ops::Sub<B>>::Output as kind_ops::Sqr>::Output
		>
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
		
		sum.when(self.cmp_order, self.time)
	}
}

/// Convenience for grouping the unmapped value & time in a [`Flux`] type.
#[derive(Copy, Clone, Debug, Default)]
pub struct FluxValue<T: Linear>(Time, T);

impl<T: Linear> FluxValue<T> {
	pub fn new(time: Time, value: T) -> Self {
		Self(time, value)
	}
	
	pub fn time(&self) -> Time {
		self.0
	}
}

impl<T: Linear> Deref for FluxValue<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.1
	}
}

impl<T: Linear> DerefMut for FluxValue<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.1
	}
}

/// Used to construct a [`Change`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(TimeUnit::Secs)` 
pub trait Per: Sized {
	fn per(&self, unit: TimeUnit) -> Change<Self> {
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
#[derive(Copy, Clone, Debug)]
pub struct Constant<T: Linear>(PhantomData<T>);

impl<T: Linear> Mul<Scalar> for Constant<T> {
	type Output = Self;
	fn mul(self, _rhs: Scalar) -> Self::Output {
		self
	}
}

impl<T: Linear> FluxKind for Constant<T> {
	type Value = T;
	type Accum<'a> = ();
	fn initial_order(&self) -> Option<Ordering> where T: PartialOrd {
		Some(Ordering::Equal)
	}
	fn zero() -> Self {
		Self(PhantomData)
	}
}

impl<T: Linear, K: FluxKind<Value=T>> Add<K> for Constant<T> {
	type Output = K;
	fn add(self, rhs: K) -> K {
		rhs
	}
}

impl<T: Linear> Mul<Constant<T>> for Constant<T> {
	type Output = Self;
	fn mul(self, _rhs: Self) -> Self {
		self
	}
}

impl<T: Linear> Mul<T> for Constant<T> {
	type Output = Self;
	fn mul(self, _rhs: T) -> Self {
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
				value: 0.0,
				fric: Fric { value: 3.5 },
				accel: Accel {
					value: 0.3,
					jerk: Jerk {
						value: 0.4,
						snap: Snap { value: -0.01 },
					},
				},
			},
			misc: Vec::new(),
		};
		pos.to_flux(Time::default())
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		
		 // Values:
		assert_eq!(pos.at(0*Secs).round(), 32.);
		assert_eq!(pos.at(10*Secs).round(), -63.);
		assert_eq!(pos.at(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
		
		 // Update:
		assert_eq!(pos.at_mut(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
		assert_eq!(pos.at_mut(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_eq!(
			pos.poly(pos.time()),
			Poly(32.0, Sum([
				-1.4691666666666667,
				-1.4045833333333335,
				0.06416666666666666,
				-0.00041666666666666664,
			]))
		);
		for _ in 0..2 {
			pos.at_mut(20*Secs);
			assert_eq!(
				pos.poly(pos.time()),
				Poly(-112.55000000000007, Sum([
					6.014166666666661,
					1.44541666666666666,
					0.030833333333333334,
					-0.00041666666666666664,
				]))
			);
		}
		pos.at_mut(0*Secs);
		assert_eq!(
			pos.poly(pos.time()),
			Poly(32.00000000000006, Sum([
				-1.4691666666666667,
				-1.4045833333333335,
				0.06416666666666666,
				-0.00041666666666666664,
			]))
		);
	}
	
	#[test]
	fn when() {
		let pos = position();
		let acc = position().spd.accel;
		
		let pos_when: Vec<(Time, Time)> = pos.when(Ordering::Greater, &acc).collect();
		let pos_when_eq: Vec<Time> = pos.when_eq(&acc).collect();
		
		assert_eq!(pos_when, [
			(0*Nanosecs, 4560099744*Nanosecs),
			(26912076291*Nanosecs, 127394131312*Nanosecs)
		]);
		assert_eq!(pos_when_eq, [
			4560099744*Nanosecs,
			26912076291*Nanosecs,
			127394131312*Nanosecs
		]);
		
		// pos.at_mut(20*Secs);
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
	fn when_at_mut() {
		let mut a_pos = position();
		let mut b_pos = position();
		
		 // Check Before:
		let vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		assert_eq!(vec, []);
		assert_eq!(vec_eq, [0*Nanosecs]);
		a_pos.at_mut(20*Secs);
		
		 // Apply Changes:
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 12.25,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.at_mut(10*Secs).value -= 100.0;
		
		 // Check After:
		let mut vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let mut vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		for (a, b) in &mut vec {
			*a = Time::from_secs(a.as_secs_f64().round() as u64);
			*b = Time::from_secs(b.as_secs_f64().round() as u64);
		}
		for t in &mut vec_eq {
			*t = Time::from_secs(t.as_secs_f64().round() as u64);
		}
		assert_eq!(vec, [(0*Secs, 8*Secs), (50*Secs, vec[1].1)]);
		assert_eq!(vec_eq[0..2], [8*Secs, 50*Secs]);
		assert!(vec[1].1 > 2_u64.pow(53)*Secs);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq)]
		#[flux(Sum<f64, 2> = {value} + spd.per(Secs), crate = "crate")]
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
		
		let a_pos = vec![
			Pos { value: 3, spd: Spd { value: 7, acc: Acc { value: -4 } } },
			Pos { value: -4, spd: Spd { value: -7, acc: Acc { value: 18 } } }
		].to_flux(Time::ZERO);
		let b_pos = vec![
			Pos { value: 8, spd: Spd { value: 8, acc: Acc { value: -5 } } },
			Pos { value: 4, spd: Spd { value: 4, acc: Acc { value: 12 } } }
		].to_flux(Time::ZERO);
		
		let dis = Spd { value: 10, acc: Acc { value: 0 } }
			.to_flux(Time::ZERO);
		
		assert_eq!(
			a_pos.when_dis(&b_pos, Ordering::Less, &dis)
				.into_iter()
				.collect::<Vec<(Time, Time)>>(),
			[
				(Time::ZERO, Time::from_secs_f64(0.0823337)),
				(Time::from_secs_f64(2.46704544), Time::from_secs_f64(4.116193987))
			]
		);
		
		let b_pos = b_pos.at(Time::ZERO).to_flux(1*Secs);
		assert_eq!(
			a_pos.when_dis(&*b_pos, Ordering::Less, &2.)
				.into_iter()
				.collect::<Vec<(Time, Time)>>(),
			[(Time::from_secs_f64(0.414068993), Time::from_secs_f64(0.84545191))]
		);
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}