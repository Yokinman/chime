//! Utilities for describing how a type changes over time.

pub mod time;
pub mod kind;
mod impls;

use std::array;
use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::ops::{Add, Deref, DerefMut, Mul};

use crate::{
	linear::*,
	kind::*,
};

use self::time::{Time, /*Times,*/ TimeRanges};

pub use self::impls::*;

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
	
	type Moment: Moment;
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// An evaluation of this flux at some point in time.
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value;
	
	/// The time of [`Flux::base_value`].
	fn base_time(&self) -> Time;
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>;
	
	/// A moment in the timeline.
	fn to_moment(&self, time: Time) -> Self::Moment;
	
	/// Sets a moment in the timeline (affects all moments).
	fn set_moment(&mut self, time: Time, moment: Self::Moment)
	where
		Self: Sized,
		<Self as Flux>::Moment: Moment<Flux=Self>,
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
	fn at_mut(&mut self, time: Time) -> MomentRefMut<Self::Moment>
	where
		<Self as Flux>::Moment: Moment<Flux=Self>
	{
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
		let base_time = self.base_time();
		if time == base_time {
			return self.base_value()
		}
		self.poly(base_time).at(time)
	}
	
	/// A polynomial description of this flux at the given time.
	fn poly(&self, time: Time) -> Poly<Self::Kind> {
		let mut poly = Self::Kind::from(self.value(time));
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(
			&mut poly,
			0,
			time,
		));
		Poly::new(poly, time)
	}
	
	/// Ranges when this is above/below/equal to another flux.
	fn when<T: Flux>(&self, order: Ordering, other: &T) -> TimeRanges
	where
		Poly<Self::Kind>: When<T::Kind>
	{
		let time = self.base_time();
		self.poly(time).when(order, other.poly(time))
	}
	
	/// Times when this is equal to another flux.
	fn when_eq<T: Flux>(&self, other: &T) -> TimeRanges
	where
		Poly<Self::Kind>: WhenEq<T::Kind>
	{
		let time = self.base_time();
		self.poly(time).when_eq(other.poly(time))
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
pub trait FluxVec<const SIZE: usize> {
	type Kind: FluxKind;
	fn times(&self) -> [Time; SIZE];
	fn polys(&self, time: Time) -> [Poly<Self::Kind>; SIZE];
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<T: FluxVec<SIZE> + ?Sized, D: Flux>(
		&self,
		other: &T,
		order: Ordering,
		dis: &D,
	) -> TimeRanges
	where
		[Poly<Self::Kind>; SIZE]: WhenDis<SIZE, T::Kind, D::Kind>
	{
		let time = self.times().into_iter()
			.chain(other.times())
			.max()
			.unwrap_or_default();
		
		self.polys(time)
			.when_dis(other.polys(time), order, dis.poly(time))
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis_eq<T: FluxVec<SIZE> + ?Sized, D: Flux>(
		&self,
		other: &T,
		dis: &D,
	) -> TimeRanges
	where
		[Poly<Self::Kind>; SIZE]: WhenDisEq<SIZE, T::Kind, D::Kind>
	{
		let time = self.times().into_iter()
			.chain(other.times())
			.max()
			.unwrap_or_default();
		
		self.polys(time)
			.when_dis_eq(other.polys(time), dis.poly(time))
	}
	
	// !!!
	// - for distance to an axis-aligned line, use normal prediction on 1 axis.
	// - to rotate line by a fixed angle, multiply axial polynomials by Scalar?
	// - to fit a line segment, filter predicted times.
	// - for non-fixed angle line segments, use a double distance check?
	// - rotating point-line may be handleable iteratively, find the bounds in
	//   which the roots may be and iterate through it.
}

impl<T: Flux, const SIZE: usize> FluxVec<SIZE> for [T; SIZE] {
	type Kind = T::Kind;
	fn times(&self) -> [Time; SIZE] {
		array::from_fn(|i| self[i].base_time())
	}
	fn polys(&self, time: Time) -> [Poly<Self::Kind>; SIZE] {
		array::from_fn(|i| self[i].poly(time))
	}
}

// !!! impl<A: Flux, B: Flux> FluxVec for (A, B)

/// Used to construct a [`Change`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(time_unit::SEC)` 
pub trait Per: Sized {
	fn per(&self, unit: Time) -> Change<&Self> {
		Change {
			rate: self,
			unit
		}
	}
}

impl<T: Flux> Per for T {}

/// A description of a change over time for use with arithmetic operators.
#[derive(Copy, Clone, Debug, Default, Ord, PartialOrd, Eq, PartialEq)]
pub struct Change<T> {
	pub rate: T,
	pub unit: Time,
}

impl<T> Change<T> {
	pub fn as_ref(&self) -> Change<&T> {
		Change {
			rate: &self.rate,
			unit: self.unit,
		}
	}
}

impl<T: Moment> Moment for Change<T> {
	type Flux = Change<T::Flux>;
	fn to_flux(&self, time: Time) -> Self::Flux {
		Change {
			rate: self.rate.to_flux(time),
			unit: self.unit,
		}
	}
}

impl<T: Flux> Flux for Change<T> {
	type Moment = Change<T::Moment>;
	type Kind = T::Kind;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		self.rate.base_value()
	}
	fn base_time(&self) -> Time {
		self.rate.base_time()
	}
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{
		self.rate.change(accum)
	}
	fn to_moment(&self, time: Time) -> Self::Moment {
		Change {
			rate: self.rate.to_moment(time),
			unit: self.unit,
		}
	}
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
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		self.value
	}
	fn base_time(&self) -> Time {
		self.time
	}
	fn change<'a>(&self, _accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{}
	fn to_moment(&self, _time: Time) -> Self::Moment {
		self.value
	}
}

impl<T: Linear> FluxKind for Constant<T> {
	type Value = T;
	type Accum<'a> = ();
	type OutAccum<'a> = ();
	fn at(&self, _time: Scalar) -> Self::Value {
		self.value
	}
	fn rate_at(&self, _time: Scalar) -> Self::Value {
		T::zero()
	}
	fn to_time(self, _time: Scalar) -> Self {
		self
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
	use super::time::SEC;
	use crate::sum::Sum;
	
	#[flux(Sum<f64, 4> = {value} + spd.per(SEC) + misc.per(SEC), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Pos {
		value: f64,
		spd: Spd,
		misc: Vec<Spd>,
	}
	
	#[flux(Sum<f64, 3> = {value} - fric.per(SEC) + accel.per(SEC), crate = "crate")]
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
	
	#[flux(Sum<f64, 2> = {value} + jerk.per(SEC), crate = "crate")]
	#[derive(Clone, Debug, Default)]
	struct Accel {
		value: f64,
		jerk: Jerk,
	}
	
	#[flux(Sum<f64, 1> = {value} + snap.per(SEC), crate = "crate")]
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
		let value = pos.value(10*SEC);
		pos.value.set_moment(10*SEC, value);
		pos.spd.accel.at_mut(20*SEC);
		pos.spd.accel.jerk.at_mut(10*SEC);
		pos
	}
	
	macro_rules! assert_times {
		($times:expr, []) => {
			assert_time_ranges!($times, [(Time::ZERO, Time::MAX)])
		};
		($times:expr, $cmp_times:expr) => {
			assert_time_ranges!(
				$times,
				$cmp_times.into_iter()
					.map(|x| (x, x))
					.collect::<Vec<_>>()
			)
			// let times: Times = $times;
			// let cmp_times = Times::new($cmp_times);
			// assert_eq!(
			// 	times.clone().count(),
			// 	cmp_times.clone().count(),
			// 	"a: {:?}, b: {:?}",
			// 	times.collect::<Box<[_]>>(),
			// 	cmp_times.collect::<Box<[_]>>(),
			// );
			// for (a, b) in times.zip(cmp_times) {
			// 	let time = ((a.max(b) - b.min(a)).as_secs_f64() * 60.).floor();
			// 	assert_eq!(time, 0., "a: {:?}, b: {:?}", a, b);
			// }
		};
	}
	
	macro_rules! assert_time_ranges {
		($ranges:expr, $cmp_ranges:expr) => {{
			let ranges: TimeRanges = $ranges;
			let cmp_ranges = Vec::from_iter($cmp_ranges).into_iter();
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
		assert_eq!(pos.base_time(), 10*SEC);
		assert_eq!(pos.spd.base_time(), 0*SEC);
		assert_eq!(pos.spd.accel.base_time(), 20*SEC);
		assert_eq!(pos.spd.accel.jerk.base_time(), 10*SEC);
		
		 // Values:
		assert_eq!(pos.at(0*SEC).round(), 32.);
		assert_eq!(pos.at(10*SEC).round(), -63.);
		assert_eq!(pos.at(20*SEC).round(), -113.);
		assert_eq!(pos.at(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
		
		 // Update:
		assert_eq!(pos.at_mut(20*SEC).round(), -113.);
		assert_eq!(pos.at(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
		assert_eq!(pos.at_mut(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_poly!(
			pos.poly(10*SEC),
			Poly::new(Sum::new(-63.15, [
				-11.9775,
				0.270416666666666,
				0.0475,
				-0.000416666666666,
			]), 10*SEC)
		);
		for _ in 0..2 {
			pos.at_mut(20*SEC);
			assert_poly!(
				pos.poly(pos.base_time()),
				Poly::new(Sum::new(-112.55, [
					 6.0141666666666666666,
					 1.4454166666666666666,
					 0.0308333333333333333,
					-0.0004166666666666666,
				]), 20*SEC)
			);
		}
		pos.at_mut(0*SEC);
		assert_poly!(
			pos.poly(pos.base_time()),
			Poly::new(Sum::new(32., [
				-1.4691666666666666666,
				-1.4045833333333333333,
				 0.0641666666666666666,
				-0.0004166666666666666,
			]), 0*SEC)
		);
		
		assert_poly!(
			Poly::new(Sum::new(-0.8_f64, [-2.7, 3.4, 2.8, 0.3]), 4000*time::MILLISEC).to_time(320*time::MILLISEC),
			Poly::new(Sum::new(-29.341750272_f64, [26.2289216, -3.13568, -1.616, 0.3]), 320*time::MILLISEC)
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
		a_pos.at_mut(20*SEC);
		
		 // Apply Changes:
		b_pos.at_mut(0*SEC).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.at_mut(0*SEC).misc.push(Spd {
			value: 12.,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.at_mut(10*SEC).value -= 100.0;
		
		 // Check After:
		assert_eq!(a_pos.when(Ordering::Greater, &b_pos).collect::<Vec<_>>(), [
			(0*SEC, 8*SEC - time::NANOSEC),
			(50*SEC, Time::MAX)
		]);
		assert_times!(a_pos.when_eq(&b_pos), [8*SEC, 50*SEC]);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq)]
		#[flux(Sum<f64, 2> = {value} + spd.per(time::MINUTE), crate = "crate")]
		#[derive(Debug)]
		struct Pos {
			value: i64,
			spd: Spd,
		}
		
		#[derive(PartialOrd, PartialEq)]
		#[flux(Sum<f64, 1> = {value} + acc.per(SEC), crate = "crate")]
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
		
		let b_pos = b_pos.to_moment(Time::ZERO).to_flux(SEC);
		assert_times!(
			a_pos.when_dis_eq(&b_pos, &2.),
			[Time::from_secs_f64(0.414068993), Time::from_secs_f64(0.84545191)]
		);
		assert_time_ranges!(
			a_pos.when_dis(&b_pos, Ordering::Equal, &2.),
			[
				(Time::from_secs_f64(0.414068993), Time::from_secs_f64(0.414068993)),
				(Time::from_secs_f64(0.84545191), Time::from_secs_f64(0.84545191))
			]
		);
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}