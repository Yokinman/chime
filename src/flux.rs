//! Utilities for describing how a type changes over time.

pub mod time;
pub mod kind;
mod impls;

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut, Mul};

use crate::{
	linear::*,
	kind::*,
};

use self::time::{Time, TimeRangeIter, TimeRanges};

pub use self::impls::*;

pub use chime_flux_proc_macro::flux;

/// A discrete interface for a value that changes over time.
pub trait Moment {
	type Flux: Flux<Moment=Self>;
	type Value: LinearIso<<<Self::Flux as Flux>::Kind as FluxKind>::Value>;
	
	/// Constructs the entirety of a [`Flux`] from a single moment.
	fn to_flux(self, time: Time) -> Self::Flux;
}

#[allow(type_alias_bounds)]
pub type FluxOf<T: Moment> = <T as Moment>::Flux;

/// The continuous interface for a value that changes over time.
pub trait Flux {
	// !!! Deriving PartialEq, Eq should count `f(t) = 1 + 2t` and
	// `g(t) = 3 + 2(t-base_time)` as the same Flux if `base_time = 1`.
	
	type Moment: Moment<Flux=Self>;
	
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
	fn to_moment(self, time: Time) -> Self::Moment;
	
	/// Sets a moment in the timeline (affects all moments).
	fn set_moment(&mut self, time: Time, moment: Self::Moment)
	where
		Self: Sized
	{
		*self = moment.to_flux(time);
	}
	
	/// A reference to a moment in the timeline.
	fn at(&self, time: Time) -> MomentRef<Self::Moment>
	where
		Self: Clone
	{
		let moment = self.clone().to_moment(time);
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
		Self: Clone
	{
		let moment = Some(self.clone().to_moment(time));
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
	fn poly(&self, time: Time) -> Poly<Self::Kind, <Self::Moment as Moment>::Value> {
		let mut poly = Self::Kind::from_value(self.value(time));
		self.change(poly.as_accum(0, time));
		Poly::new(poly, time).with_iso()
	}
	
	/// Ranges when this is above/below/equal to another flux.
	fn when<T>(&self, order: Ordering, other: &T) -> TimeRanges<impl TimeRangeIter>
	where
		T: Flux,
		Poly<Self::Kind, <Self::Moment as Moment>::Value>:
			When<T::Kind>
	{
		let time = self.base_time();
		self.poly(time).when(order, other.poly(time))
	}
	
	/// Times when this is equal to another flux.
	fn when_eq<T>(&self, other: &T) -> TimeRanges<impl TimeRangeIter>
	where
		T: Flux,
		Poly<Self::Kind, <Self::Moment as Moment>::Value>:
			WhenEq<T::Kind>
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

/// Multidimensional interface for a vector that changes over time.
pub trait MomentVec<const SIZE: usize> {
	type Flux: FluxVec<SIZE, Moment=Self>;
	type Value: LinearIsoVec<SIZE, <<Self::Flux as FluxVec<SIZE>>::Kind as FluxKindVec<SIZE>>::Value>;
	
	/// Constructs the entirety of a [`FluxVec`] from a single moment.
	fn to_flux_vec(self, time: Time) -> Self::Flux;
}

/// Multidimensional change over time.
pub trait FluxVec<const SIZE: usize> {
	type Moment: MomentVec<SIZE, Flux=Self>;
	type Kind: FluxKindVec<SIZE>;
	
	fn index_base_time(&self, index: usize) -> Time;
	fn max_base_time(&self) -> Time {
		let mut time = Time::ZERO;
		for i in 0..SIZE {
			time = time.max(self.index_base_time(i));
		}
		time
	}
	
	fn index_poly(&self, index: usize, time: Time) -> Poly<
		<Self::Kind as FluxKindVec<SIZE>>::Kind,
		<<Self::Moment as MomentVec<SIZE>>::Value as LinearIsoVec<SIZE, <Self::Kind as FluxKindVec<SIZE>>::Value>>::Value
	>;
	
	fn poly_vec(&self, time: Time)
		-> PolyVec<SIZE, Self::Kind, <Self::Moment as MomentVec<SIZE>>::Value>;
	
	fn to_moment_vec(self, time: Time) -> Self::Moment;
	
	fn set_moment_vec(&mut self, time: Time, moment: Self::Moment)
	where
		Self: Sized
	{
		*self = moment.to_flux_vec(time);
	}
	
	fn at_vec(&self, time: Time) -> MomentVecRef<SIZE, Self::Moment>
	where
		Self: Clone
	{
		let moment = self.clone().to_moment_vec(time);
		MomentVecRef {
			moment,
			borrow: PhantomData,
		}
	}
	
	fn at_vec_mut(&mut self, time: Time) -> MomentVecRefMut<SIZE, Self::Moment>
	where
		Self: Clone
	{
		let moment = Some(self.clone().to_moment_vec(time));
		MomentVecRefMut {
			moment,
			time,
			borrow: self,
		}
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<T, D>(&self, other: &T, order: Ordering, dis: &D) -> TimeRanges<impl TimeRangeIter>
	where
		T: FluxVec<SIZE> + ?Sized,
		D: Flux,
		PolyVec<SIZE, Self::Kind, <Self::Moment as MomentVec<SIZE>>::Value>:
			WhenDis<SIZE, T::Kind, D::Kind>,
	{
		let time = self.max_base_time();
		self.poly_vec(time)
			.when_dis(other.poly_vec(time), order, dis.poly(time))
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis_eq<T, D>(&self, other: &T, dis: &D) -> TimeRanges<impl TimeRangeIter>
	where
		T: FluxVec<SIZE> + ?Sized,
		D: Flux,
		PolyVec<SIZE, Self::Kind, <Self::Moment as MomentVec<SIZE>>::Value>:
			WhenDisEq<SIZE, T::Kind, D::Kind>,
	{
		let time = self.max_base_time();
		self.poly_vec(time)
			.when_dis_eq(other.poly_vec(time), dis.poly(time))
	}
	
	/// Ranges when a component is above/below/equal to another flux.
	fn when_index<T>(&self, index: usize, order: Ordering, other: &T) -> TimeRanges<impl TimeRangeIter>
	where
		T: Flux,
		Poly<<Self::Kind as FluxKindVec<SIZE>>::Kind, <<Self::Moment as MomentVec<SIZE>>::Value as LinearIsoVec<SIZE, <Self::Kind as FluxKindVec<SIZE>>::Value>>::Value>:
			When<T::Kind>
	{
		let time = self.index_base_time(index);
		self.index_poly(index, time)
			.when(order, other.poly(time))
	}
	
	/// Times when a component is equal to another flux.
	fn when_index_eq<T>(&self, index: usize, other: &T) -> TimeRanges<impl TimeRangeIter>
	where
		T: Flux,
		Poly<<Self::Kind as FluxKindVec<SIZE>>::Kind, <<Self::Moment as MomentVec<SIZE>>::Value as LinearIsoVec<SIZE, <Self::Kind as FluxKindVec<SIZE>>::Value>>::Value>:
			WhenEq<T::Kind>
	{
		let time = self.index_base_time(index);
		self.index_poly(index, time)
			.when_eq(other.poly(time))
	}
	
	// !!!
	// - to rotate line by a fixed angle, multiply axial polynomials by Scalar?
	// - to fit a line segment, filter predicted times.
	// - for non-fixed angle line segments, use a double distance check?
	// - rotating point-line may be handleable iteratively, find the bounds in
	//   which the roots may be and iterate through it.
}

/// Immutable moment-in-time interface for [`FluxVec::at_vec`].
pub struct MomentVecRef<'b, const SIZE: usize, M: MomentVec<SIZE>> {
	moment: M,
	borrow: PhantomData<&'b M::Flux>,
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Deref for MomentVecRef<'_, SIZE, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		&self.moment
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Debug for MomentVecRef<'_, SIZE, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Display for MomentVecRef<'_, SIZE, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

/// Mutable moment-in-time interface for [`FluxVec::at_vec_mut`].
pub struct MomentVecRefMut<'b, const SIZE: usize, M: MomentVec<SIZE>> {
	moment: Option<M>,
	time: Time,
	borrow: &'b mut M::Flux,
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Drop for MomentVecRefMut<'_, SIZE, M> {
	fn drop(&mut self) {
		if let Some(moment) = std::mem::take(&mut self.moment) {
			self.borrow.set_moment_vec(self.time, moment);
		}
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Deref for MomentVecRefMut<'_, SIZE, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		if let Some(moment) = self.moment.as_ref() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> DerefMut for MomentVecRefMut<'_, SIZE, M> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		if let Some(moment) = self.moment.as_mut() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Debug for MomentVecRefMut<'_, SIZE, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<const SIZE: usize, M: MomentVec<SIZE>> Display for MomentVecRefMut<'_, SIZE, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

/// Used to construct a [`Change`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(time_unit::SEC)` 
pub trait Per: Sized {
	fn per(self, unit: Time) -> Change<Self> {
		Change {
			rate: self,
			unit
		}
	}
}

impl<T> Per for T {}

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
	type Value = T::Value;
	fn to_flux(self, time: Time) -> Self::Flux {
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
	fn to_moment(self, time: Time) -> Self::Moment {
		Change {
			rate: self.rate.to_moment(time),
			unit: self.unit,
		}
	}
}

/// Wrapper for partial [`Flux`] types created by the [`flux`] macro.
#[derive(Copy, Clone, Debug, Default)]
#[cfg_attr(feature = "bevy", derive(
	bevy_ecs::component::Component,
	bevy_ecs::system::Resource,
))]
pub struct FluxValue<T> {
	inner: T,
	time: Time,
}

impl<T> FluxValue<T> {
	pub fn new(inner: T, time: Time) -> Self {
		Self { inner, time }
	}
	pub fn time(&self) -> Time {
		self.time
	}
}

impl<T> Deref for FluxValue<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.inner
	}
}

impl<T> DerefMut for FluxValue<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.inner
	}
}

impl<T: _hidden::InnerFlux> Flux for FluxValue<T>
where
	T::Moment: Moment<Flux=Self>
{
	type Moment = T::Moment;
	type Kind = T::Kind;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		self.inner.base_value()
	}
	fn base_time(&self) -> Time {
		self.time
	}
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{
		self.inner.change(accum)
	}
	fn to_moment(self, time: Time) -> Self::Moment {
		let value = self.value(time);
		self.inner.to_moment(time, value)
	}
}

#[doc(hidden)]
pub mod _hidden {
	use super::*;
	
	/// Intermediary for the [`FluxValue`] generic [`Flux`] implementation.
	/// Implemented automatically by the [`flux`] macro, not manually.
	pub trait InnerFlux {
		type Moment: Moment;
		type Kind: FluxKind;
		fn base_value(&self) -> <Self::Kind as FluxKind>::Value;
		fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
			-> <Self::Kind as FluxKind>::OutAccum<'a>;
		fn to_moment(self, time: Time, base_value: <Self::Kind as FluxKind>::Value)
			-> Self::Moment;
	}
}

/// No change over time.
/// 
/// Equivalent "constant" flux kinds should implement both `Into<Constant<T>>`
/// and `From<Constant<T>>` (e.g. `Sum<T, 0>`).
#[derive(Copy, Clone, Debug, Default)]
pub struct Constant<T>(T);

impl<T> Deref for Constant<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}

impl<T> DerefMut for Constant<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<T: Linear> From<T> for Constant<T> {
	fn from(value: T) -> Self {
		Constant(value)
	}
}

impl<T: Linear> FluxKind for Constant<T> {
	type Value = T;
	type Accum<'a> = ();
	type OutAccum<'a> = ();
	fn from_value(value: Self::Value) -> Self {
		Constant(value)
	}
	fn as_accum(&mut self, _depth: usize, _time: Time) -> Self::Accum<'_> {}
	fn at(&self, _time: Scalar) -> Self::Value {
		self.0
	}
	fn rate_at(&self, _time: Scalar) -> Self::Value {
		T::zero()
	}
	fn to_time(self, _time: Scalar) -> Self {
		self
	}
	fn initial_order(&self, _time: Scalar) -> Option<Ordering> where T: PartialOrd {
		self.0.partial_cmp(&T::zero())
	}
	fn zero() -> Self {
		Self::from(T::zero())
	}
}

impl<const SIZE: usize, T: LinearVec<SIZE>> FluxKindVec<SIZE> for Constant<T> {
	type Kind = Constant<T::Value>;
	type Value = T;
	fn index_kind(&self, index: usize) -> Self::Kind {
		Constant(self.0.index(index))
	}
}

impl<T: Linear> Mul<Scalar> for Constant<T> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.0 * rhs;
		self
	}
}

impl<T> Mul<Constant<T>> for Constant<T>
where
	T: Mul<Output = T>
{
	type Output = Self;
	fn mul(mut self, rhs: Self) -> Self {
		self.0 = self.0 * rhs.0;
		self
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use super::time::SEC;
	use crate::sum::Sum;
	
	#[flux(
		kind = "Sum<f64, 4>",
		value = value,
		change = |c| c + spd.per(SEC) + misc.per(SEC),
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Pos {
		value: f64,
		spd: Spd,
		misc: Vec<Spd>,
	}
	
	#[flux(
		kind = "Sum<f64, 3>",
		value = value,
		change = |c| c - fric.per(SEC) + accel.per(SEC),
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Spd {
		value: f64,
		fric: Fric,
		accel: Accel,
	}
	
	#[flux(
		kind = "Sum<f64, 0>",
		value = value,
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Fric {
		value: f64,
	}
	
	#[flux(
		kind = "Sum<f64, 2>",
		value = value,
		change = |c| c + jerk.per(SEC),
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Accel {
		value: f64,
		jerk: Jerk,
	}
	
	#[flux(
		kind = "Sum<f64, 1>",
		value = value,
		change = |c| c + snap.per(SEC),
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Jerk {
		value: f64,
		snap: Snap,
	}
	
	#[flux(
		kind = "Constant<f64>",
		value = value,
		crate = "crate",
	)]
	#[derive(Clone, Debug, Default)]
	struct Snap {
		value: f64,
	}
	
	impl Deref for Pos {
		type Target = f64;
		fn deref(&self) -> &Self::Target {
			&self.value
		}
	}
	
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
		pos.value = pos.value(10*SEC);
		pos.time = 10*SEC;
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
			let ranges: TimeRanges<_> = $ranges;
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
			for (degree, coeff) in dif_poly.into_inner().into_iter().enumerate() {
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
			(50*SEC + time::NANOSEC, Time::MAX)
		]);
		assert_times!(a_pos.when_eq(&b_pos), [8*SEC, 50*SEC]);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = "Sum<f64, 2>",
			value = value,
			change = |c| c + spd.per(time::MINUTE),
			crate = "crate",
		)]
		#[derive(Debug)]
		struct Pos {
			value: i64,
			spd: Spd,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = "Sum<f64, 1>",
			value = value,
			change = |c| c + acc.per(SEC),
			crate = "crate",
		)]
		#[derive(Debug)]
		struct Spd {
			value: i64,
			acc: Acc,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = "Constant<f64>",
			value = value,
			crate = "crate",
		)]
		#[derive(Debug)]
		struct Acc {
			value: i64,
		}
		
		let a_pos = [
			Pos { value: 3, spd: Spd { value: 300, acc: Acc { value: -240 } } },
			Pos { value: -4, spd: Spd { value: 120, acc: Acc { value: 1080 } } }
		].to_flux_vec(Time::ZERO);
		let b_pos = [
			Pos { value: 8, spd: Spd { value: 330, acc: Acc { value: -300 } } },
			Pos { value: 4, spd: Spd { value: 600, acc: Acc { value: 720 } } }
		].to_flux_vec(Time::ZERO);
		
		let dis = Spd { value: 10, acc: Acc { value: 0 } }
			.to_flux(Time::ZERO);
		
		assert_time_ranges!(
			a_pos.when_dis(&b_pos, Ordering::Less, &dis),
			// https://www.desmos.com/calculator/spxyoloyx9
			[
				(Time::ZERO, Time::from_secs_f64(0.0823337)),
				(Time::from_secs_f64(2.246506069), Time::from_secs_f64(4.116193987)),
			]
		);
		
		let b_pos = b_pos.to_moment_vec(Time::ZERO).to_flux_vec(SEC);
		assert_time_ranges!(
			a_pos.when_dis_eq(&b_pos, &2.),
			[
				(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
				(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191)),
			]
		);
		assert_time_ranges!(
			a_pos.when_dis(&b_pos, Ordering::Equal, &2.),
			[
				(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
				(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191))
			]
		);
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}