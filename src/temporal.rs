//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Sub};

use crate::linear::{Basis, Linear, Vector};
use crate::time::Time;
use crate::{Flux, Moment, MomentMut, ToMoment, ToMomentMut};
use crate::change::Change;
use crate::poly::{Constant, Poly, ops as kind_ops, Roots, Deriv, Translate};
use crate::pred::{When, WhenDis, WhenDisEq, WhenEq};

/// A [`Poly`] paired with a basis time.
/// e.g. `Temporal<1 + 2x>` => `1 + 2(x-time)`.
#[derive(Copy, Clone, Debug, Default)]
pub struct Temporal<T> {
	pub(crate) inner: T,
	pub(crate) time: Time,
}

/// ... [`<Temporal as IntoIterator>::IntoIter`]
pub struct TemporalIter<T> {
	iter: T,
	time: Time,
}

impl<K> PartialEq for Temporal<K>
where
	K: PartialEq
{
	fn eq(&self, other: &Self) -> bool {
		// ??? Deriving PartialEq, Eq could count `f(t) = 1 + 2t` and
 		// `g(t) = 3 + 2(t-basis_time)` as the same if `basis_time = 1`.
		self.inner == other.inner && self.time == other.time
	}
}
	
impl<T> Deref for Temporal<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.inner
	}
}

impl<T> DerefMut for Temporal<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.inner
	}
}

impl<T> Temporal<T> {
	pub fn new(inner: T, time: Time) -> Self {
		Self { inner, time }
	}
	
	pub fn map_ref<U>(&self, f: impl FnOnce(&T) -> &U) -> Temporal<&U> {
		Temporal {
			inner: f(&self.inner),
			time: self.time,
		}
	}
	
	pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Temporal<U> {
		Temporal {
			inner: f(self.inner),
			time: self.time,
		}
	}
	
	pub fn as_ref(&self) -> Temporal<&T> {
		Temporal {
			inner: &self.inner,
			time: self.time,
		}
	}
	
	fn secs(&self, time: Time) -> f64 {
		if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		}
	}
}

impl<T: Flux> Temporal<T> {
	/// An evaluation of this flux at some point in time.
	pub fn basis(&self) -> T::Basis {
		self.inner.basis()
	}
	
	/// The time of [`Flux::basis`].
	pub fn basis_time(&self) -> Time {
		self.time
	}
	
	/// Conversion into a standard representation.
	pub fn to_poly(&self) -> Temporal<<T::Change as Change<T::Basis>>::Poly> {
		Temporal {
			inner: self.inner.to_poly(),
			time: self.time,
		}
	}
}

impl<T> Temporal<T> {
	/// Ranges when this is above/below/equal to another flux.
	pub fn when<U, B>(self, cmp: Ordering, other: Temporal<U>) -> T::Pred
	where
		T: When<U::Output, B>,
		U: Translate<B>,
		B: Basis,
	{
		let time = self.time;
		T::when(self, cmp, other.at_time(time))
	}
	
	/// Times when this is equal to another flux.
	pub fn when_eq<U, B>(self, other: Temporal<U>) -> T::Pred
	where
		T: WhenEq<U::Output, B>,
		U: Translate<B>,
		B: Basis,
	{
		let time = self.time;
		T::when_eq(self, other.at_time(time))
	}
	
	/// Ranges when this is above/below/equal to a constant.
	pub fn when_constant<B>(self, cmp: Ordering, other: B) -> T::Pred
	where
		T: When<symb_poly::Invar<Constant<B>>, B>,
		B: Basis,
	{
		self.when(cmp, Temporal::from(symb_poly::Invar(Constant(other))))
	}
	
	/// Times when this is equal to a constant.
	pub fn when_eq_constant<B>(self, other: B) -> T::Pred
	where
		T: WhenEq<symb_poly::Invar<Constant<B>>, B>,
		B: Basis,
	{
		self.when_eq(Temporal::from(symb_poly::Invar(Constant(other))))
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	pub fn when_dis<U, D, A, B, const SIZE: usize>(
		self,
		other: Temporal<U>,
		cmp: Ordering,
		dis: Temporal<D>,
	) -> T::Pred
	where
		T: WhenDis<U::Output, D::Output, B, SIZE>,
		U: Translate<A>,
		D: Translate<B>,
		A: Basis,
		B: Basis,
	{
		let time = self.time;
		T::when_dis(self, other.at_time(time), cmp, dis.at_time(time))
	}
	
	/// Ranges when the distance to another vector is equal to X.
	pub fn when_dis_eq<U, D, A, B, const SIZE: usize>(
		self,
		other: Temporal<U>,
		dis: Temporal<D>,
	) -> T::Pred
	where
		T: WhenDisEq<U::Output, D::Output, B, SIZE>,
		U: Translate<A>,
		D: Translate<B>,
		A: Basis,
		B: Basis,
	{
		let time = self.time;
		T::when_dis_eq(self, other.at_time(time), dis.at_time(time))
	}
	
	/// Ranges when the distance to another vector is above/below/equal to a constant.
	pub fn when_dis_constant<U, B, const SIZE: usize>(
		self,
		other: Temporal<U>,
		cmp: Ordering,
		dis: B,
	) -> T::Pred
	where
		T: Vector<SIZE, Output: Poly<B>>
			+ WhenDis<U::Output, symb_poly::Invar<Constant<B>>, B, SIZE>,
		U: Translate<B>,
		B: Basis,
	{
		self.when_dis(other, cmp, Temporal::from(symb_poly::Invar(Constant(dis))))
	}
	
	/// Ranges when the distance to another vector is equal to a constant.
	pub fn when_dis_eq_constant<U, B, const SIZE: usize>(
		self,
		other: Temporal<U>,
		dis: B,
	) -> T::Pred
	where
		T: Vector<SIZE, Output: Poly<B>>
			+ WhenDisEq<U::Output, symb_poly::Invar<Constant<B>>, B, SIZE>,
		U: Translate<B>,
		B: Basis,
	{
		self.when_dis_eq(other, Temporal::from(symb_poly::Invar(Constant(dis))))
	}
	
	/// Ranges when a component is above/below/equal to another flux.
	pub fn when_index<U, B, const SIZE: usize>(
		self,
		index: usize,
		cmp: Ordering,
		other: Temporal<U>,
	) -> <T::Output as When<U::Output, B>>::Pred
	where
		T: Vector<SIZE, Output: When<U::Output, B>>,
		U: Translate<B>,
		B: Basis,
	{
		let time = self.time;
		<T::Output as When<U::Output, B>>::when(self.index(index), cmp, other.at_time(time))
	}
	
	/// Times when a component is equal to another flux.
	pub fn when_index_eq<U, B, const SIZE: usize>(
		self,
		index: usize,
		other: Temporal<U>,
	) -> <T::Output as WhenEq<U::Output, B>>::Pred
	where
		T: Vector<SIZE, Output: WhenEq<U::Output, B>>,
		U: Translate<B>,
		B: Basis,
	{
		let time = self.time;
		<T::Output as WhenEq<U::Output, B>>::when_eq(self.index(index), other.at_time(time))
	}
	
	/// Ranges when a component is above/below/equal to a constant.
	pub fn when_index_constant<B, const SIZE: usize>(
		self,
		index: usize,
		cmp: Ordering,
		other: B,
	) -> <T::Output as When<symb_poly::Invar<Constant<B>>, B>>::Pred
	where
		T: Vector<SIZE, Output: When<symb_poly::Invar<Constant<B>>, B>>,
		B: Basis,
	{
		self.when_index(index, cmp, Temporal::from(symb_poly::Invar(Constant(other))))
	}
	
	/// Times when a component is equal to a constant.
	pub fn when_index_eq_constant<B, const SIZE: usize>(
		self,
		index: usize,
		other: B,
	) -> <T::Output as WhenEq<symb_poly::Invar<Constant<B>>, B>>::Pred
	where
		T: Vector<SIZE, Output: WhenEq<symb_poly::Invar<Constant<B>>, B>>,
		B: Basis,
	{
		self.when_index_eq(index, Temporal::from(symb_poly::Invar(Constant(other))))
	}
	
	// !!!
	// - to rotate line by a fixed angle, multiply axial polynomials by Scalar?
	// - to fit a line segment, filter predicted times.
	// - for non-fixed angle line segments, use a double distance check?
	// - rotating point-line may be handleable iteratively, find the bounds in
	//   which the roots may be and iterate through it.
}

impl<K> Temporal<K> {
	pub fn eval<B>(&self, time: Time) -> B
	where
		K: Poly<B>,
		B: Basis,
	{
		self.inner.eval(Linear::from_f64(self.secs(time)))
	}
	
	pub fn initial_order<B>(self, time: Time) -> Option<Ordering>
	where
		K: Deriv<B>,
		B: Basis + PartialOrd,
	{
		let secs = self.secs(time);
		self.inner.initial_order(Linear::from_f64(secs))
	}
	
	pub fn deriv<B>(self) -> Temporal<K::Deriv>
	where
		K: Deriv<B>,
		B: Basis,
	{
		Temporal {
			inner: <K as Deriv<B>>::deriv(self.inner),
			time: self.time,
		}
	}
	
	pub fn at_time<B>(self, time: Time) -> Temporal<K::Output>
	where
		K: Translate<B>,
		B: Basis,
	{
		let secs = self.secs(time);
		Temporal {
			inner: self.inner.translate(Linear::from_f64(secs)),
			time,
		}
	}
	
	pub fn sqr(self) -> Temporal<K::Output>
	where
		K: kind_ops::Sqr
	{
		Temporal::new(self.inner.sqr(), self.time)
	}
	
	pub fn add_poly<P, B>(self, other: Temporal<P>) -> Temporal<K::Output>
	where
		K: Add<P::Output>,
		P: Translate<B>,
		B: Basis,
	{
		Temporal {
			inner: self.inner + other.at_time(self.time).inner,
			time: self.time,
		}
	}
	
	pub fn sub_poly<P, B>(self, other: Temporal<P>) -> Temporal<K::Output>
	where
		K: Sub<P::Output>,
		P: Translate<B>,
		B: Basis,
	{
		Temporal {
			inner: self.inner - other.at_time(self.time).inner,
			time: self.time,
		}
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	pub(crate) fn when_sign<F, B>(self, order: Ordering, filter: F) -> crate::pred::PredFilter<crate::pred::Pred<K, B>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots<B> + PartialEq,
		B: Basis + PartialOrd,
	{
		let pred = crate::pred::Pred {
			poly: self,
			order,
			basis: std::marker::PhantomData,
		};
		crate::pred::PredFilter { pred, filter }
	}
	
	/// Times when the value is equal to zero.
	pub(crate) fn when_zero<F, B>(self, filter: F) -> crate::pred::PredFilter<crate::pred::PredEq<K, B>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots<B> + PartialEq,
		B: Basis + PartialEq,
	{
		crate::pred::PredFilter {
			pred: crate::pred::PredEq { poly: self, basis: std::marker::PhantomData },
			filter,
		}
	}
}

impl<T: ToMoment> Temporal<T> {
	/// See [`ToMoment::to_moment`].
	pub fn to_moment(&self, time: Time) -> T::Moment<'_> {
		self.inner.to_moment(self.secs(time))
	}
	
	pub fn moment(&self, time: Time) -> Moment<T> {
		Moment {
			moment: self.to_moment(time),
			borrow: std::marker::PhantomData,
		}
	}
}

impl<T: ToMomentMut> Temporal<T> {
	/// See [`ToMomentMut::to_moment_mut`].
	pub fn to_moment_mut(&mut self, time: Time) -> T::MomentMut<'_> {
		let secs = self.secs(time);
		self.time = time;
		self.inner.to_moment_mut(secs)
	}
	
	pub fn moment_mut(&mut self, time: Time) -> MomentMut<T> {
		MomentMut {
			moment: self.to_moment_mut(time),
			borrow: std::marker::PhantomData,
		}
	}
}

impl<T> From<T> for Temporal<T> {
	fn from(value: T) -> Self {
		Self::new(value, Time::ZERO)
	}
}

impl<K, const SIZE: usize> Vector<SIZE> for Temporal<K>
where
	K: Vector<SIZE>,
{
	type Output = Temporal<K::Output>;
	fn index(&self, index: usize) -> Self::Output {
		Temporal {
			inner: self.inner.index(index),
			time: self.time,
		}
	}
}

impl<K> IntoIterator for Temporal<K>
where
	K: IntoIterator
{
	type Item = Temporal<K::Item>;
	type IntoIter = TemporalIter<K::IntoIter>;
	fn into_iter(self) -> Self::IntoIter {
		TemporalIter {
			iter: self.inner.into_iter(),
			time: self.time,
		}
	}
}

impl<T> Iterator for TemporalIter<T>
where
	T: Iterator
{
	type Item = Temporal<T::Item>;
	fn next(&mut self) -> Option<Self::Item> {
		self.iter.next()
			.map(|x| Temporal::new(x, self.time))
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.iter.size_hint()
	}
}

#[cfg(feature = "bevy")]
mod bevy_items {
	use bevy_ecs::component::{Component, ComponentHooks, StorageType};
	use bevy_ecs::system::Resource;
	use super::Temporal;
	
	/// Implemented by `T` to blanket impl [`Component`] for `Temporal<T>`.
	pub trait TemporalComponent: Send + Sync + 'static {}
	
	impl<T> Component for Temporal<T>
	where
		T: TemporalComponent
	{
		const STORAGE_TYPE: StorageType = StorageType::Table;
		
		fn register_component_hooks(hooks: &mut ComponentHooks) {
			hooks.on_add(|mut world, entity, _| {
				let time = world.get_resource::<bevy_time::Time<crate::Chime>>()
					.expect("`Time<Chime>` must exist in world")
					.elapsed();
				
				let mut temporal = world.get_mut::<Self>(entity)
					.expect("temporal entity should exist");
				
				temporal.time = time;
			});
		}
	}
	
	/// Implemented by `T` to blanket impl [`Resource`] for `Temporal<T>`.
	pub trait TemporalResource: Send + Sync + 'static {}
	
	impl<T> Resource for Temporal<T>
	where
		T: TemporalResource
	{}
}

#[cfg(feature = "bevy")]
pub use bevy_items::*;