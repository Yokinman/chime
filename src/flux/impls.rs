//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::{Flux, ToMoment, ToMomentMut};
use crate::kind::{EmptyFluxAccum, FluxAccum, FluxKind, Poly};
use crate::linear::Scalar;

impl<'t, T> Flux for &'t T
where
	T: Flux,
{
	type Kind = T::Kind;
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
}

impl<'t, T> ToMoment for &'t T
where
	T: ToMoment,
{
	type Moment<'a> = T::Moment<'a> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		T::to_moment(self, time)
	}
}

impl<'t, T> Flux for &'t mut T
where
	T: Flux,
{
	type Kind = T::Kind;
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
}

impl<'t, T> ToMoment for &'t mut T
where
	T: ToMoment,
{
	type Moment<'a> = T::Moment<'a> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		T::to_moment(self, time)
	}
}

impl<'t, T> ToMomentMut for &'t mut T
where
	T: ToMomentMut,
{
	type MomentMut<'a> = T::MomentMut<'a> where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		T::to_moment_mut(self, time)
	}
}

impl<T: Flux, const SIZE: usize> Flux for [T; SIZE] {
	type Kind = [T::Kind; SIZE];
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		crate::linear::BasisArray(self.each_ref().map(T::basis))
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		let mut accums = accum.0.into_iter();
		let poly = self.each_ref()
			.map(|x| x.change(FluxAccum(accums.next().unwrap())).0);
		FluxAccum(poly)
	}
}

impl<T: ToMoment, const SIZE: usize> ToMoment for [T; SIZE] {
	type Moment<'a> = [T::Moment<'a>; SIZE] where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		self.each_ref().map(|x| x.to_moment(time))
	}
}

impl<T: ToMomentMut, const SIZE: usize> ToMomentMut for [T; SIZE] {
	type MomentMut<'a> = [T::MomentMut<'a>; SIZE] where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		self.each_mut().map(|x| x.to_moment_mut(time))
	}
}

impl<T: ToMoment> ToMoment for Vec<T> {
	type Moment<'a> = Vec<T::Moment<'a>> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		self.iter()
			.map(|x| x.to_moment(time))
			.collect()
	}
}

impl<T: ToMomentMut> ToMomentMut for Vec<T> {
	type MomentMut<'a> = Vec<T::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|x| x.to_moment_mut(time))
			.collect()
	}
}

impl<K, V> ToMoment for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: ToMoment,
{
	type Moment<'a> = HashMap<&'a K, V::Moment<'a>> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		self.iter()
			.map(|(k, v)| (k, v.to_moment(time)))
			.collect()
	}
}

impl<K, V> ToMomentMut for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: ToMomentMut,
{
	type MomentMut<'a> = HashMap<&'a K, V::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|(k, v)| (k, v.to_moment_mut(time)))
			.collect()
	}
}