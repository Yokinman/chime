//! Convenient [`FluxChange`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::{FluxChange, FluxMoment, FluxMomentMut};
use crate::time::Time;
use crate::kind::{EmptyFluxAccum, FluxAccum, FluxKind, Poly};

impl<'t, T> FluxChange for &'t T
where
	T: FluxChange,
{
	type Kind = T::Kind;
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
}

impl<'t, T> FluxMoment for &'t T
where
	T: FluxMoment,
{
	type Moment<'a> = T::Moment<'a> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		T::to_moment(self, from_time, to_time)
	}
}

impl<'t, T> FluxChange for &'t mut T
where
	T: FluxChange,
{
	type Kind = T::Kind;
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
}

impl<'t, T> FluxMoment for &'t mut T
where
	T: FluxMoment,
{
	type Moment<'a> = T::Moment<'a> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		T::to_moment(self, from_time, to_time)
	}
}

impl<'t, T> FluxMomentMut for &'t mut T
where
	T: FluxMomentMut,
{
	type MomentMut<'a> = T::MomentMut<'a> where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		T::to_moment_mut(self, from_time, to_time)
	}
}

impl<T: FluxChange, const SIZE: usize> FluxChange for [T; SIZE] {
	type Kind = [T::Kind; SIZE];
	fn basis(&self) -> <Self::Kind as FluxKind>::Basis {
		crate::linear::BasisArray(self.each_ref().map(T::basis))
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		let FluxAccum { poly, time } = accum;
		let poly_time = poly.time;
		let mut accums = poly.into_iter()
			.map(|poly| FluxAccum { poly, time });
		let poly = self.each_ref()
			.map(|x| x.change(accums.next().unwrap()).poly.kind);
		FluxAccum {
			poly: Poly::new(poly, poly_time),
			time,
		}
	}
}

impl<T: FluxMoment, const SIZE: usize> FluxMoment for [T; SIZE] {
	type Moment<'a> = [T::Moment<'a>; SIZE] where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.each_ref().map(|x| x.to_moment(from_time, to_time))
	}
}

impl<T: FluxMomentMut, const SIZE: usize> FluxMomentMut for [T; SIZE] {
	type MomentMut<'a> = [T::MomentMut<'a>; SIZE] where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.each_mut().map(|x| x.to_moment_mut(from_time, to_time))
	}
}

impl<T: FluxMoment> FluxMoment for Vec<T> {
	type Moment<'a> = Vec<T::Moment<'a>> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.iter()
			.map(|x| x.to_moment(from_time, to_time))
			.collect()
	}
}

impl<T: FluxMomentMut> FluxMomentMut for Vec<T> {
	type MomentMut<'a> = Vec<T::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|x| x.to_moment_mut(from_time, to_time))
			.collect()
	}
}

impl<K, V> FluxMoment for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: FluxMoment,
{
	type Moment<'a> = HashMap<&'a K, V::Moment<'a>> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.iter()
			.map(|(k, v)| (k, v.to_moment(from_time, to_time)))
			.collect()
	}
}

impl<K, V> FluxMomentMut for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: FluxMomentMut,
{
	type MomentMut<'a> = HashMap<&'a K, V::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|(k, v)| (k, v.to_moment_mut(from_time, to_time)))
			.collect()
	}
}