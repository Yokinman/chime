//! Convenient [`FluxChange`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::{FluxChange, ToMoment, ToMomentMut};
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

impl<'t, T> ToMoment for &'t T
where
	T: ToMoment,
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

impl<'t, T> ToMoment for &'t mut T
where
	T: ToMoment,
{
	type Moment<'a> = T::Moment<'a> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		T::to_moment(self, from_time, to_time)
	}
}

impl<'t, T> ToMomentMut for &'t mut T
where
	T: ToMomentMut,
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

impl<T: ToMoment, const SIZE: usize> ToMoment for [T; SIZE] {
	type Moment<'a> = [T::Moment<'a>; SIZE] where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.each_ref().map(|x| x.to_moment(from_time, to_time))
	}
}

impl<T: ToMomentMut, const SIZE: usize> ToMomentMut for [T; SIZE] {
	type MomentMut<'a> = [T::MomentMut<'a>; SIZE] where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.each_mut().map(|x| x.to_moment_mut(from_time, to_time))
	}
}

impl<T: ToMoment> ToMoment for Vec<T> {
	type Moment<'a> = Vec<T::Moment<'a>> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.iter()
			.map(|x| x.to_moment(from_time, to_time))
			.collect()
	}
}

impl<T: ToMomentMut> ToMomentMut for Vec<T> {
	type MomentMut<'a> = Vec<T::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|x| x.to_moment_mut(from_time, to_time))
			.collect()
	}
}

impl<K, V> ToMoment for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: ToMoment,
{
	type Moment<'a> = HashMap<&'a K, V::Moment<'a>> where Self: 'a;
	fn to_moment(&self, from_time: Time, to_time: Time) -> Self::Moment<'_> {
		self.iter()
			.map(|(k, v)| (k, v.to_moment(from_time, to_time)))
			.collect()
	}
}

impl<K, V> ToMomentMut for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: ToMomentMut,
{
	type MomentMut<'a> = HashMap<&'a K, V::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, from_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|(k, v)| (k, v.to_moment_mut(from_time, to_time)))
			.collect()
	}
}