//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::Moment;
use crate::Flux;
use crate::time::Time;
use crate::kind::FluxKind;
use crate::linear::{Linear, LinearPlus};

impl<'t, T> Moment for &'t T
where
	T: Flux,
{
	type Flux = &'t T;
	fn to_flux(self, _time: Time) -> Self::Flux {
		unimplemented!()
	}
}

impl<'t, T> Flux for &'t T
where
	T: Flux,
{
	type Moment = &'t T;
	type Kind = T::Kind;
	fn base_value(&self, base_time: Time) -> <<Self::Kind as FluxKind>::Value as LinearPlus>::Inner {
		T::base_value(self, base_time)
	}
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> <Self::Kind as FluxKind>::OutAccum<'a> {
		T::change(self, accum)
	}
	fn to_moment(self, _base_time: Time, _time: Time) -> Self::Moment {
		unimplemented!("References to Flux types don't support `at` or `at_mut`")
	}
}


impl<T: Moment, const SIZE: usize> Moment for [T; SIZE] {
	type Flux = [T::Flux; SIZE];
	fn to_flux(self, time: Time) -> Self::Flux {
		self.map(|x| x.to_flux(time))
	}
}

impl<T: Flux, const SIZE: usize> Flux for [T; SIZE] {
	type Moment = [T::Moment; SIZE];
	type Kind = [T::Kind; SIZE];
	fn base_value(&self, base_time: Time)
		-> <<Self::Kind as FluxKind>::Value as LinearPlus>::Inner
	{
		self.each_ref().map(|x| x.base_value(base_time))
	}
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{
		let mut i = 0;
		accum.map(|a| {
			i += 1;
			self[i - 1].change(a)
		})
	}
	fn to_moment(self, base_time: Time, time: Time) -> Self::Moment {
		self.map(|x| x.to_moment(base_time, time))
	}
}

// !!! impl<A: Flux, B: Flux> FluxVec for (A, B)

// !!! Remove this impl, replace with a FluxVec impl:
impl<T: Flux> Flux for Vec<T>
where
	T::Kind: for<'a> FluxKind<Accum<'a> = <T::Kind as FluxKind>::OutAccum<'a>>,
	Vec<T::Moment>: Moment<Flux = Self>,
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	fn base_value(&self, base_time: Time)
		-> <<Self::Kind as FluxKind>::Value as LinearPlus>::Inner
	{
		let mut value: <<T::Kind as FluxKind>::Value as LinearPlus>::Inner = Linear::zero();
		for item in self {
			value = value.add(item.base_value(base_time));
		}
		value
	}
	fn change<'a>(&self, mut changes: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{
		for item in self {
			changes = item.change(changes);
		}
		changes
	}
	fn to_moment(self, base_time: Time, time: Time) -> Self::Moment {
		self.into_iter()
			.map(|x| x.to_moment(base_time, time))
			.collect()
	}
}

impl<T: Moment> Moment for Vec<T>
where
	<T::Flux as Flux>::Kind:
		for<'a> FluxKind<Accum<'a> = <<T::Flux as Flux>::Kind as FluxKind>::OutAccum<'a>>,
{
	type Flux = Vec<T::Flux>;
	fn to_flux(self, time: Time) -> Self::Flux {
		self.into_iter()
			.map(|x| x.to_flux(time))
			.collect()
	}
}


impl<K, V> Flux for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: Flux,
{
	type Moment = HashMap<K, V::Moment>;
	type Kind = V::Kind;
	fn base_value(&self, _base_time: Time) -> <<Self::Kind as FluxKind>::Value as LinearPlus>::Inner {
		todo!()
	}
	fn change<'a>(&self, _accum: <Self::Kind as FluxKind>::Accum<'a>) -> <Self::Kind as FluxKind>::OutAccum<'a> {
		todo!()
	}
	fn to_moment(self, base_time: Time, time: Time) -> Self::Moment {
		self.into_iter()
			.map(|(k, v)| (k, v.to_moment(base_time, time)))
			.collect()
	}
}

impl<K, V> Moment for HashMap<K, V>
where
	K: std::hash::Hash + Eq,
	V: Moment,
{
	type Flux = HashMap<K, V::Flux>;
	fn to_flux(self, time: Time) -> Self::Flux {
		self.into_iter()
			.map(|(k, v)| (k, v.to_flux(time)))
			.collect()
	}
}