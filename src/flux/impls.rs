//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::Moment;
use crate::Flux;
use crate::time::Time;
use crate::kind::{EmptyFluxAccum, FluxAccum, FluxKind, Poly};
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
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
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
	fn to_moment(self, base_time: Time, time: Time) -> Self::Moment {
		self.map(|x| x.to_moment(base_time, time))
	}
}

// !!! impl<A: Flux, B: Flux> FluxVec for (A, B)

// !!! Make underlying linear type a Vec, instead of this just being a convenience.
impl<T: Flux> Flux for Vec<T>
where
	// T::Kind: for<'a> FluxKind<Accum<'a> = <T::Kind as FluxKind>::OutAccum<'a>>,
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
	fn change(&self, _accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		// let mut accum = accum.into();
		// for item in self {
		// 	changes = item.change(changes);
		// }
		// changes
		todo!()
	}
	fn to_moment(self, base_time: Time, time: Time) -> Self::Moment {
		self.into_iter()
			.map(|x| x.to_moment(base_time, time))
			.collect()
	}
}

impl<T: Moment> Moment for Vec<T>
where
	// <T::Flux as Flux>::Kind:
	// 	for<'a> FluxKind<Accum<'a> = <<T::Flux as Flux>::Kind as FluxKind>::OutAccum<'a>>,
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
	fn change(&self, _accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
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