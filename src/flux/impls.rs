//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::vec::Vec;

use crate::{Flux, Moment};
use crate::time::Time;
use crate::kind::FluxKind;
use crate::linear::Linear;

// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<T: Flux, const S: usize> Flux for [T; S]
where
	T::Kind: for<'a> FluxKind<Accum<'a> = T::OutAccum<'a>>
{
	type Moment = [T::Moment; S];
	type Kind = T::Kind;
	type OutAccum<'a> = T::OutAccum<'a>;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		let mut value = <Self::Kind as FluxKind>::Value::zero();
		let time = self.base_time();
		for item in self {
			value = value + item.value(time);
		}
		value
	}
	fn base_time(&self) -> Time {
		self.first()
			.map(|x| x.base_time())
			.unwrap_or_default()
	}
	fn change<'a>(&self, mut changes: <Self::Kind as FluxKind>::Accum<'a>)
		-> Self::OutAccum<'a>
	{
		for item in self {
			changes = item.change(changes);
		}
		changes
	}
	fn to_moment(&self, time: Time) -> Self::Moment {
		std::array::from_fn(|i| self[i].to_moment(time))
	}
}

impl<T: Moment, const S: usize> Moment for [T; S]
where
	[T::Flux; S]: Flux<Moment = [T; S]>
{
	type Flux = [T::Flux; S];
	fn to_flux(&self, time: Time) -> Self::Flux {
		std::array::from_fn(|i| self[i].to_flux(time))
	}
}

impl<T: Flux> Flux for Vec<T>
where
	T::Kind: for<'a> FluxKind<Accum<'a> = T::OutAccum<'a>>
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	type OutAccum<'a> = T::OutAccum<'a>;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		let mut value = <Self::Kind as FluxKind>::Value::zero();
		let time = self.base_time();
		for item in self {
			value = value + item.value(time);
		}
		value
	}
	fn base_time(&self) -> Time {
		self.first()
			.map(|x| x.base_time())
			.unwrap_or_default()
	}
	fn change<'a>(&self, mut changes: <Self::Kind as FluxKind>::Accum<'a>)
		-> Self::OutAccum<'a>
	{
		for item in self {
			changes = item.change(changes);
		}
		changes
	}
	fn to_moment(&self, time: Time) -> Self::Moment {
		self.iter()
			.map(|x| x.to_moment(time))
			.collect()
	}
}

impl<T: Moment> Moment for Vec<T>
where
	Vec<T::Flux>: Flux<Moment = Vec<T>>
{
	type Flux = Vec<T::Flux>;
	fn to_flux(&self, time: Time) -> Self::Flux {
		self.iter()
			.map(|x| x.to_flux(time))
			.collect()
	}
}