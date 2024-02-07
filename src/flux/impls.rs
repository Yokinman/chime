//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::vec::Vec;

use crate::{Constant, Flux, Moment};
use crate::time::Time;
use crate::kind::FluxKind;
use crate::linear::Linear;

// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<T: Linear> Flux for T {
	type Moment = Self;
	type Kind = Constant<Self>;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		*self
	}
	fn base_time(&self) -> Time {
		Time::ZERO
	}
	fn change<'a>(&self, _accum: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
	{}
	fn to_moment(&self, _time: Time) -> Self::Moment {
		*self
	}
}

impl<T: Linear> Moment for T {
	type Flux = Self;
	fn to_flux(&self, _time: Time) -> Self::Flux {
		*self
	}
}

impl<T: Flux> Flux for Vec<T>
where
	T::Kind: for<'a> FluxKind<Accum<'a> = <T::Kind as FluxKind>::OutAccum<'a>>,
	Vec<T::Moment>: Moment<Flux = Vec<T>>,
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	fn base_value(&self) -> <Self::Kind as FluxKind>::Value {
		let mut value = <Self::Kind as FluxKind>::Value::zero();
		let time = self.base_time();
		for item in self {
			value = value + item.value(time);
		}
		value
	}
	fn base_time(&self) -> Time {
		self.iter()
			.map(|x| x.base_time())
			.max()
			.unwrap_or_default()
	}
	fn change<'a>(&self, mut changes: <Self::Kind as FluxKind>::Accum<'a>)
		-> <Self::Kind as FluxKind>::OutAccum<'a>
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
	<T::Flux as Flux>::Kind:
		for<'a> FluxKind<Accum<'a> = <<T::Flux as Flux>::Kind as FluxKind>::OutAccum<'a>>,
{
	type Flux = Vec<T::Flux>;
	fn to_flux(&self, time: Time) -> Self::Flux {
		self.iter()
			.map(|x| x.to_flux(time))
			.collect()
	}
}