//! Convenient [`Flux`] implementations (Vec<T>, ??? tuples, etc.).

use std::vec::Vec;

use crate::{Flux, Time, Moment};
use crate::kind::FluxKind;
use crate::linear::Linear;

// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<T: Flux> Flux for Vec<T>
where
	T::Kind: for<'a> FluxKind<Accum<'a> = T::OutAccum<'a>>
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	type OutAccum<'a> = T::OutAccum<'a>;
	fn value(&self) -> <Self::Kind as FluxKind>::Value {
		let mut value = <Self::Kind as FluxKind>::Value::zero();
		for item in self {
			value = value + item.value();
		}
		value
	}
	fn time(&self) -> Time {
		if let Some(item) = self.first() {
			item.time()
		} else {
			Time::default()
		}
	}
	fn change<'a>(&self, mut changes: <Self::Kind as FluxKind>::Accum<'a>)
		-> Self::OutAccum<'a>
	{
		for item in self {
			changes = item.change(changes);
		}
		changes
	}
	fn at(&self, time: Time) -> Self::Moment {
		self.iter()
			.map(|x| x.at(time))
			.collect()
	}
}

impl<T: Moment> Moment for Vec<T>
where
	Vec<T::Flux>: Flux<Moment=Vec<T>>
{
	type Flux = Vec<T::Flux>;
	fn to_flux(self, time: Time) -> Self::Flux {
		self.into_iter()
			.map(|x| x.to_flux(time))
			.collect()
	}
}