//! Convenient [`Flux`] implementations (Vec<T>, ??? tuples, etc.).

use super::*;

use std::vec::Vec;

use time::Time;
use crate::linear::Linear;

// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<T: Flux> Flux for Vec<T>
where
	for<'t> T::Kind: FluxKind<Accum<'t> = T::OutAccum<'t>>,
	for<'t> <T::Kind as FluxKind>::Linear: 't, // !!! I think this gives T::Linear a static lifetime, not good
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	type OutAccum<'a> = T::OutAccum<'a> where <Self::Kind as FluxKind>::Linear: 'a;
	fn value(&self) -> <Self::Kind as FluxKind>::Linear {
		let mut value = <Self::Kind as FluxKind>::Linear::zero();
		for item in self {
			value = value + item.value();
		}
		value
	}
	fn time(&self) -> Time {
		self.iter()
			.map(|x| x.time())
			.min() // ??? or maximum
			.unwrap_or_default()
	}
	fn change<'a>(&self, mut changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
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