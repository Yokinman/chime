//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::array;
use std::vec::Vec;

use crate::{Flux, FluxValue, FluxVec, Moment, MomentVec};
use crate::_hidden::InnerFlux;
use crate::time::Time;
use crate::kind::{ArrayFluxKindValue, FluxKind, Poly, PolyVec};
use crate::linear::{Linear, LinearPlus, Vector};

impl<T: Moment, const SIZE: usize> Moment for [T; SIZE] {
	type Flux = FluxValue<[<T::Flux as Flux>::Inner; SIZE]>;
	fn to_flux(self, time: Time) -> Self::Flux {
		FluxValue::new(self.map(|x| x.to_flux(time).into_inner_flux()), time)
	}
}

impl<T: InnerFlux, const SIZE: usize> InnerFlux for [T; SIZE] {
	type Moment = [T::Moment; SIZE];
	type Kind = [T::Kind; SIZE];
	fn base_value(&self, base_time: Time) -> <Self::Kind as FluxKind>::Value {
		ArrayFluxKindValue(self.each_ref().map(|x| x.base_value(base_time)))
	}
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> <Self::Kind as FluxKind>::OutAccum<'a> {
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
impl<T: InnerFlux> InnerFlux for Vec<T>
where
	T::Kind: for<'a> FluxKind<Accum<'a> = <T::Kind as FluxKind>::OutAccum<'a>>,
	Vec<T::Moment>: Moment<Flux = FluxValue<Vec<T>>>,
{
	type Moment = Vec<T::Moment>;
	type Kind = T::Kind;
	fn base_value(&self, base_time: Time) -> <Self::Kind as FluxKind>::Value {
		let mut value: <<T::Kind as FluxKind>::Value as LinearPlus>::Inner = Linear::zero();
		for item in self {
			value = value.add(item.base_value(base_time).into_inner());
		}
		LinearPlus::from_inner(value)
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
	type Flux = FluxValue<Vec<<T::Flux as Flux>::Inner>>;
	fn to_flux(self, time: Time) -> Self::Flux {
		FluxValue::new(
			self.into_iter()
				.map(|x| x.to_flux(time).into_inner_flux())
				.collect(),
			time,
		)
	}
}