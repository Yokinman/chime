//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::vec::Vec;

use crate::{Constant, Flux, FluxVec, Moment, MomentVec};
use crate::time::Time;
use crate::kind::{FluxKind, FluxKindVec, Poly};
use crate::linear::{Linear, LinearVec};

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


impl<T: Moment, const SIZE: usize> MomentVec<SIZE> for [T; SIZE] {
	type Flux = [T::Flux; SIZE];
	fn to_flux_vec(&self, time: Time) -> Self::Flux {
		std::array::from_fn(|i| self[i].to_flux(time))
	}
}

impl<T: Flux, const SIZE: usize> FluxVec<SIZE> for [T; SIZE] {
	type Moment = [T::Moment; SIZE];
	type Kind = T::Kind;
	fn index_base_time(&self, index: usize) -> Time {
		self[index].base_time()
	}
	fn index_poly(&self, index: usize, time: Time) -> Poly<Self::Kind> {
		self[index].poly(time)
	}
	fn to_moment_vec(&self, time: Time) -> Self::Moment {
		std::array::from_fn(|i| self[i].to_moment(time))
	}
}

impl<T: Moment, const SIZE: usize> MomentVec<SIZE> for T
where
	<<T::Flux as Flux>::Kind as FluxKind>::Value: LinearVec<SIZE>,
	<T::Flux as Flux>::Kind: FluxKindVec<SIZE>,
{
	type Flux = T::Flux;
	fn to_flux_vec(&self, time: Time) -> Self::Flux {
		self.to_flux(time)
	}
}

impl<T: Flux, const SIZE: usize> FluxVec<SIZE> for T
where
	<T::Kind as FluxKind>::Value: LinearVec<SIZE>,
	T::Kind: FluxKindVec<SIZE>,
	T::Moment: MomentVec<SIZE>,
{
	type Moment = T::Moment;
	type Kind = <T::Kind as FluxKindVec<SIZE>>::Kind;
	fn index_base_time(&self, _index: usize) -> Time {
		T::base_time(self)
	}
	fn index_poly(&self, index: usize, time: Time) -> Poly<Self::Kind> {
		Poly::new(T::poly(self, time).index_kind(index), time)
	}
	fn to_moment_vec(&self, time: Time) -> Self::Moment {
		self.to_moment(time)
	}
}

// !!! impl<A: Flux, B: Flux> FluxVec for (A, B)

// !!! Removed this impl, replace with a FluxVec impl:
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