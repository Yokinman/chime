//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::array;
use std::vec::Vec;

use crate::{Flux, FluxValue, FluxVec, Moment, MomentVec};
use crate::_hidden::InnerFlux;
use crate::time::Time;
use crate::kind::{FluxKind, FluxKindVec, Poly, PolyVec};
use crate::linear::{Linear, LinearPlus, Vector};

impl<T: Moment, const SIZE: usize> MomentVec<SIZE> for [T; SIZE] {
	type Flux = [T::Flux; SIZE];
	fn to_flux_vec(self, time: Time) -> Self::Flux {
		self.map(|x| x.to_flux(time))
	}
}

impl<T: Flux, const SIZE: usize> FluxVec<SIZE> for [T; SIZE] {
	type Moment = [T::Moment; SIZE];
	type Kind = [T::Kind; SIZE];
	fn index_base_time(&self, index: usize) -> Time {
		self[index].base_time()
	}
	fn index_poly(&self, index: usize, time: Time) -> Poly<<Self::Kind as FluxKindVec<SIZE>>::Output> {
		self[index].poly(time)
	}
	fn poly_vec(&self, time: Time) -> PolyVec<SIZE, Self::Kind> {
		PolyVec::new(array::from_fn(|i| self.index_poly(i, time).into_inner()), time)
	}
	fn to_moment_vec(self, time: Time) -> Self::Moment {
		self.map(|x| x.to_moment(time))
	}
}

impl<T: Moment, const SIZE: usize> MomentVec<SIZE> for T
where
	<<T::Flux as Flux>::Kind as FluxKind>::Value: Vector<SIZE, Output: LinearPlus>,
	<T::Flux as Flux>::Kind: FluxKindVec<SIZE>,
{
	type Flux = T::Flux;
	fn to_flux_vec(self, time: Time) -> Self::Flux {
		self.to_flux(time)
	}
}

impl<T: Flux, const SIZE: usize> FluxVec<SIZE> for T
where
	<T::Kind as FluxKind>::Value: Vector<SIZE, Output: LinearPlus>,
	T::Kind: FluxKindVec<SIZE>,
	T::Moment: MomentVec<SIZE, Flux=Self>,
{
	type Moment = T::Moment;
	type Kind = T::Kind;
	fn index_base_time(&self, _index: usize) -> Time {
		T::base_time(self)
	}
	fn index_poly(&self, index: usize, time: Time) -> Poly<<Self::Kind as FluxKindVec<SIZE>>::Output> {
		Poly::new(self.poly(time).into_inner().index(index), time)
	}
	fn poly_vec(&self, time: Time) -> PolyVec<SIZE, Self::Kind> {
		PolyVec::new(self.poly(time).into_inner(), time)
	}
	fn to_moment_vec(self, time: Time) -> Self::Moment {
		self.to_moment(time)
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
		let mut value = Linear::zero();
		for item in self {
			value = value + item.base_value(base_time).into_inner();
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