//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::Flux;
use crate::time::Time;
use crate::kind::{EmptyFluxAccum, FluxAccum, FluxKind, Poly};
use crate::linear::{Linear, Basis};

impl<'t, T> Flux for &'t T
where
	T: Flux,
{
	type Moment = T::Moment;
	type MomentMut<'a> = () where Self: 'a;
	type Kind = T::Kind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
	fn to_moment(&self, basis_time: Time, to_time: Time) -> Self::Moment {
		T::to_moment(self, basis_time, to_time)
	}
	fn to_moment_mut(&mut self, _basis_time: Time, _to_time: Time) -> Self::MomentMut<'_> {}
}

impl<T: Flux, const SIZE: usize> Flux for [T; SIZE] {
	type Moment = [T::Moment; SIZE];
	type MomentMut<'a> = [T::MomentMut<'a>; SIZE] where Self: 'a;
	type Kind = [T::Kind; SIZE];
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		self.each_ref().map(T::basis)
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
	fn to_moment(&self, basis_time: Time, to_time: Time) -> Self::Moment {
		self.each_ref().map(|x| x.to_moment(basis_time, to_time))
	}
	fn to_moment_mut(&mut self, basis_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.each_mut().map(|x| x.to_moment_mut(basis_time, to_time))
	}
}

// !!! impl<A: Flux, B: Flux> FluxVec for (A, B)

// !!! Make underlying linear type a Vec, instead of this just being a convenience.
impl<T: Flux> Flux for Vec<T> {
	type Moment = Vec<T::Moment>;
	type MomentMut<'a> = Vec<T::MomentMut<'a>> where Self: 'a;
	type Kind = T::Kind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		let mut value: <<T::Kind as FluxKind>::Basis as Basis>::Inner = Linear::zero();
		for item in self {
			value = value.add(item.basis());
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
	fn to_moment(&self, basis_time: Time, to_time: Time) -> Self::Moment {
		self.iter()
			.map(|x| x.to_moment(basis_time, to_time))
			.collect()
	}
	fn to_moment_mut(&mut self, basis_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|x| x.to_moment_mut(basis_time, to_time))
			.collect()
	}
}

impl<K, V> Flux for HashMap<K, V>
where
	K: std::hash::Hash + Eq + Clone,
	V: Flux,
{
	type Moment = HashMap<K, V::Moment>;
	type MomentMut<'a> = HashMap<&'a K, V::MomentMut<'a>> where Self: 'a;
	type Kind = V::Kind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		todo!()
	}
	fn change(&self, _accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		todo!()
	}
	fn to_moment(&self, basis_time: Time, to_time: Time) -> Self::Moment {
		self.iter()
			.map(|(k, v)| (k.clone(), v.to_moment(basis_time, to_time)))
			.collect()
	}
	fn to_moment_mut(&mut self, basis_time: Time, to_time: Time) -> Self::MomentMut<'_> {
		self.iter_mut()
			.map(|(k, v)| (k, v.to_moment_mut(basis_time, to_time)))
			.collect()
	}
}