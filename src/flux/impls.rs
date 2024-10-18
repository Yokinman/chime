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
	type MomentMut<'a> = T::MomentMut<'a> where Self: 'a;
	type Kind = T::Kind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		T::basis(self)
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		T::change(self, accum)
	}
	fn to_moment(&self, _base_time: Time, _time: Time) -> Self::Moment {
		unimplemented!("References to Flux types don't support `at` or `at_mut`")
	}
	fn set_moment(&mut self, _moment: Self::Moment) {
		unimplemented!()
	}
	fn from_moment(_moment: Self::Moment) -> Self {
		unimplemented!()
	}
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
	fn to_moment(&self, basis_time: Time, time: Time) -> Self::Moment {
		self.each_ref().map(|x| x.to_moment(basis_time, time))
	}
	fn set_moment(&mut self, moment: Self::Moment) {
		for (x, moment) in self.iter_mut().zip(moment) {
			x.set_moment(moment);
		}
	}
	fn from_moment(moment: Self::Moment) -> Self {
		moment.map(T::from_moment)
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
	fn to_moment(&self, basis_time: Time, time: Time) -> Self::Moment {
		self.iter()
			.map(|x| x.to_moment(basis_time, time))
			.collect()
	}
	fn set_moment(&mut self, moment: Self::Moment) {
		debug_assert_eq!(self.len(), moment.len());
		for (x, moment) in self.iter_mut().zip(moment) {
			x.set_moment(moment);
		}
	}
	fn from_moment(moment: Self::Moment) -> Self {
		moment.into_iter()
			.map(T::from_moment)
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
	fn to_moment(&self, basis_time: Time, time: Time) -> Self::Moment {
		self.iter()
			.map(|(k, v)| (k.clone(), v.to_moment(basis_time, time)))
			.collect()
	}
	fn set_moment(&mut self, moment: Self::Moment) {
		debug_assert_eq!(self.len(), moment.len());
		for (k, v) in moment {
			self.get_mut(&k)
				.expect("key should exist")
				.set_moment(v);
		}
	}
	fn from_moment(moment: Self::Moment) -> Self {
		moment.into_iter()
			.map(|(k, v)| (k, V::from_moment(v)))
			.collect()
	}
}