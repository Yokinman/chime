use super::{FluxValue, Changes};
use crate::change::LinearValue;
use time::Time;

use std::vec::Vec;
use crate::degree::FluxKind;
// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<A: FluxValue> FluxValue for &A {
	type Value = A::Value;
	type Kind = A::Kind;
	type OutAccum<'a> = A::OutAccum<'a> where <Self::Kind as FluxKind>::Linear: 'a;
	fn value(&self) -> <Self::Kind as FluxKind>::Linear {
		A::value(self)
	}
	fn set_value(&mut self, _value: <Self::Kind as FluxKind>::Linear) {
		unreachable!() // ??? I think not
	}
	fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
		A::change(self, changes)
	}
	fn advance(&mut self, _time: Time) {
		unreachable!() // ??? Is it?
	}
}

impl<A: FluxValue> FluxValue for &mut A {
	type Value = A::Value;
	type Kind = A::Kind;
	type OutAccum<'a> = A::OutAccum<'a> where <Self::Kind as FluxKind>::Linear: 'a;
	fn value(&self) -> <Self::Kind as FluxKind>::Linear {
		A::value(self)
	}
	fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
		A::set_value(self, value)
	}
	fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
		A::change(self, changes)
	}
	fn advance(&mut self, time: Time) {
		A::advance(self, time)
	}
}

impl<T: FluxValue> FluxValue for Vec<T>
where
	for<'t> T::Kind: FluxKind<Accum<'t> = T::OutAccum<'t>>,
	for<'t> <T::Kind as FluxKind>::Linear: 't, // !!! I think this gives T::Linear a static lifetime, not good
{
	type Value = T::Value;
	type Kind = T::Kind;
	type OutAccum<'a> = T::OutAccum<'a> where <Self::Kind as FluxKind>::Linear: 'a;
	fn value(&self) -> <Self::Kind as FluxKind>::Linear {
		let mut value = <Self::Kind as FluxKind>::Linear::zero();
		for item in self {
			value = value + item.value();
		}
		value
	}
	fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
		for item in self {
			item.set_value(value); // ???
		}
	}
	fn change<'a>(&self, mut changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
		for item in self {
			changes = item.change(changes);
		}
		changes
	}
	fn advance(&mut self, time: Time) {
		for item in self {
			item.advance(time);
		}
	}
}