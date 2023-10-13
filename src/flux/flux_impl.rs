use super::{FluxValue, ChangeAccum};
use crate::change::{LinearIso, LinearValue};
use time::Time;

use std::vec::Vec;
// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<A: FluxValue> FluxValue for &A {
	type Iso = A::Iso;
	type Degree = A::Degree;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		A::value(self)
	}
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
		A::change(self, changes)
	}
	fn update(&mut self, _time: Time) {
		unreachable!() // ??? Is it?
	}
}

impl<A: FluxValue> FluxValue for &mut A {
	type Iso = A::Iso;
	type Degree = A::Degree;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		A::value(self)
	}
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
		A::change(self, changes)
	}
	fn update(&mut self, time: Time) {
		A::update(self, time)
	}
}

impl<T: FluxValue> FluxValue for Vec<T> {
	type Iso = T::Iso;
	type Degree = T::Degree;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		let mut value = <Self::Iso as LinearIso>::Linear::ZERO;
		for item in self {
			value = value + item.value();
		}
		value
	}
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
		for item in self {
			item.change(changes);
		}
	}
	fn update(&mut self, time: Time) {
		for item in self {
			item.update(time);
		}
	}
}