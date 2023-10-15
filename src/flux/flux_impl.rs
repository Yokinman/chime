use super::{FluxValue, Changes};
use crate::change::LinearValue;
use time::Time;

use std::vec::Vec;
// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

impl<A: FluxValue> FluxValue for &A {
	type Value = A::Value;
	type Linear = A::Linear;
	type Degree = A::Degree;
	fn value(&self) -> Self::Linear {
		A::value(self)
	}
	fn set_value(&mut self, _value: Self::Linear) {
		unreachable!() // ??? I think not
	}
	fn change(&self, changes: &mut Changes<Self>) {
		A::change(self, changes)
	}
	fn update(&mut self, _time: Time) {
		unreachable!() // ??? Is it?
	}
}

impl<A: FluxValue> FluxValue for &mut A {
	type Value = A::Value;
	type Linear = A::Linear;
	type Degree = A::Degree;
	fn value(&self) -> Self::Linear {
		A::value(self)
	}
	fn set_value(&mut self, value: Self::Linear) {
		A::set_value(self, value);
	}
	fn change(&self, changes: &mut Changes<Self>) {
		A::change(self, changes)
	}
	fn update(&mut self, time: Time) {
		A::update(self, time)
	}
}

impl<T: FluxValue> FluxValue for Vec<T> {
	type Value = T::Value;
	type Linear = T::Linear;
	type Degree = T::Degree;
	fn value(&self) -> Self::Linear {
		let mut value = Self::Linear::zero();
		for item in self {
			value = value + item.value();
		}
		value
	}
	fn set_value(&mut self, value: Self::Linear) {
		for item in self {
			item.set_value(value); // ???
		}
	}
	fn change(&self, changes: &mut Changes<Self>) {
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