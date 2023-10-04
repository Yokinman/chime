use super::{FluxValue, ChangeAccum};
use crate::change::{LinearIso, LinearValue};
use time::Time;

use std::vec::Vec;
use std::cell::{Cell, RefCell};
// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

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

impl<T: FluxValue> FluxValue for Cell<T>
where
	T: Copy + Clone
{
	type Iso = T::Iso;
	type Degree = T::Degree;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		self.get().value()
	}
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
		self.get().change(changes)
	}
	fn update(&mut self, time: Time) {
		self.get_mut().update(time)
	}
}

impl<T: FluxValue> FluxValue for RefCell<T> {
	type Iso = T::Iso;
	type Degree = T::Degree;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		self.borrow().value()
	}
	fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
		self.borrow().change(changes)
	}
	fn update(&mut self, time: Time) {
		self.borrow_mut().update(time)
	}
}