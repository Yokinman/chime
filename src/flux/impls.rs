//! Convenient [`FluxValue`] implementations (Vec<T>, ??? tuples, etc.).

use super::*;

use std::vec::Vec;

use time::Time;
use crate::linear::LinearValue;

// !!! Arrays
// !!! Slices
// ??? Other collections
// ??? Tuples
// ??? Function pointers

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
	fn time(&self) -> Time {
		self.iter()
			.map(|x| x.time())
			.min()
			.unwrap_or_default()
	}
	fn set_time(&mut self, time: Time) {
		for item in self {
			item.set_time(time);
		}
	}
}