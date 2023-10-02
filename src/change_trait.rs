//! ...

use std::slice::IterMut;
use std::vec::IntoIter;
use time::{Time, TimeUnit};
use crate::change::{Linear, LinearIso, Scalar};
use crate::degree::Degree;

pub trait FluxValue: Sized {
	type Iso: LinearIso;
	const DEGREE: Degree;
	
	/// Initial value.
	fn value(&self) -> <Self::Iso as LinearIso>::Linear; // ??? Unsure if this should return a reference or not
	
	/// Accumulate change over time.
	fn change(&self, changes: &mut ChangeAccum<Self::Iso>);
	
	/// Apply change over time.
	fn update(&mut self, time: Time); 
	
	fn calc(&self, time: Time) -> <Self::Iso as LinearIso>::Linear {
		let mut changes = ChangeAccum::new(time);
		self.change(&mut changes);
		match changes.sum {
			Some(change_sum) => self.value() + change_sum,
			None => self.value(),
		}
	}
	
	fn at(&self, time: Time) -> Self::Iso {
		Self::Iso::map(self.calc(time))
	}
}

/// Change accumulator.
#[derive(Debug)]
pub struct ChangeAccum<T: LinearIso> {
	sum: Option<T::Linear>,
	time: f64,
	depth: f64,
}

impl<T: LinearIso> ChangeAccum<T> {
	fn new(time: Time) -> Self {
		Self {
			sum: None,
			time: (time >> TimeUnit::Nanosecs) as f64,
			depth: 0.0,
		}
	}
	
	pub fn add<U: FluxValue<Iso=T>>(&mut self, value: &U, unit: TimeUnit) {
		self.accum(Scalar(1.0), value, unit);
	}
	
	pub fn sub<U: FluxValue<Iso=T>>(&mut self, value: &U, unit: TimeUnit) {
		self.accum(Scalar(-1.0), value, unit);
	}
	
	fn accum<U: FluxValue<Iso=T>>(&mut self, scalar: Scalar, change: &U, unit: TimeUnit) {
		let sum = std::mem::take(&mut self.sum);
		
		 // Fetch Value of Change:
		self.depth += 1.0;
		change.change(self);
		self.depth -= 1.0;
		let value = match self.sum {
			Some(change_sum) => change.value() + change_sum,
			None => change.value(),
		};
		
		 // Add to Sum:
		let time = self.time * ((unit >> TimeUnit::Nanosecs) as f64).recip();
		let addend = value * Scalar((time + self.depth) / (self.depth + 1.0)) * scalar;
		self.sum = match sum {
			Some(sum) => Some(sum + addend),
			None => Some(addend),
		};
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	
	#[derive(Debug)] struct Pos   { value: f64, spd: Spd, }
	#[derive(Debug)] struct Spd   { value: f64, fric: Fric, accel: Accel, }
	#[derive(Debug)] struct Fric  { value: f64 }
	#[derive(Debug)] struct Accel { value: f64, jerk: Jerk, }
	#[derive(Debug)] struct Jerk  { value: f64, snap: Snap, }
	#[derive(Debug)] struct Snap  { value: f64 }
	
	impl FluxValue for Pos {
		type Iso = Linear<i64>; // ??? Iso could be generic, allowing different types of the same isomorphism to be compared
		const DEGREE: Degree = 4;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {
			changes.add(&self.spd, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.spd.update(time);
		}
	}
	
	impl FluxValue for Spd {
		type Iso = Linear<i64>;
		const DEGREE: Degree = 3;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {
			changes.sub(&self.fric, TimeUnit::Secs);
			changes.add(&self.accel, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.fric.update(time);
			self.accel.update(time);
		}
	}
	
	impl FluxValue for Fric {
		type Iso = Linear<i64>;
		const DEGREE: Degree = 0;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	impl FluxValue for Accel {
		type Iso = Linear<i64>;
		const DEGREE: Degree = 2;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {
			changes.add(&self.jerk, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.jerk.update(time);
		}
	}
	
	impl FluxValue for Jerk {
		type Iso = Linear<i64>;
		const DEGREE: Degree = 1;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {
			changes.add(&self.snap, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.snap.update(time);
		}
	}
	
	impl FluxValue for Snap {
		type Iso = Linear<i64>;
		const DEGREE: Degree = 0;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self::Iso>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	#[test]
	fn test() {
		let mut position = Pos {
			value: 32.0, 
			spd: Spd {
				value: 0.0,
				fric: Fric {
					value: 3.5,
				},
				accel: Accel {
					value: 0.3,
					jerk: Jerk {
						value: 0.4,
						snap: Snap { value: -0.01 },
					},
				},
			},
		};
		
		 // Values:
		assert_eq!(position.at(0*Secs), Linear::from(32));
		assert_eq!(position.at(10*Secs), Linear::from(-63));
		assert_eq!(position.at(20*Secs), Linear::from(-113));
		assert_eq!(position.at(100*Secs), Linear::from(8339));
		assert_eq!(position.at(200*Secs), Linear::from(-209779));
		
		 // Update:
		position.update(20*Secs);
		assert_eq!(position.at(0*Secs), Linear::from(-113));
		assert_eq!(position.at(80*Secs), Linear::from(8339));
		assert_eq!(position.at(180*Secs), Linear::from(-209779));
		position.update(55*Secs);
		assert_eq!(position.at(25*Secs), Linear::from(8339));
		assert_eq!(position.at(125*Secs), Linear::from(-209779));
	}
}