//! ...

use std::slice::IterMut;
use std::vec::IntoIter;
use time::{Time, TimeUnit};
use crate::change::{Linear, LinearIso, LinearValue, Scalar};
use crate::degree::{Degree, IsBelowDeg, IsDeg};
use crate::polynomial::Polynomial;

pub trait FluxValue: Sized {
	type Iso: LinearIso;
	type Degree: IsDeg;
	
	/// Initial value.
	fn value(&self) -> <Self::Iso as LinearIso>::Linear; // ??? Unsure if this should return a reference or not
	
	/// Accumulate change over time.
	fn change(&self, changes: &mut ChangeAccum<Self>);
	
	/// Apply change over time.
	fn update(&mut self, time: Time); 
	
	fn calc(&self, time: Time) -> <Self::Iso as LinearIso>::Linear {
		let mut changes = ChangeAccum::new(ChangeAccumArgs::Sum((time >> TimeUnit::Nanosecs) as f64));
		self.change(&mut changes);
		self.value() + changes.sum
	}
	
	fn at(&self, time: Time) -> Self::Iso {
		Self::Iso::map(self.calc(time))
	}
	
	fn coeff_calc(&self, depth: f64, degree: f64) -> <Self::Iso as LinearIso>::Linear {
		if depth == 0.0 && degree == 0.0 {
			self.value()
		} else {
			let mut accum = ChangeAccum::new(ChangeAccumArgs::Coeff(degree));
			accum.depth = depth;
			self.change(&mut accum);
			let coeff_sum = accum.sum * Scalar(1.0 / (depth + 1.0));
			if depth >= degree {
				self.value() + coeff_sum
			} else {
				coeff_sum
			}
		}
		// https://www.desmos.com/calculator/tn0udn7swy
	}
	
	fn polynomial(&self) -> Polynomial<<Self::Iso as LinearIso>::Linear> {
		match Self::Degree::USIZE {
			0 => Polynomial::Constant ([self.coeff_calc(0.0, 0.0)]),
			1 => Polynomial::Linear   ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0)]),
			2 => Polynomial::Quadratic([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0)]),
			3 => Polynomial::Cubic    ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0), self.coeff_calc(0.0, 3.0)]),
			4 => Polynomial::Quartic  ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0), self.coeff_calc(0.0, 3.0), self.coeff_calc(0.0, 4.0)]),
			_ => unreachable!()
		}
	}
}

#[derive(Debug)]
enum ChangeAccumArgs {
	Sum(f64),
	Coeff(f64),
}

/// Change accumulator.
#[derive(Debug)]
pub struct ChangeAccum<T: FluxValue> {
	sum: <T::Iso as LinearIso>::Linear,
	args: ChangeAccumArgs,
	depth: f64,
}

impl<T: FluxValue> ChangeAccum<T> {
	fn new(args: ChangeAccumArgs) -> Self {
		Self {
			sum: <<T::Iso as LinearIso>::Linear as LinearValue>::ZERO,
			args,
			depth: 0.0,
		}
	}
	
	pub fn add<U: FluxValue<Iso=T::Iso>>(&mut self, value: &U, unit: TimeUnit)
	where
		U::Degree: IsBelowDeg<T::Degree>
	{
		self.accum(Scalar(1.0), value, unit);
	}
	
	pub fn sub<U: FluxValue<Iso=T::Iso>>(&mut self, value: &U, unit: TimeUnit)
	where
		U::Degree: IsBelowDeg<T::Degree>
	{
		self.accum(Scalar(-1.0), value, unit);
	}
	
	fn accum<U: FluxValue<Iso=T::Iso>>(&mut self, scalar: Scalar, change: &U, unit: TimeUnit)
	where
		U::Degree: IsBelowDeg<T::Degree>
	{
		self.sum = self.sum + match &mut self.args {
			ChangeAccumArgs::Sum(time) => {
				 // Fetch Value of Change:
				let mut changes = ChangeAccum::new(ChangeAccumArgs::Sum(*time));
				changes.depth = self.depth + 1.0;
				change.change(&mut changes);
				let value = change.value() + changes.sum;
				
				let time = *time * ((unit >> TimeUnit::Nanosecs) as f64).recip();
				value * Scalar((time + self.depth) / (self.depth + 1.0)) * scalar
			},
			ChangeAccumArgs::Coeff(degree) => {
				let coeff = change.coeff_calc(self.depth + 1.0, *degree);
				if self.depth < *degree {
					let coeff = coeff * Scalar(((unit >> TimeUnit::Nanosecs) as f64).recip());
					if self.depth == 0.0 {
						coeff * scalar
					} else {
						(coeff + (change.coeff_calc(self.depth + 1.0, *degree + 1.0) * Scalar(self.depth))) * scalar
					}
				} else {
					coeff * Scalar(self.depth) * scalar
				}
			}
		};
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	use crate::degree::Deg;
	
	#[derive(Debug)] struct Pos   { value: f64, spd: Spd, }
	#[derive(Debug)] struct Spd   { value: f64, fric: Fric, accel: Accel, }
	#[derive(Debug)] struct Fric  { value: f64 }
	#[derive(Debug)] struct Accel { value: f64, jerk: Jerk, }
	#[derive(Debug)] struct Jerk  { value: f64, snap: Snap, }
	#[derive(Debug)] struct Snap  { value: f64 }
	
	impl FluxValue for Pos {
		type Iso = Linear<i64>; // ??? Iso could be generic, allowing different types of the same isomorphism to be compared
		type Degree = Deg<4>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {
			changes.add(&self.spd, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.spd.update(time);
		}
	}
	
	impl FluxValue for Spd {
		type Iso = Linear<i64>;
		type Degree = Deg<3>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {
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
		type Degree = Deg<0>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	impl FluxValue for Accel {
		type Iso = Linear<i64>;
		type Degree = Deg<2>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {
			changes.add(&self.jerk, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.jerk.update(time);
		}
	}
	
	impl FluxValue for Jerk {
		type Iso = Linear<i64>;
		type Degree = Deg<1>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {
			changes.add(&self.snap, TimeUnit::Secs);
		}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
			self.snap.update(time);
		}
	}
	
	impl FluxValue for Snap {
		type Iso = Linear<i64>;
		type Degree = Deg<0>;
		fn value(&self) -> <Self::Iso as LinearIso>::Linear {
			self.value
		}
		fn change(&self, changes: &mut ChangeAccum<Self>) {}
		fn update(&mut self, time: Time) {
			self.value = self.calc(time);
		}
	}
	
	fn pos() -> Pos {
		Pos {
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
		}
	}
	
	#[test]
	fn value() {
		let mut position = pos();
		
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
	
	#[test]
	fn poly() {
		let mut position = pos();
		assert_eq!(
			position.polynomial(),
			Polynomial::Quartic([32.0, -1.469166666666667e-9, -1.4045833333333335e-18, 0.06416666666666669e-27, -0.0004166666666666668e-36])
		);
		position.update(20*Secs);
		assert_eq!(
			position.polynomial(),
			Polynomial::Quartic([-112.55000000000007, 6.0141666666666615e-9, 1.4454166666666668e-18, 0.030833333333333334e-27, -0.0004166666666666668e-36])
		);
	}
}