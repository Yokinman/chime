//! A dynamically-definable value of change-over-time.

use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub, Div, Shr};

use time::{Time, TimeUnit};
use crate::flux::*;
use crate::linear::*;

/// A value with a dynamically-defined change over time.
#[derive(Clone, Debug)]
pub struct Flux<V, K: FluxKind> {
	value: LinearFlux<K::Linear>,
	phantom: PhantomData<V>,
}

impl<V, I: LinearIso<V>> Flux<V, Deg<I, 0>> {
	pub fn new(initial_value: V) -> Self {
		Self {
			value: LinearFlux::new(I::inv_map(initial_value)),
			phantom: PhantomData,
		}
	}
}

impl<V, K: FluxKind> Flux<V, K> {
	pub fn add_change<T>(mut self, rate: impl Into<Flux<V, T>>, unit: TimeUnit)
		-> Flux::<V, <K as AppendFluxKind<T>>::Output>
	where
		T: FluxKind<Linear=K::Linear>,
		<K as AppendFluxKind<T>>::Output: FluxKind<Linear=K::Linear>,
		K: AppendFluxKind<T>,
	{
		let change = Change::new(Scalar(1.0), rate.into().value, unit);
		self.value.add_change(change);
		Flux {
			value: self.value,
			phantom: PhantomData,
		}
	}
	
	pub fn sub_change<T>(mut self, rate: impl Into<Flux<V, T>>, unit: TimeUnit)
		-> Flux::<V, <K as AppendFluxKind<T>>::Output>
	where
		T: FluxKind<Linear=K::Linear>,
		<K as AppendFluxKind<T>>::Output: FluxKind<Linear=K::Linear>,
		K: AppendFluxKind<T>,
	{
		let change = Change::new(Scalar(-1.0), rate.into().value, unit);
		self.value.add_change(change);
		Flux {
			value: self.value,
			phantom: PhantomData,
		}
	}
}

impl<V, I: LinearIso<V>> FluxValue for Flux<V, Deg<I, 0>> {
	type Value = V;
	type Kind = Deg<I, 0>;
	type OutAccum<'v> = <Self::Kind as FluxKind>::Accum<'v> where I: 'v;
	fn value(&self) -> <Self::Kind as FluxKind>::Linear {
		self.value.initial_value
	}
	fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
		self.value.initial_value = value;
	}
	fn change<'v>(&self, changes: Changes<'v, Self>) -> Self::OutAccum<'v> {
		changes
	}
	fn advance(&mut self, _time: Time) {}
}

macro_rules! impl_flux_value_for_value {
	($($degree:literal),+) => {$(
		impl<V, I: LinearIso<V>> FluxValue for Flux<V, Deg<I, $degree>> {
			type Value = V;
			type Kind = Deg<I, $degree>;
			type OutAccum<'v> = <Self::Kind as FluxKind>::Accum<'v> where I: 'v;
			fn value(&self) -> <Self::Kind as FluxKind>::Linear {
				self.value.initial_value
			}
			fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
				self.value.initial_value = value;
			}
			fn change<'v>(&self, mut changes: Changes<'v, Self>) -> Self::OutAccum<'v> {
				for change in self.value.change_list.clone() {
					changes = if change.scalar.0 == 1.0 {
						changes.add(
							&Flux::<V, Deg<I, { $degree - 1 }>> {
								value: change.rate,
								phantom: PhantomData,
							},
							change.unit
						)
					} else {
						changes.sub(
							&Flux::<V, Deg<I, { $degree - 1 }>> {
								value: change.rate,
								phantom: PhantomData,
							},
							change.unit
						)
					};
				}
				changes
			}
			fn advance(&mut self, time: Time) {
				self.value.initial_value = self.calc_value(time);
				for change in &mut self.value.change_list {
					let mut v = Flux::<V, Deg<I, { $degree - 1 }>> {
						value: change.rate.clone(),
						phantom: PhantomData,
					};
					v.advance(time);
					change.rate = v.value;
				}
				// ??? Could probably reuse most of the initial_value calculation for
				// updating each sub-change's initial_value.
			}
		}
	)+};
}
impl_flux_value_for_value!(1, 2, 3, 4, 5, 6, 7);

/// A value that linearly changes over time.
#[derive(Clone, Debug)]
struct LinearFlux<T: LinearValue> {
	initial_value: T,
	change_list: Vec<Change<T>>,
}

impl<T: LinearValue> LinearFlux<T> {
	fn new(initial_value: T) -> Self {
		Self {
			initial_value,
			change_list: vec![],
		}
	}
	
	fn add_change(&mut self, change: Change<T>) {
		// ??? Maybe use independent time offsets instead of a shared time
		// offset, so less has to be updated when a change is made to a value.
		// Downside - this may complicate value calculation & root finding.
		self.change_list.push(change);
	}
}

/// A linear change over time.
#[derive(Clone, Debug)]
struct Change<T: LinearValue> {
	scalar: Scalar,
	rate: LinearFlux<T>,
	unit: TimeUnit,
}

impl<T: LinearValue> Change<T> {
	pub fn new(scalar: Scalar, rate: LinearFlux<T>, unit: TimeUnit) -> Self {
		Self {
			scalar,
			rate,
			unit,
		}
	}
}

/// Produces the degree of the resulting change, à la: `max(A, B+1)`. 
pub trait AppendFluxKind<K: FluxKind>: FluxKind {
	type Output: FluxKind;
}

impl<A, B> AppendFluxKind<B> for A
where
	A: FluxKind,
	B: FluxKind + Shr<DegShift>,
	<B as Shr<DegShift>>::Output: FluxKind + Add<A>,
	<<B as Shr<DegShift>>::Output as Add<A>>::Output: FluxKind,
{
	type Output = <<B as Shr<DegShift>>::Output as Add<A>>::Output;
}

/// Used to construct a [`FluxChange`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(TimeUnit::Secs)` 
pub trait Per: Sized {
	fn per(self, unit: TimeUnit) -> FluxChange<Self> {
		FluxChange {
			rate: self,
			unit
		}
	}
}
impl Per for f64 {}
impl Per for u64 {}
impl Per for i64 {}
impl<V, K: FluxKind> Per for Flux<V, K> {}

/// A change over time used with arithmetic operators to construct [`Flux`] values.
pub struct FluxChange<T> {
	rate: T,
	unit: TimeUnit,
}

macro_rules! impl_flux_ops {
	($op_trait:ident::$op_fn:ident, $change_fn:ident, $map:ty, $iso:ty) => {
		 // Num + Num.per(Unit)
		impl $op_trait<FluxChange<$map>> for $map
		where
			$iso: LinearIso<$map>
		{
			type Output = Flux<$map, Deg<$iso, 1>>;
			fn $op_fn(self, rhs: FluxChange<$map>) -> Self::Output {
				Flux::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Num + Flux.per(Unit)
		impl<K> $op_trait<FluxChange<Flux<$map, K>>> for $map
		where
			K: FluxKind<Linear=$iso>,
			Deg<$iso, 0>: AppendFluxKind<K>,
		    <Deg<$iso, 0> as AppendFluxKind<K>>::Output: FluxKind<Linear=K::Linear>,
			$iso: LinearIso<$map>,
		{
			type Output = Flux<$map, <Deg<$iso, 0> as AppendFluxKind<K>>::Output>;
			fn $op_fn(self, rhs: FluxChange<Flux<$map, K>>) -> Self::Output {
				Flux::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Flux + Num.per(Unit)
		impl<K> $op_trait<FluxChange<$map >> for Flux<$map, K>
		where
			K: FluxKind<Linear=$iso> + AppendFluxKind<Deg<$iso, 0>>,
		    <K as AppendFluxKind<Deg<$iso, 0>>>::Output: FluxKind<Linear=K::Linear>,
			$iso: LinearIso<$map>,
		{
			type Output = Flux<$map, <K as AppendFluxKind<Deg<$iso, 0>>>::Output>;
			fn $op_fn(self, rhs: FluxChange<$map>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Flux + Flux.per(Unit)
		impl<K, T> $op_trait<FluxChange<Flux<$map, T>>> for Flux<$map, K>
		where
			K: FluxKind + AppendFluxKind<T>,
			T: FluxKind<Linear=K::Linear>,
		    <K as AppendFluxKind<T>>::Output: FluxKind<Linear=K::Linear>,
			$iso: LinearIso<$map>,
		{
			type Output = Flux<$map, <K as AppendFluxKind<T>>::Output>;
			fn $op_fn(self, rhs: FluxChange<Flux<$map, T>>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.unit)
			}
		}
	}
}
impl_flux_ops!(Add::add, add_change, f64, f64);
impl_flux_ops!(Add::add, add_change, u64, f64);
impl_flux_ops!(Add::add, add_change, i64, f64);
impl_flux_ops!(Sub::sub, sub_change, f64, f64);
impl_flux_ops!(Sub::sub, sub_change, u64, f64);
impl_flux_ops!(Sub::sub, sub_change, i64, f64);
impl_flux_ops!(Mul::mul, add_change, f64, Exp<f64>);
impl_flux_ops!(Mul::mul, add_change, u64, Exp<f64>);
impl_flux_ops!(Div::div, sub_change, f64, Exp<f64>);
impl_flux_ops!(Div::div, sub_change, u64, Exp<f64>);

impl<V, I: LinearIso<V>> From<V> for Flux<V, Deg<I, 0>> {
	fn from(value: V) -> Self {
		Self::new(value)
	}
}

#[cfg(test)]
mod value_tests {
	use super::*;
	use TimeUnit::*;
	
	fn linear() -> (
		Flux<i64, Deg<f64, 1>>,
		Flux<i64, Deg<f64, 2>>,
		Flux<i64, Deg<f64, 3>>,
		Flux<i64, Deg<f64, 4>>,
	) {
		// !!! Change these polynomials to something more practical
		let v0 = 3_i64 - 5.per(Nanosecs);
		let v1 = 20 + v0.clone().per(Nanosecs) + v0.clone().per(Nanosecs);
		let v2 = 52 + v1.clone().per(Nanosecs);
		let v3 = 150 + v2.clone().per(Nanosecs) + v0.clone().per(Nanosecs);
		(v0, v1, v2, v3)
	}
	
	// #[test]
	// fn linear_polynomial() {
	// 	let v = linear();
	// 	assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([3.0, -5.0]));
	// 	assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([20.0, 6.0, -10.0]));
	// 	assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([52.0, 20.0, 6.0, -10.0]));
	// 	assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([150.0, 55.0, 15.0, 6.0, -10.0]));
	// }
	
	#[test]
	fn linear_value() {
		let v = linear();
		let n0 = [3, -2, -7, -12, -17, -22, -27, -32, -37, -42];
		let n1 = [20, 16, 2, -22, -56, -100, -154, -218, -292, -376];
		let n2 = [52, 68, 70, 48, -8, -108, -262, -480, -772, -1148];
		let n3 = [150, 216, 279, 315, 290, 160, -129, -641, -1450, -2640];
		for t in 0..10_u64 {
			assert_eq!(v.0.at(t*Nanosecs), n0[t as usize].into());
			assert_eq!(v.1.at(t*Nanosecs), n1[t as usize].into());
			assert_eq!(v.2.at(t*Nanosecs), n2[t as usize].into());
			assert_eq!(v.3.at(t*Nanosecs), n3[t as usize].into());
		}
	}
	
	#[test]
	fn calc_test() {
		let d = 1.0 + 2.0.per(Mins);
		let c = 1.0 + d.per(Mins) + 7.0.per(Secs);
		let c1 = 3.0 + 9.0.per(Millisecs);
		let b = 10.0 + c.per(Microsecs) + c1.per(Mins);
		let a = 30.0 + b.per(Microsecs);
		
		assert_eq!(a.poly(), crate::poly::Poly(
			30.0, [
			Deg(0.018166666666666664),
			Deg(3.668167918055556e-6),
			Deg(1.1764138889120372e-15),
			Deg(2.3148148148148152e-29),
		]));
	}
	
	#[test]
	fn linear_roots() {
		let (v0, v1, v2, v3) = linear();
		let (mut v0, mut v1, mut v2, mut v3) = (
			PredValue::new(v0),
			PredValue::new(v1),
			PredValue::new(v2),
			PredValue::new(v3),
		);
		assert_eq!(v0.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]); // 0.6
		assert_eq!(v1.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]); // -1.902, 2.102
		assert_eq!(v2.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [3*Nanosecs]); // 3.892
		assert_eq!(v3.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [5*Nanosecs]); // -3.687, 5.631
		v0.borrow_mut().advance(1*Nanosecs);
		v1.borrow_mut().advance(1*Nanosecs);
		v2.borrow_mut().advance(1*Nanosecs);
		v3.borrow_mut().advance(1*Nanosecs);
		assert_eq!(v0.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [1*Nanosecs]);
		assert_eq!(v2.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
		assert_eq!(v3.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [4*Nanosecs]);
		v0.borrow_mut().advance(2*Nanosecs);
		v1.borrow_mut().advance(2*Nanosecs);
		v2.borrow_mut().advance(2*Nanosecs);
		v3.borrow_mut().advance(2*Nanosecs);
		assert_eq!(v0.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]);
		assert_eq!(v3.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
	}
	
	fn exponential() -> (
		Flux<f64, Deg<Exp<f64>, 1>>,
		Flux<f64, Deg<Exp<f64>, 2>>,
		Flux<f64, Deg<Exp<f64>, 3>>,
		Flux<f64, Deg<Exp<f64>, 4>>,
	) {
		// !!! Change polynomial to something more practical
		let v0 = 3.2 / 1.1.per(Nanosecs);
		let v1 = 14.515 * v0.clone().per(Nanosecs) * v0.clone().per(Nanosecs);
		let v2 = 30.4 * v1.clone().per(Nanosecs);
		let v3 = 300.0 * v2.clone().per(Nanosecs) * v0.clone().per(Nanosecs);
		(v0, v1, v2, v3)
	}
	
	// #[test]
	// fn exponential_polynomial() {
	// 	let v = exponential();
	// 	fn ln(n: f64) -> f64 { n.ln() }
	// 	assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([ln(3.2), -ln(1.1)]));
	// 	assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// 	assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([ln(30.4), ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// 	assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([ln(300.0), ln(30.4) - ln(3.2), ln(14.515) + ln(1.1), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// }
	
	#[test]
	fn exponential_value() {
		let (v0, v1, v2, v3) = exponential();
		let n0 = [3.2, 2.9090909, 2.6446281, 2.4042074, 2.1856431];
		let n1 = [14.515, 122.83769, 859.13387, 4965.97682, 23722.647933];
		let n2 = [30.4, 3734.26565, 3208234.11489, 15932016253.16265, 377949612447409.6];
		for t in 0..5_u64 {
			assert!((LinearIso::<f64>::map(v0.at(t*Nanosecs)) - n0[t as usize]).abs() < 0.00001);
			assert!((LinearIso::<f64>::map(v1.at(t*Nanosecs)) - n1[t as usize]).abs() < 0.00001);
			assert!((LinearIso::<f64>::map(v2.at(t*Nanosecs)) - n2[t as usize]).abs() < 0.00001);
		}
		assert!(LinearIso::<f64>::map(v3.at(100*Nanosecs)).abs() < 0.00001);
	}
	
	#[test]
	fn exponential_roots() {
		let (v0, v1, v2, v3) = exponential();
		let (mut v0, mut v1, mut v2, mut v3) = (
			PredValue::new(v0),
			PredValue::new(v1),
			PredValue::new(v2),
			PredValue::new(v3),
		);
		assert_eq!(v0.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [12*Nanosecs]); // 12.204
		assert_eq!(v1.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [24*Nanosecs]); // -1.143, 24.551
		assert_eq!(v2.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [36*Nanosecs]); // 36.91 
		assert_eq!(v3.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [49*Nanosecs]); // -4.752, 49.329
		v0.borrow_mut().advance(6*Nanosecs);
		v1.borrow_mut().advance(6*Nanosecs);
		v2.borrow_mut().advance(6*Nanosecs);
		v3.borrow_mut().advance(6*Nanosecs);
		assert_eq!(v0.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [6*Nanosecs]);
		assert_eq!(v1.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [18*Nanosecs]);
		assert_eq!(v2.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [30*Nanosecs]);
		assert_eq!(v3.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [43*Nanosecs]);
		v0.borrow_mut().advance(20*Nanosecs);
		v1.borrow_mut().advance(20*Nanosecs);
		v2.borrow_mut().advance(20*Nanosecs);
		v3.borrow_mut().advance(20*Nanosecs);
		assert_eq!(v0.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [10*Nanosecs]);
		assert_eq!(v3.when_eq(&Flux::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]);
	}
	
	#[test]
	fn exponential_non_pos() {
		let v = 0.0 * 2.0.per(Nanosecs);
		
		// assert_eq!(v.polynomial(), FluxPolynomial::Linear([f64::NEG_INFINITY, f64::ln(2.0)]));
		assert_eq!(v.at(0*Nanosecs), 0.0);
		assert_eq!(v.at(u64::MAX*Nanosecs), 0.0);
		
		let v = -1.5 * 2.0.per(Nanosecs);
		
		// match v.polynomial() {
		// 	FluxPolynomial::Linear([a, _]) => assert!(a.is_nan()),
		// 	_ => panic!()
		// }
		assert!(v.at(0*Nanosecs).is_nan());
		assert!(v.at(u64::MAX*Nanosecs).is_nan());
	}
	
	#[test]
	fn long_term() {
		let h = 10.0 + 1.0.per(Hours);
		let n = 10.0 + 1.0.per(Nanosecs);
		assert_eq!(h.at(1*Hours), n.at(1*Nanosecs));
		assert_eq!(h.at(10*Hours), n.at(10*Nanosecs));
		
		let h1 = 30.0 + h.per(Hours);
		let n1 = 30.0 + n.per(Nanosecs);
		assert_eq!(h1.at(1*Hours), n1.at(1*Nanosecs));
		assert_eq!(h1.at(10*Hours), n1.at(10*Nanosecs));
	}
	
	#[test]
	fn discrete_update() {
		let (mut v0, _, _, mut v3) = linear();
		
		v0.advance(10*Nanosecs);
		let new_v0 = PredValue::new(v0.clone() + 7.per(Nanosecs));
		assert_eq!(v0.at(0*Nanosecs), new_v0.at(0*Nanosecs));
		assert_eq!(new_v0.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]); // 23.5
		
		v3.advance(6*Nanosecs);
		let new_v3 = PredValue::new(v3.clone() + 1000.per(Nanosecs));
		assert_eq!(v3.at(0*Nanosecs), new_v3.at(0*Nanosecs));
		assert_eq!(new_v3.when_eq(&Flux::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs, 3*Nanosecs]); // 0.22, 3.684
	}
	
	#[test]
	fn when() {
		use std::cmp::Ordering::*;
		let (v0, v1, v2, v3) = linear();
		assert_eq!(PredValue::new(v0.clone()).when(Equal,   &Flux::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v0.clone()).when(Less,    &Flux::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v0.clone()).when(Greater, &Flux::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(0*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Equal,   &Flux::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, 2*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Less,    &Flux::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Greater, &Flux::new(-25_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 3*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Equal,   &Flux::new(71_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs), (1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Less,    &Flux::new(20_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(3*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Greater, &Flux::new(69_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 2*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Equal,   &Flux::new(140_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(5*Nanosecs, 5*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Less,    &Flux::new(160_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 0*Nanosecs), (5*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Greater, &Flux::new(280_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(2*Nanosecs, 4*Nanosecs)]);
	}
	
	#[test]
	fn when_eq() {
		use std::cmp::Ordering::*;
		let (v0, v1, v2, v3) = linear();
		let (v0, v1, v2, v3) = (
			PredValue::new(v0),
			PredValue::new(v1),
			PredValue::new(v2),
			PredValue::new(v3),
		);
		assert_eq!(
			v0.when_eq(&Flux::new(3_i64)).into_iter().collect::<Vec<Time>>(),
			v0.when(Equal, &Flux::new(3_i64)).into_iter()
				.map(|(a, _)| a)
				.collect::<Vec<Time>>()
		);
		assert_eq!(
			v1.when_eq(&Flux::new(3_i64)).into_iter().collect::<Vec<Time>>(),
			v1.when(Equal, &Flux::new(3_i64)).into_iter()
				.map(|(a, _)| a)
				.collect::<Vec<Time>>()
		);
		assert_eq!(
			v2.when_eq(&Flux::new(3_i64)).into_iter().collect::<Vec<Time>>(),
			v2.when(Equal, &Flux::new(3_i64)).into_iter()
				.map(|(a, _)| a)
				.collect::<Vec<Time>>()
		);
		assert_eq!(
			v3.when_eq(&Flux::new(3_i64)).into_iter().collect::<Vec<Time>>(),
			v3.when(Equal, &Flux::new(3_i64)).into_iter()
				.map(|(a, _)| a)
				.collect::<Vec<Time>>()
		);
	}
}