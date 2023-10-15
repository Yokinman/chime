//! Change describing utilities ...

use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub, Div};

use impl_op::impl_op;

use time::{Time, TimeUnit};
use crate::flux::{Changes, FluxValue};

use crate::degree::{Deg, IsDeg, IsBelowDeg, HasUpDeg, MaxDeg};
use crate::polynomial::{Poly, RootList, Roots};

/// A value with a dynamically-defined change over time.
#[derive(Clone, Debug)]
pub struct Value<V, I: LinearIso<V>, D: IsDeg> {
	value: LinearFluxValue<I>,
	degree: PhantomData<D>,
	mapped: PhantomData<V>,
}

impl<V, I: LinearIso<V>> Value<V, I, Deg<0>> {
	pub fn new(initial_value: V) -> Self {
		Self {
			value: LinearFluxValue::new(I::inv_map(initial_value)),
			degree: PhantomData,
			mapped: PhantomData,
		}
	}
}

impl<V, I: LinearIso<V>, D: IsDeg> Value<V, I, D> {
	pub fn add_change<T>(mut self, rate: impl Into<Value<V, I, T>>, unit: TimeUnit)
		-> Value::<V, I, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<V, I, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData,
			mapped: PhantomData,
		}
	}
	
	pub fn sub_change<T>(mut self, rate: impl Into<Value<V, I, T>>, unit: TimeUnit)
		-> Value::<V, I, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(-1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<V, I, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData,
			mapped: PhantomData,
		}
	}
}

impl<V, I: LinearIso<V>> FluxValue for Value<V, I, Deg<0>> {
	type Value = V;
	type Linear = I;
	type Degree = Deg<0>;
	fn value(&self) -> Self::Linear {
		self.value.initial_value
	}
	fn set_value(&mut self, value: Self::Linear) {
		self.value.initial_value = value;
	}
	fn change(&self, _changes: &mut Changes<Self>) {}
	fn update(&mut self, _time: Time) {}
}

macro_rules! impl_flux_value_for_value {
	($($degree:literal),+) => {$(
		impl<V, I: LinearIso<V>> FluxValue for Value<V, I, Deg<$degree>> {
			type Value = V;
			type Linear = I;
			type Degree = Deg<$degree>;
			fn value(&self) -> Self::Linear {
				self.value.initial_value
			}
			fn set_value(&mut self, value: Self::Linear) {
				self.value.initial_value = value;
			}
			fn change(&self, changes: &mut Changes<Self>) {
				for change in self.value.change_list.clone() {
					if change.scalar.0 == 1.0 {
						changes.add(
							&Value::<V, I, Deg<{ $degree - 1 }>> {
								value: change.rate,
								degree: PhantomData,
								mapped: PhantomData,
							},
							change.unit
						);
					} else {
						changes.sub(
							&Value::<V, I, Deg<{ $degree - 1 }>> {
								value: change.rate,
								degree: PhantomData,
								mapped: PhantomData,
							},
							change.unit
						);
					}
				}
			}
			fn update(&mut self, time: Time) {
				self.value.initial_value = self.calc_value(time);
				for change in &mut self.value.change_list {
					let mut v = Value::<V, I, Deg<{ $degree - 1 }>> {
						value: change.rate.clone(),
						degree: PhantomData,
						mapped: PhantomData,
					};
					v.update(time);
					change.rate = v.value;
				}
				// ??? Could probably reuse most of the initial_value calculation for
				// updating each sub-change's initial_value.
			}
		}
	)+};
}
impl_flux_value_for_value!(1, 2, 3, 4, 5, 6, 7);

/// Implements `From` for casting a [`Value`] to a higher degree.
macro_rules! impl_cast_up {
	($($degree:literal),+) => {$(
		impl<V, I, D> From<Value<V, I, Deg<{ $degree }>>> for Value<V, I, D>
		where
			I: LinearIso<V>,
			D: IsDeg,
			Deg<{ $degree }>: IsBelowDeg<D>,
		{
			fn from(value: Value<V, I, Deg<{ $degree }>>) -> Self {
				Self {
					value: value.value,
					degree: PhantomData,
					mapped: PhantomData,
				}
			}
		}
	)+}
}
impl_cast_up!(0, 1, 2, 3, 4, 5, 6, 7);

/// A value that linearly changes over time.
#[derive(Clone, Debug)]
struct LinearFluxValue<T: LinearValue> {
	initial_value: T,
	change_list: Vec<Change<T>>,
}

impl<T: LinearValue> LinearFluxValue<T> {
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
	rate: LinearFluxValue<T>,
	unit: TimeUnit,
}

impl<T: LinearValue> Change<T> {
	pub fn new(scalar: Scalar, rate: LinearFluxValue<T>, unit: TimeUnit) -> Self {
		Self {
			scalar,
			rate,
			unit,
		}
	}
}

/// Produces the degree of the resulting change, à la: `max(A, B+1)`. 
pub trait AddChange<D: IsDeg>: IsDeg {
	type Sum: IsDeg;
}

impl<A, B> AddChange<B> for A
where
	A: IsDeg,
	B: IsDeg + HasUpDeg,
	B::Up: MaxDeg<A>,
{
	type Sum = <B::Up as MaxDeg<A>>::Max;
}

/// A scalar value, used for multiplication with any [`LinearValue`].
#[derive(Copy, Clone, Debug)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) => Scalar(a.0 * b.0),
	(f64, Scalar)['commut] => a * b.0,
	(f32, Scalar)['commut] => a * (b.0 as f32),
	// (u64, Scalar)['commut] |
	// (u32, Scalar)['commut] |
	// (u16, Scalar)['commut] |
	// (u8,  Scalar)['commut] |
	// (i64, Scalar)['commut] |
	// (i32, Scalar)['commut] |
	// (i16, Scalar)['commut] |
	// (i8,  Scalar)['commut] => ((a as f64) * b.0).round() as Self,
}}

/// Any vector type that has addition and [`Scalar`] multiplication.
pub trait LinearValue:
	Sized + Copy + Clone + Debug
	+ Add<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	fn zero() -> Self;
}

impl<T> LinearValue for T
where T:
	Sized + Copy + Clone + Debug
	+ Add<Output=T>
	+ Mul<Scalar, Output=T>
	+ Default
{
	fn zero() -> Self {
		Self::default()
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// ??? Maybe rename to LinearMap, since it's not *really* isomorphic.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<Mapped>: LinearValue {
	fn map(self) -> Mapped;
	fn inv_map(value: Mapped) -> Self;
	fn set(&mut self, value: Mapped) {
		*self = Self::inv_map(value);
	}
}

impl<T: LinearValue> LinearIso<T> for T {
	fn map(self) -> T {
		self
	}
	fn inv_map(value: T) -> Self {
		value
	}
}

macro_rules! impl_iso_for_int {
	($a:ty: $($b:ty),+) => {$(
		impl LinearIso<$b> for $a {
			fn map(self) -> $b {
				self.round() as $b
			}
			fn inv_map(value: $b) -> Self {
				value as $a
			}
		}
	)+}
}
impl_iso_for_int!(f32: u8, u16, u32, u64, i8, i16, i32, i64);
impl_iso_for_int!(f64: u8, u16, u32, u64, i8, i16, i32, i64);

/// Used to construct a [`Flux`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(TimeUnit::Secs)` 
pub trait Per: Sized {
	fn per(self, unit: TimeUnit) -> Flux<Self> {
		Flux {
			rate: self,
			unit
		}
	}
}
impl Per for f64 {}
impl Per for u64 {}
impl Per for i64 {}
impl<V, I: LinearIso<V>, D: IsDeg> Per for Value<V, I, D> {}

/// A change over time used with arithmetic operators to construct [`Value`]s.
pub struct Flux<T> {
	rate: T,
	unit: TimeUnit,
}

macro_rules! impl_flux_ops {
	($op_trait:ident::$op_fn:ident, $change_fn:ident, $map:ty, $iso:ty) => {
		 // Num + Num.per(Unit)
		impl $op_trait<Flux<$map >> for $map
		where
			$iso: LinearIso<$map>
		{
			type Output = Value<$map, $iso, Deg<1>>;
			fn $op_fn(self, rhs: Flux<$map>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Num + Value.per(Unit)
		impl<D> $op_trait<Flux<Value<$map, $iso, D>>> for $map
		where
			D: IsDeg,
			Deg<0>: AddChange<D>,
			$iso: LinearIso<$map>,
		{
			type Output = Value<$map, $iso, <Deg<0> as AddChange<D>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$map, $iso, D>>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Value + Num.per(Unit)
		impl<D> $op_trait<Flux<$map >> for Value<$map, $iso, D>
		where
			D: IsDeg + AddChange<Deg<0>>,
			$iso: LinearIso<$map>,
		{
			type Output = Value<$map, $iso, <D as AddChange<Deg<0>>>::Sum>;
			fn $op_fn(self, rhs: Flux<$map>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Value + Value.per(Unit)
		impl<D, T> $op_trait<Flux<Value<$map, $iso, T>>> for Value<$map, $iso, D>
		where
			D: IsDeg + AddChange<T>,
			T: IsDeg,
			$iso: LinearIso<$map>,
		{
			type Output = Value<$map, $iso, <D as AddChange<T>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$map, $iso, T>>) -> Self::Output {
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

impl<V, I: LinearIso<V>> From<V> for Value<V, I, Deg<0>> {
	fn from(value: V) -> Self {
		Self::new(value)
	}
}

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: LinearValue>(T);

impl<T: LinearValue> Add for Exp<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}

impl<T: LinearValue> Mul<Scalar> for Exp<T> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self {
		Self(self.0 * rhs)
	}
}

impl<T: LinearValue> Default for Exp<T> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: LinearValue> From<T> for Exp<T> {
	fn from(value: T) -> Exp<T> {
		Exp(value)
	}
}

impl<T: LinearValue, D: IsDeg> Roots<D> for Exp<T>
where
	T: Roots<D>
{
	fn roots(poly: Poly<Exp<T>, D>) -> Result<RootList, RootList> {
		let mut b_poly = Poly::<T, D>::default();
		let mut b_coeff_iter = b_poly.coeff_iter_mut();
		for coeff in poly.coeff_iter() {
			*b_coeff_iter.next().unwrap() = coeff.0;
		}
		b_poly.0 = poly.constant().0;
		T::roots(b_poly)
	}
}

impl LinearIso<f64> for Exp<f64> {
	fn map(self) -> f64 {
		self.0.exp()
	}
	fn inv_map(value: f64) -> Self {
		Self(value.ln())
	}
}

impl LinearIso<u64> for Exp<f64> {
	fn map(self) -> u64 {
		self.0.exp().round() as u64
	}
	fn inv_map(value: u64) -> Self {
		Self((value as f64).ln())
	}
}

#[cfg(test)]
mod value_tests {
	use super::*;
	use TimeUnit::*;
	use std::cmp::Ordering;
	use crate::flux::PredValue;
	use crate::polynomial::*;
	
	fn linear() -> (
		Value<i64, f64, Deg<1>>,
		Value<i64, f64, Deg<2>>,
		Value<i64, f64, Deg<3>>,
		Value<i64, f64, Deg<4>>,
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
		
		assert_eq!(a.poly(), Poly(
			30.0, [
			0.018166666666666664,
			3.668167918055556e-6,
			1.1764138889120372e-15,
			2.314814814814815e-29
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
		assert_eq!(v0.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]); // 0.6
		assert_eq!(v1.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]); // -1.902, 2.102
		assert_eq!(v2.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [3*Nanosecs]); // 3.892
		assert_eq!(v3.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [5*Nanosecs]); // -3.687, 5.631
		v0.borrow_mut().update(1*Nanosecs);
		v1.borrow_mut().update(1*Nanosecs);
		v2.borrow_mut().update(1*Nanosecs);
		v3.borrow_mut().update(1*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [1*Nanosecs]);
		assert_eq!(v2.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [4*Nanosecs]);
		v0.borrow_mut().update(2*Nanosecs);
		v1.borrow_mut().update(2*Nanosecs);
		v2.borrow_mut().update(2*Nanosecs);
		v3.borrow_mut().update(2*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
	}
	
	fn exponential() -> (
		Value<f64, Exp<f64>, Deg<1>>,
		Value<f64, Exp<f64>, Deg<2>>,
		Value<f64, Exp<f64>, Deg<3>>,
		Value<f64, Exp<f64>, Deg<4>>,
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
		assert_eq!(v0.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [12*Nanosecs]); // 12.204
		assert_eq!(v1.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [24*Nanosecs]); // -1.143, 24.551
		assert_eq!(v2.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [36*Nanosecs]); // 36.91 
		assert_eq!(v3.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [49*Nanosecs]); // -4.752, 49.329
		v0.borrow_mut().update(6*Nanosecs);
		v1.borrow_mut().update(6*Nanosecs);
		v2.borrow_mut().update(6*Nanosecs);
		v3.borrow_mut().update(6*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [6*Nanosecs]);
		assert_eq!(v1.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [18*Nanosecs]);
		assert_eq!(v2.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [30*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [43*Nanosecs]);
		v0.borrow_mut().update(20*Nanosecs);
		v1.borrow_mut().update(20*Nanosecs);
		v2.borrow_mut().update(20*Nanosecs);
		v3.borrow_mut().update(20*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [10*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(1.0_f64)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]);
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
		
		v0.update(10*Nanosecs);
		let new_v0 = PredValue::new(v0.clone() + 7.per(Nanosecs));
		assert_eq!(v0.at(0*Nanosecs), new_v0.at(0*Nanosecs));
		assert_eq!(new_v0.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]); // 23.5
		
		v3.update(6*Nanosecs);
		let new_v3 = PredValue::new(v3.clone() + 1000.per(Nanosecs));
		assert_eq!(v3.at(0*Nanosecs), new_v3.at(0*Nanosecs));
		assert_eq!(new_v3.when_eq(&Value::new(0_i64)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs, 3*Nanosecs]); // 0.22, 3.684
	}
	
	#[test]
	fn when() {
		use Ordering::*;
		let (v0, v1, v2, v3) = linear();
		assert_eq!(PredValue::new(v0.clone()).when(Equal,   &Value::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v0.clone()).when(Less,    &Value::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v0.clone()).when(Greater, &Value::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(0*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Equal,   &Value::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, 2*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Less,    &Value::new(-5_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v1.clone()).when(Greater, &Value::new(-25_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 3*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Equal,   &Value::new(71_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs), (1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Less,    &Value::new(20_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(3*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v2.clone()).when(Greater, &Value::new(69_i64)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 2*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Equal,   &Value::new(140_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(5*Nanosecs, 5*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Less,    &Value::new(160_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 0*Nanosecs), (5*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(PredValue::new(v3.clone()).when(Greater, &Value::new(280_i64)).into_iter().collect::<Vec<(Time, Time)>>(), [(2*Nanosecs, 4*Nanosecs)]);
	}
	
	#[test]
	fn when_eq() {
		let (v0, v1, v2, v3) = linear();
		let v: [PredValue<Value<i64, f64, Deg<4>>>; 4] = [
			PredValue::new(v0.into()),
			PredValue::new(v1.into()),
			PredValue::new(v2.into()),
			PredValue::new(v3.into())
		];
		for value in v {
			assert_eq!(
				value.when_eq(&Value::new(3_i64)).into_iter().collect::<Vec<Time>>(),
				value.when(Ordering::Equal, &Value::new(3_i64)).into_iter()
					.map(|(a, _)| a)
					.collect::<Vec<Time>>()
			)
		}
	}
}