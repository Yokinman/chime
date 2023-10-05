//! Change describing utilities ...

use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub, Div};

use impl_op::impl_op;

use time::{Time, TimeUnit};
use crate::flux::{ChangeAccum, FluxValue};

use crate::degree::{Deg, IsDeg, IsBelowDeg, HasUpDeg, MaxDeg};

/// A value with a dynamically-defined change over time.
#[derive(Clone, Debug)]
pub struct Value<I: LinearIso, D: IsDeg> {
	value: LinearFluxValue<I::Linear>,
	degree: PhantomData<D>,
}

impl<I: LinearIso> Value<I, Deg<0>> {
	pub fn new(initial_value: impl Into<I>) -> Self {
		Self {
			value: LinearFluxValue::new(initial_value.into().inv_map()),
			degree: PhantomData,
		}
	}
}

impl<I: LinearIso, D: IsDeg> Value<I, D> {
	pub fn add_change<T>(mut self, rate: impl Into<Value<I, T>>, unit: TimeUnit)
		-> Value::<I, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<I, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData
		}
	}
	
	pub fn sub_change<T>(mut self, rate: impl Into<Value<I, T>>, unit: TimeUnit)
		-> Value::<I, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(-1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<I, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData
		}
	}
}

// impl<Iso: LinearIso> AsRef<Value<Iso::Linear>> for FluxValue<Iso> {
// 	fn as_ref(&self) -> &Value<Iso::Linear> {
// 		&self.value
// 	}
// }

impl<I: LinearIso> FluxValue for Value<I, Deg<0>> {
	type Iso = I;
	type Degree = Deg<0>;
	fn value(&self) -> <Self::Iso as LinearIso>::Linear {
		self.value.initial_value
	}
	fn change(&self, _changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {}
	fn update(&mut self, _time: Time) {}
}

macro_rules! impl_flux_value_for_value {
	($($degree:literal),+) => {$(
		impl<I: LinearIso> FluxValue for Value<I, Deg<$degree>> {
			type Iso = I;
			type Degree = Deg<$degree>;
			fn value(&self) -> <Self::Iso as LinearIso>::Linear {
				self.value.initial_value
			}
			fn change(&self, changes: &mut ChangeAccum<Self::Iso, Self::Degree>) {
				for change in self.value.change_list.clone() {
					if change.scalar.0 == 1.0 {
						changes.add(&Value::<I, Deg<{ $degree - 1 }>> { value: change.rate, degree: PhantomData }, change.unit);
					} else {
						changes.sub(&Value::<I, Deg<{ $degree - 1 }>> { value: change.rate, degree: PhantomData }, change.unit);
					}
				}
			}
			fn update(&mut self, time: Time) {
				self.value.initial_value = self.calc(time);
				for change in &mut self.value.change_list {
					let mut v = Value::<I, Deg<{ $degree - 1 }>> { value: change.rate.clone(), degree: PhantomData };
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
		impl<I, D> From<Value<I, Deg<{ $degree }>>> for Value<I, D>
		where
			I: LinearIso,
			D: IsDeg,
			Deg<{ $degree }>: IsBelowDeg<D>,
		{
			fn from(value: Value<I, Deg<{ $degree }>>) -> Self {
				Self {
					value: value.value,
					degree: PhantomData
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
	(Scalar, Scalar) as (Scalar(a), Scalar(b)) => {
		if a == 0. || b == 0. {
			Scalar(0.)
		} else {
			Scalar(a * b)
		}
	},
	(f64, Scalar)['commut] => if b.0 == 0. { 0. } else { a * b.0 },
	(f32, Scalar)['commut] => if b.0 == 0. { 0. } else { a * (b.0 as f32) },
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
pub trait LinearValue
where Self:
	Sized
	+ Copy + Clone + Debug
	+ Add<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	const ZERO: Self; // ??? Could use Default, but it isn't a guaranteed zero
}

impl LinearValue for f64 {
	const ZERO: Self = 0.0;
}

/// A mapping of a vector space that preserves addition & multiplication.
/// ??? Maybe rename to LinearMap, since it's not *really* isomorphic.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso {
	type Linear: LinearValue;
	fn map(value: Self::Linear) -> Self;
	fn inv_map(self) -> Self::Linear;
}

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
impl<Iso: LinearIso, D: IsDeg> Per for Value<Iso, D> {}

/// A change over time used with arithmetic operators to construct [`Value`]s.
pub struct Flux<T> {
	rate: T,
	unit: TimeUnit,
}

macro_rules! impl_flux_ops {
	($op_trait:ident::$op_fn:ident, $change_fn:ident, $change_type:ident<$t:ty>) => {
		 // Num + Num
		impl $op_trait<Flux<$t>> for $t {
			type Output = Value<$change_type<$t>, Deg<1>>;
			fn $op_fn(self, rhs: Flux<$t>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Num + Value
		impl<D> $op_trait<Flux<Value<$change_type<$t>, D>>> for $t
		where
			D: IsDeg,
			Deg<0>: AddChange<D>,
		{
			type Output = Value<$change_type<$t>, <Deg<0> as AddChange<D>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$change_type<$t>, D>>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Value + Num
		impl<D> $op_trait<Flux<$t>> for Value<$change_type<$t>, D>
		where
			D: IsDeg + AddChange<Deg<0>>,
		{
			type Output = Value<$change_type<$t>, <D as AddChange<Deg<0>>>::Sum>;
			fn $op_fn(self, rhs: Flux<$t>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.unit)
			}
		}
		
		 // Value + Value
		impl<D, T> $op_trait<Flux<Value<$change_type<$t>, T>>> for Value<$change_type<$t>, D>
		where
			D: IsDeg + AddChange<T>,
			T: IsDeg,
		{
			type Output = Value<$change_type<$t>, <D as AddChange<T>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$change_type<$t>, T>>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.unit)
			}
		}
	}
}
impl_flux_ops!(Add::add, add_change, Linear<f64>);
impl_flux_ops!(Add::add, add_change, Linear<u64>);
impl_flux_ops!(Add::add, add_change, Linear<i64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<f64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<u64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<i64>);
impl_flux_ops!(Mul::mul, add_change, Exponential<f64>);
impl_flux_ops!(Mul::mul, add_change, Exponential<u64>);
impl_flux_ops!(Div::div, sub_change, Exponential<f64>);
impl_flux_ops!(Div::div, sub_change, Exponential<u64>);

/// A linear map that performs no change - pure addition.
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Linear<T>(T);

impl<T> From<T> for Linear<T> {
	fn from(value: T) -> Linear<T> {
		Linear(value)
	}
}

impl<T> From<T> for Value<Linear<T>, Deg<0>>
where
	Linear<T>: LinearIso
{
	fn from(value: T) -> Self {
		Self {
			value: LinearFluxValue::new(Linear(value).inv_map()),
			degree: PhantomData,
		}
	}
}

impl LinearIso for Linear<f64> {
	type Linear = f64;
	fn map(value: Self::Linear) -> Self { Self(value) }
	fn inv_map(self) -> Self::Linear { self.0 }
}

impl LinearIso for Linear<u64> {
	type Linear = f64;
	fn map(value: Self::Linear) -> Self { Self(value.round() as u64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

impl LinearIso for Linear<i64> {
	type Linear = f64;
	fn map(value: Self::Linear) -> Self { Self(value.round() as i64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exponential(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exponential<T>(T);

impl<T> From<T> for Exponential<T> {
	fn from(value: T) -> Exponential<T> {
		Exponential(value)
	}
}

impl<T> From<T> for Value<Exponential<T>, Deg<0>>
where
	Exponential<T>: LinearIso
{
	fn from(value: T) -> Self {
		Self {
			value: LinearFluxValue::new(Exponential(value).inv_map()),
			degree: PhantomData,
		}
	}
}

impl LinearIso for Exponential<f64> {
	type Linear = f64;
	fn map(value: Self::Linear) -> Self { Self(value.exp()) }
	fn inv_map(self) -> Self::Linear { self.0.ln() }
}

impl LinearIso for Exponential<u64> {
	type Linear = f64;
	fn map(value: Self::Linear) -> Self { Self(value.exp().round() as u64) }
	fn inv_map(self) -> Self::Linear { (self.0 as Self::Linear).ln() }
}

#[cfg(test)]
mod value_tests {
	use super::*;
	use TimeUnit::*;
	use std::cmp::Ordering;
	use crate::polynomial::*;
	
	fn linear() -> (
		Value<Linear<i64>, Deg<1>>,
		Value<Linear<i64>, Deg<2>>,
		Value<Linear<i64>, Deg<3>>,
		Value<Linear<i64>, Deg<4>>,
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
		let (mut v0, mut v1, mut v2, mut v3) = linear();
		assert_eq!(v0.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]); // 0.6
		assert_eq!(v1.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]); // -1.902, 2.102
		assert_eq!(v2.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [3*Nanosecs]); // 3.892
		assert_eq!(v3.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [5*Nanosecs]); // -3.687, 5.631
		v0.update(1*Nanosecs);
		v1.update(1*Nanosecs);
		v2.update(1*Nanosecs);
		v3.update(1*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [1*Nanosecs]);
		assert_eq!(v2.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [4*Nanosecs]);
		v0.update(2*Nanosecs);
		v1.update(2*Nanosecs);
		v2.update(2*Nanosecs);
		v3.update(2*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [2*Nanosecs]);
	}
	
	fn exponential() -> (
		Value<Exponential<f64>, Deg<1>>,
		Value<Exponential<f64>, Deg<2>>,
		Value<Exponential<f64>, Deg<3>>,
		Value<Exponential<f64>, Deg<4>>,
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
			assert!((v0.at(t*Nanosecs).0 - n0[t as usize]).abs() < 0.00001);
			assert!((v1.at(t*Nanosecs).0 - n1[t as usize]).abs() < 0.00001);
			assert!((v2.at(t*Nanosecs).0 - n2[t as usize]).abs() < 0.00001);
		}
		assert!(v3.at(100*Nanosecs).0.abs() < 0.00001);
	}
	
	#[test]
	fn exponential_roots() {
		let (mut v0, mut v1, mut v2, mut v3) = exponential();
		assert_eq!(v0.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [12*Nanosecs]); // 12.204
		assert_eq!(v1.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [24*Nanosecs]); // -1.143, 24.551
		assert_eq!(v2.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [36*Nanosecs]); // 36.91 
		assert_eq!(v3.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [49*Nanosecs]); // -4.752, 49.329
		v0.update(6*Nanosecs);
		v1.update(6*Nanosecs);
		v2.update(6*Nanosecs);
		v3.update(6*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [6*Nanosecs]);
		assert_eq!(v1.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [18*Nanosecs]);
		assert_eq!(v2.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [30*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [43*Nanosecs]);
		v0.update(20*Nanosecs);
		v1.update(20*Nanosecs);
		v2.update(20*Nanosecs);
		v3.update(20*Nanosecs);
		assert_eq!(v0.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v1.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), []);
		assert_eq!(v2.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [10*Nanosecs]);
		assert_eq!(v3.when_eq(&Value::new(1.0)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]);
	}
	
	#[test]
	fn exponential_non_pos() {
		let v = 0.0 * 2.0.per(Nanosecs);
		
		// assert_eq!(v.polynomial(), FluxPolynomial::Linear([f64::NEG_INFINITY, f64::ln(2.0)]));
		assert_eq!(v.at(0*Nanosecs).0, 0.0);
		assert_eq!(v.at(u64::MAX*Nanosecs).0, 0.0);
		
		let v = -1.5 * 2.0.per(Nanosecs);
		
		// match v.polynomial() {
		// 	FluxPolynomial::Linear([a, _]) => assert!(a.is_nan()),
		// 	_ => panic!()
		// }
		assert!(v.at(0*Nanosecs).0.is_nan());
		assert!(v.at(u64::MAX*Nanosecs).0.is_nan());
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
		let new_v0 = v0.clone() + 7.per(Nanosecs);
		assert_eq!(v0.at(0*Nanosecs), new_v0.at(0*Nanosecs));
		assert_eq!(new_v0.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [23*Nanosecs]); // 23.5
		
		v3.update(6*Nanosecs);
		let new_v3 = v3.clone() + 1000.per(Nanosecs);
		assert_eq!(v3.at(0*Nanosecs), new_v3.at(0*Nanosecs));
		assert_eq!(new_v3.when_eq(&Value::new(0)).into_iter().collect::<Vec<Time>>(), [0*Nanosecs, 3*Nanosecs]); // 0.22, 3.684
	}
	
	#[test]
	fn when() {
		use Ordering::*;
		let (v0, v1, v2, v3) = linear();
		assert_eq!(v0.when(Equal,   &Value::new(-5)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v0.when(Less,    &Value::new(-5)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v0.when(Greater, &Value::new(-5)).into_iter().collect::<Vec<(Time, Time)>>(),  [(0*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v1.when(Equal,   &Value::new(-5)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, 2*Nanosecs)]);
		assert_eq!(v1.when(Less,    &Value::new(-5)).into_iter().collect::<Vec<(Time, Time)>>(),  [(2*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v1.when(Greater, &Value::new(-25)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 3*Nanosecs)]);
		assert_eq!(v2.when(Equal,   &Value::new(71)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 1*Nanosecs), (1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v2.when(Less,    &Value::new(20)).into_iter().collect::<Vec<(Time, Time)>>(),  [(3*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v2.when(Greater, &Value::new(69)).into_iter().collect::<Vec<(Time, Time)>>(),  [(1*Nanosecs, 2*Nanosecs)]);
		assert_eq!(v3.when(Equal,   &Value::new(140)).into_iter().collect::<Vec<(Time, Time)>>(), [(5*Nanosecs, 5*Nanosecs)]);
		assert_eq!(v3.when(Less,    &Value::new(160)).into_iter().collect::<Vec<(Time, Time)>>(), [(0*Nanosecs, 0*Nanosecs), (5*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v3.when(Greater, &Value::new(280)).into_iter().collect::<Vec<(Time, Time)>>(), [(2*Nanosecs, 4*Nanosecs)]);
	}
	
	#[test]
	fn when_eq() {
		let (v0, v1, v2, v3) = linear();
		let v: [Value<Linear<i64>, Deg<4>>; 4] = [
			v0.into(),
			v1.into(),
			v2.into(),
			v3.into()
		];
		for value in v {
			assert_eq!(
				value.when_eq(&Value::new(3)).into_iter().collect::<Vec<Time>>(),
				value.when(Ordering::Equal, &Value::new(3)).into_iter()
					.map(|(a, _)| a)
					.collect::<Vec<Time>>()
			)
		}
	}
}