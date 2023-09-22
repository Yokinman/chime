//! Change describing utilities ...

use std::marker::PhantomData;
use std::ops::{Add, Sub, Mul, Div};

use impl_op::impl_op;

use time::TimeUnit;

use crate::polynomial::FluxPolynomial;

#[test]
fn test() {
	use TimeUnit::*;
	let mut v = Value::new(3_f64);
	v.add_change(Change::new(ExponentialOp::Mul, Value::new(3_f64), Secs));
	let mut v = Value::new(3_u8);
	v.add_change(Change::new(LinearOp::Add, Value::new(3_u8), Secs));
	let c = Change::new(LinearOp::Add, Value::new(3_u64), Secs);
	let c = Change::new(ExponentialOp::Mul, Value::new(3_f64), Secs);
	// assert_eq!(2_u64 + 3|Secs, Value::with_change(2_u64, Change::new(LinearOp::Add, Value::new(3), Secs)));
	// assert_eq!(2_u64 - 3|Secs, Value::with_change(2_u64, Change::new(LinearOp::Sub, Value::new(3), Secs)));
}

/// ...
struct Value<T: LinearIso, O: Operator<Term<T::Linear> = T>> {
	polynomial: FluxPolynomial<T::Linear>,
	change: Vec<Change<T, O>>,
}

impl<T: LinearIso, O: Operator<Term<T::Linear> = T>> Value<T, O> {
	pub fn new(value: impl Into<T>) -> Self {
		Self {
			polynomial: FluxPolynomial::Constant([value.into().inv_map()]),
			change: vec![],
		}
	}
	
	pub fn add_change(&mut self, change: Change<T, O>) {
		self.change.push(change)
	}
}

/// ...
struct Change<T: LinearIso, O: Operator<Term<T::Linear> = T>> {
	operator: O,
	rate: Value<T, O>,
	time_unit: TimeUnit,
}

impl<T: LinearIso, O: Operator<Term<T::Linear> = T>> Change<T, O> {
	pub fn new(operator: O, rate: Value<T, O>, time_unit: TimeUnit) -> Self {
		Self {
			operator,
			rate,
			time_unit,
		}
	}
}

/// ...
#[derive(Debug, Copy, Clone)]
pub struct Scalar(pub u64);

impl_op!{ a * b {
	(u64, Scalar)['commut] as (_, Scalar(b)) => a * b,
	(u32, Scalar)['commut] as (_, Scalar(b)) => a * (b as u32),
	(u16, Scalar)['commut] as (_, Scalar(b)) => a * (b as u16),
	(u8,  Scalar)['commut] as (_, Scalar(b)) => a * (b as u8),
	(f64, Scalar)['commut] as (_, Scalar(b)) => a * (b as f64),
	(f32, Scalar)['commut] as (_, Scalar(b)) => a * (b as f32),
}}

/*
pub enum Scalar {
	Pos(u64),
	Neg(u64),
}

impl_op!{ a * b {
	(u64, Scalar)['commut] => match b {
		Scalar::Pos(b) =>  a * b,
		Scalar::Neg(b) => -a * b,
	},
	(f64, Scalar)['commut] => match b {
		Scalar::Pos(b) =>  a * (b as f64),
		Scalar::Neg(b) => -a * (b as f64),
	},
}}
*/

/// Any vector type that has addition & scalar multiplication.
pub trait LinearValue
where
	Self:
		Sized + Copy + Clone
		+ Add<Output=Self> + Mul<Scalar, Output=Self>
		// + Sub<Output=Self> + Div<Scalar, Output=Self>
	// Scalar:
	// 	Mul<Self, Output=Self>,
{}

impl<T> LinearValue for T
where
	T:
		Copy + Clone
		+ Add<Output=T> + Mul<Scalar, Output=T>
		// + Sub<Output=T> + Div<Scalar, Output=T>
	// Scalar:
	// 	Mul<T, Output=T>,
{}

/// The binary operators for addition & subtraction.
enum LinearOp {
	Add,
	Sub,
}

impl Operator for LinearOp {
	type Term<T: LinearValue> = Linear<T>;
}

/// The binary operators for multiplication & division.
enum ExponentialOp {
	Mul,
	Div,
}

impl Operator for ExponentialOp {
	type Term<T: LinearValue> = Exponential;
}

/// A binary operator described by a [`LinearIso`].
trait Operator {
	type Term<T: LinearValue>: LinearIso;
	
	// fn scalar(&self) -> Scalar;
	
	fn sum<T: LinearValue>(&self, augend: Self::Term<T>, addend: Self::Term<T>) -> Self::Term<T> {
		//! Returns the result of the given terms combined by this operator.
		let augend = Self::Term::<T>::inv_map(augend);
		let addend = Self::Term::<T>::inv_map(addend) /* * self.scalar()*/;
		Self::Term::<T>::map(augend + addend)
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic  - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Vector addition       - `map(A + B) = map(A) • map(B)`
/// - Scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso {
	type Linear: LinearValue;
	fn map(value: Self::Linear) -> Self;
	fn inv_map(self) -> Self::Linear;
}

/// A linear map that performs no change - pure addition.
pub struct Linear<T: LinearValue>(T);

impl<T: LinearValue> LinearIso for Linear<T> {
	type Linear = T;
	fn map(value: T) -> Self { Self(value) }
	fn inv_map(self) -> T { self.0 }
}

impl<T: LinearValue> From<T> for Linear<T> {
	fn from(value: T) -> Linear<T> {
		Linear(value)
	}
}

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exponential(A * B)`
pub struct Exponential(f64);

impl LinearIso for Exponential {
	type Linear = f64;
	fn map(value: f64) -> Self { Self(value.exp()) }
	fn inv_map(self) -> f64 { self.0.ln() }
}

impl_op!{ a -> Exponential { f64 => Exponential(a) } }