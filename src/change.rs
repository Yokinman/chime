//! Change describing utilities ...

use std::marker::PhantomData;
use std::ops::{Add, Mul};

use impl_op::impl_op;

fn test() {
	let mut v = Value::new(3_f64);
	v.add_change(Change::new(ExponentialOp::Mul, Value::new(3_f64)));
	let mut v = Value::new(3_u8);
	v.add_change(Change::new(LinearOp::Add, Value::new(3_u8)));
	let c = Change::new(LinearOp::Add, Value::new(3_u64));
	let c = Change::new(ExponentialOp::Mul, Value::new(3_f64));
}

/// ...
struct Value<T: LinearIso, O: Operator<Term<T::Linear> = T>> {
	value: T,
	change: Vec<Change<T, O>>,
}

impl<T: LinearIso, O: Operator<Term<T::Linear> = T>> Value<T, O> {
	pub fn new(value: impl Into<T>) -> Self {
		Self {
			value: value.into(),
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
	time: u64,
}

impl<T: LinearIso, O: Operator<Term<T::Linear> = T>> Change<T, O> {
	pub fn new(operator: O, rate: Value<T, O>) -> Self {
		Self {
			operator,
			rate,
			time: 0,
		}
	}
}

/// ...
#[derive(Debug, Copy, Clone)]
pub struct Scalar(u64);

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
	Self: Add<Output=Self> + Mul<Scalar, Output=Self> + Sized,
{}

impl<T> LinearValue for T
where
	T: Add<Output=Self> + Mul<Scalar, Output=Self>,
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