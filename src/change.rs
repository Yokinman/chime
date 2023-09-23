//! Change describing utilities ...

use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Add, Mul};

use impl_op::impl_op;

use time::{Time, TimeUnit};

use crate::polynomial::{FluxPolynomial, Degree};

/// ...
#[derive(Debug, Clone)] // ??? Copy may be worth deriving
struct Value<Iso: LinearIso> {
	polynomial: FluxPolynomial<Iso::Linear>,
	change_list: Vec<Change<Iso>>,
}

impl<Iso: LinearIso> Value<Iso> {
	pub fn new(value: impl Into<Iso>) -> Self {
		Self {
			polynomial: FluxPolynomial::Constant([value.into().inv_map()]),
			change_list: vec![],
		}
	}
	
	pub fn degree(&self) -> Degree {
		self.polynomial.degree()
	}
	
	fn polynomial(&self) -> FluxPolynomial<Iso::Linear> {
		self.polynomial
	}
	
	fn at(&self, time: Time) -> Iso {
		Iso::map(self.polynomial.at(time /* - time_offset */))
	}
	
	pub fn add_change(mut self, change: Change<Iso>) -> Result<Self, Self> {
		let scalar = change.operator.scalar();
		for term_index in 0..=change.rate.polynomial.degree() {
			let term = change.rate.polynomial.term(term_index) * scalar;
			match self.polynomial.add_term(term_index + 1, term) {
				Ok(polynomial) => self.polynomial = polynomial,
				Err(..) => return Err(self),
			}
		}
		self.change_list.push(change);
		Ok(self)
	}
}

/// ...
#[derive(Debug, Clone)]
struct Change<Iso: LinearIso> {
	operator: Iso::Op,
	rate: Value<Iso>,
	time_unit: TimeUnit,
}

impl<Iso: LinearIso> Change<Iso> {
	pub fn new(operator: Iso::Op, rate: Value<Iso>, time_unit: TimeUnit) -> Self {
		Self {
			operator,
			rate,
			time_unit,
		}
	}
}

/// ...
#[derive(Debug, Copy, Clone)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) => Scalar(a.0 * b.0),
	(f64, Scalar)['commut] => a * b.0,
	(f32, Scalar)['commut] => ((a as f64) * b.0) as f32,
	// (u64, Scalar)['commut] |
	// (u32, Scalar)['commut] |
	// (u16, Scalar)['commut] |
	// (u8,  Scalar)['commut] |
	// (i64, Scalar)['commut] |
	// (i32, Scalar)['commut] |
	// (i16, Scalar)['commut] |
	// (i8,  Scalar)['commut] => ((a as f64) * b.0).round() as Self,
}}

/// Any vector type that has addition & scalar multiplication.
pub trait LinearValue
where
	Self:
		Sized + Copy + Clone
		+ Add<Output=Self> + Mul<Scalar, Output=Self>
		// + Sub<Output=Self> + Div<Scalar, Output=Self>
		+ PartialEq + PartialOrd // !!! Temporary, wouldn't work for something like vectors
		+ Debug // ??? Could be optional
	// Scalar:
	// 	Mul<Self, Output=Self>,
{}

impl<T> LinearValue for T
where
	T:
		Copy + Clone
		+ Add<Output=T> + Mul<Scalar, Output=T>
		// + Sub<Output=T> + Div<Scalar, Output=T>
		+ PartialEq + PartialOrd
		+ Debug
	// Scalar:
	// 	Mul<T, Output=T>,
{}

/// The binary operators for addition & subtraction.
#[derive(Debug, Copy, Clone)]
pub enum LinearOp {
	Add,
	Sub,
}

impl Operator for LinearOp {
	fn scalar(&self) -> Scalar {
		match self {
			Self::Add => Scalar(1.0),
			Self::Sub => Scalar(-1.0),
		}
	}
}

/// The binary operators for multiplication & division.
#[derive(Debug, Copy, Clone)]
pub enum ExponentialOp {
	Mul,
	Div,
}

impl Operator for ExponentialOp {
	fn scalar(&self) -> Scalar {
		match self {
			Self::Mul => Scalar(1.0),
			Self::Div => Scalar(-1.0),
		}
	}
}

/// A binary operator described by a [`LinearIso`].
pub trait Operator: Debug + Copy + Clone {
	fn scalar(&self) -> Scalar;
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
	type Op: Operator;
	fn map(value: Self::Linear) -> Self;
	fn inv_map(self) -> Self::Linear;
	
	// /// Combine terms using an operator.
	// fn sum<T: LinearValue>(op: Self::Op, augend: Self, addend: Self) -> Self {
	// 	let augend = Self::inv_map(augend);
	// 	let addend = Self::inv_map(addend) * op.scalar();
	// 	Self::map(augend + addend)
	// }
}

/// A linear map that performs no change - pure addition.
#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub struct Linear<T>(T);

impl<T> From<T> for Linear<T> {
	fn from(value: T) -> Linear<T> {
		Linear(value)
	}
}

impl LinearIso for Linear<f64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value) }
	fn inv_map(self) -> Self::Linear { self.0 }
}

impl LinearIso for Linear<u64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value.round() as u64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

impl LinearIso for Linear<i64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value.round() as i64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exponential(A * B)`
#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub struct Exponential<T>(T);

impl<T> From<T> for Exponential<T> {
	fn from(value: T) -> Exponential<T> {
		Exponential(value)
	}
}

impl LinearIso for Exponential<f64> {
	type Linear = f64;
	type Op = ExponentialOp;
	fn map(value: Self::Linear) -> Self { Self(value.exp()) }
	fn inv_map(self) -> Self::Linear { self.0.ln() }
}

impl LinearIso for Exponential<u64> {
	type Linear = f64;
	type Op = ExponentialOp;
	fn map(value: Self::Linear) -> Self { Self(value.exp().round() as u64) }
	fn inv_map(self) -> Self::Linear { (self.0 as Self::Linear).ln() }
}

#[cfg(test)]
mod value_tests {
	use super::*;
	use TimeUnit::*;
	
	fn linear() -> [Value<Linear<i64>>; 4] {
		let mut v = Value::new(3)
			.add_change(Change::new(LinearOp::Sub, Value::new(5), Nanosecs))
			.unwrap();
		
		let mut v1 = Value::new(20)
			.add_change(Change::new(LinearOp::Add, v.clone(), Nanosecs))
			.unwrap()
			.add_change(Change::new(LinearOp::Add, v.clone(), Nanosecs))
			.unwrap();
		
		let mut v2 = Value::new(52)
			.add_change(Change::new(LinearOp::Add, v1.clone(), Nanosecs))
			.unwrap();
		
		let mut v3 = Value::new(150)
			.add_change(Change::new(LinearOp::Add, v2.clone(), Nanosecs))
			.unwrap()
			.add_change(Change::new(LinearOp::Add, v.clone(), Nanosecs))
			.unwrap();
		
		[v, v1, v2, v3]
	}
	
	#[test]
	fn degree_limit() {
		assert!(Value::new(200)
			.add_change(Change::new(LinearOp::Add, linear()[3].clone(), Nanosecs))
			.is_err());
	}
	
	#[test]
	fn linear_polynomial() {
		let v = linear();
		assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([3.0, -5.0]));
		assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([20.0, 6.0, -10.0]));
		assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([52.0, 20.0, 6.0, -10.0]));
		assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([150.0, 55.0, 15.0, 6.0, -10.0]));
		
		// assert_eq!(2_u64 + 3|Secs, Value::with_change(2_u64, Change::new(LinearOp::Add, Value::new(3), Secs)));
		// assert_eq!(2_u64 - 3|Secs, Value::with_change(2_u64, Change::new(LinearOp::Sub, Value::new(3), Secs)));
	}
	
	#[test]
	fn linear_value() {
		let v = linear();
		let n0 = [3, -2, -7, -12, -17, -22, -27, -32, -37, -42];
		let n1 = [20, 16, 2, -22, -56, -100, -154, -218, -292, -376];
		let n2 = [52, 68, 70, 48, -8, -108, -262, -480, -772, -1148];
		let n3 = [150, 216, 279, 315, 290, 160, -129, -641, -1450, -2640];
		for t in 0..10_u64 {
			assert_eq!(v[0].at(t*Nanosecs), n0[t as usize].into());
			assert_eq!(v[1].at(t*Nanosecs), n1[t as usize].into());
			assert_eq!(v[2].at(t*Nanosecs), n2[t as usize].into());
			assert_eq!(v[3].at(t*Nanosecs), n3[t as usize].into());
		}
	}
	
	fn exponential() -> [Value<Exponential<f64>>; 4] {
		let mut v = Value::new(3.2)
			.add_change(Change::new(ExponentialOp::Div, Value::new(1.1), Nanosecs))
			.unwrap();
		
		let mut v1 = Value::new(14.515)
			.add_change(Change::new(ExponentialOp::Mul, v.clone(), Nanosecs))
			.unwrap()
			.add_change(Change::new(ExponentialOp::Mul, v.clone(), Nanosecs))
			.unwrap();
		
		let mut v2 = Value::new(30.4)
			.add_change(Change::new(ExponentialOp::Mul, v1.clone(), Nanosecs))
			.unwrap();
		
		let mut v3 = Value::new(300.0)
			.add_change(Change::new(ExponentialOp::Mul, v2.clone(), Nanosecs))
			.unwrap()
			.add_change(Change::new(ExponentialOp::Div, v.clone(), Nanosecs))
			.unwrap();
		
		[v, v1, v2, v3]
	}
	
	#[test]
	fn exponential_polynomial() {
		let v = exponential();
		fn ln(n: f64) -> f64 { n.ln() }
		assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([ln(3.2), -ln(1.1)]));
		assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
		assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([ln(30.4), ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
		assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([ln(300.0), ln(30.4) - ln(3.2), ln(14.515) + ln(1.1), 2.0*ln(3.2), -2.0*ln(1.1)]));
	}
	
	#[test]
	fn exponential_value() {
		let v = exponential();
		let n0 = [3.2, 2.9090909, 2.6446281, 2.4042074, 2.1856431];
		let n1 = [14.515, 122.83769, 859.13387, 4965.97682, 23722.647933];
		let n2 = [30.4, 3734.26565, 3208234.11489, 15932016253.16265, 377949612447409.6];
		for t in 0..5_u64 {
			assert!((v[0].at(t*Nanosecs).0 - n0[t as usize]).abs() < 0.00001);
			assert!((v[1].at(t*Nanosecs).0 - n1[t as usize]).abs() < 0.00001);
			assert!((v[2].at(t*Nanosecs).0 - n2[t as usize]).abs() < 0.00001);
		}
		assert!(v[3].at(100*Nanosecs).0.abs() < 0.00001);
	}
	
	#[test]
	fn exponential_non_pos() {
		let v = Value::<Exponential<f64>>::new(0.0)
			.add_change(Change::new(ExponentialOp::Mul, Value::new(2.0), Nanosecs))
			.unwrap();
		
		assert_eq!(v.polynomial(), FluxPolynomial::Linear([f64::NEG_INFINITY, f64::ln(2.0)]));
		assert_eq!(v.at(0*Nanosecs).0, 0.0);
		assert_eq!(v.at(u64::MAX*Nanosecs).0, 0.0);
		
		let v = Value::<Exponential<f64>>::new(-1.5)
			.add_change(Change::new(ExponentialOp::Mul, Value::new(2.0), Nanosecs))
			.unwrap();
		
		match v.polynomial() {
			FluxPolynomial::Linear([a, _]) => assert!(a.is_nan()),
			_ => panic!()
		}
		assert!(v.at(0*Nanosecs).0.is_nan());
		assert!(v.at(u64::MAX*Nanosecs).0.is_nan());
	}
}