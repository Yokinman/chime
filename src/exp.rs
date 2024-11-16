//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use crate::linear::*;

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T>(pub T);

impl<T: Linear> Linear for Exp<T> {
	fn add(self, other: Self) -> Self {
		Self(self.0.add(other.0))
	}
	fn sub(self, other: Self) -> Self {
		Self(self.0.sub(other.0))
	}
	fn mul(self, other: Self) -> Self {
		Self(self.0.mul(other.0))
	}
	fn div(self, other: Self) -> Self {
		Self(self.0.div(other.0))
	}
	fn pow(self, other: Self) -> Self {
		Self(self.0.pow(other.0))
	}
	fn mul_scalar(self, scalar: Scalar) -> Self {
		Self(self.0.mul_scalar(scalar))
	}
	fn pow_scalar(self, scalar: Scalar) -> Self {
		Self(self.0.pow_scalar(scalar))
	}
	fn exp(self) -> Self {
		Self(self.0.exp())
	}
	fn ln(self) -> Self {
		Self(self.0.ln())
	}
	fn sqr(mut self) -> Self {
		self.0 = self.0.sqr();
		self
	}
	fn sqrt(mut self) -> Self {
		self.0 = self.0.sqrt();
		self
	}
	fn sign(&self) -> Self {
		Exp(self.0.sign())
	}
	fn from_f64(n: f64) -> Self {
		Exp(T::from_f64(n))
	}
	fn zero() -> Self {
		Self(<T as Linear>::zero())
	}
}

impl<T> Simple for Exp<T> {}

impl<T: Linear> Default for Exp<T> {
	fn default() -> Self {
		Self(<T as Linear>::zero())
	}
}

impl<T> From<T> for Exp<T> {
	fn from(value: T) -> Exp<T> {
		Exp(value)
	}
}

impl LinearIso<Exp<f64>> for f64 {
	fn into_linear(value: f64) -> Exp<f64> {
		Exp(value.ln())
	}
	fn from_linear(value: Exp<f64>) -> f64 {
		value.0.exp()
	}
}

impl LinearIso<Exp<f64>> for u64 {
	fn into_linear(value: u64) -> Exp<f64> {
		Exp((value as f64).ln())
	}
	fn from_linear(value: Exp<f64>) -> u64 {
		value.0.exp().round() as u64
	}
}