//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use std::ops::{Add, Mul, Sub};
use crate::linear::*;

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: Linear>(pub T);

impl<T: Linear> Linear for Exp<T> {
	fn mul(self, scalar: Scalar) -> Self {
		self * scalar
	}
	fn sqrt(mut self) -> Self {
		self.0 = self.0.sqrt();
		self
	}
	fn sign(mut self) -> Self {
		self.0 = self.0.sign();
		self
	}
	fn zero() -> Self {
		Self(T::zero())
	}
}

impl<T: Linear> Add for Exp<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}

impl<T: Linear> Sub for Exp<T> {
	type Output = Self;
	fn sub(self, rhs: Self) -> Self {
		Self(self.0 - rhs.0)
	}
}

impl<T: Linear> Mul<Scalar> for Exp<T> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self {
		Self(self.0.mul(rhs))
	}
}

impl<T: Linear> Default for Exp<T> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: Linear> From<T> for Exp<T> {
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