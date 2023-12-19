//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use std::ops::{Add, Mul};
use crate::linear::*;

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: Linear>(pub T);

impl<T: Linear> Linear for Exp<T> {
	fn sqrt(mut self) -> Self {
		self.0 = self.0.sqrt();
		self
	}
	fn next_up(mut self) -> Self {
		self.0 = self.0.next_up();
		self
	}
	fn next_down(mut self) -> Self {
		self.0 = self.0.next_down();
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

impl<T: Linear> Mul<Scalar> for Exp<T> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self {
		Self(self.0 * rhs)
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

impl LinearIso<f64> for Exp<f64> {
	fn map(self) -> f64 {
		self.0.exp()
	}
}

impl LinearIsoInv<Exp<f64>> for f64 {
	fn inv_map(self) -> Exp<f64> {
		Exp(self.ln())
	}
}

impl LinearIso<u64> for Exp<f64> {
	fn map(self) -> u64 {
		self.0.exp().round() as u64
	}
}

impl LinearIsoInv<Exp<f64>> for u64 {
	fn inv_map(self) -> Exp<f64> {
		Exp((self as f64).ln())
	}
}