//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use crate::{Constant, Flux};
use crate::linear::*;

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: Linear>(pub T);

impl<T: Linear> Linear for Exp<T> {
	fn add(self, other: Self) -> Self {
		Self(self.0.add(other.0))
	}
	fn sub(self, other: Self) -> Self {
		Self(self.0.sub(other.0))
	}
	fn mul_scalar(self, scalar: Scalar) -> Self {
		Self(self.0.mul_scalar(scalar))
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
	fn zero() -> Self {
		Self(<T as Linear>::zero())
	}
}

impl<T: Linear> Basis for Exp<T> {
	type Inner = Self;
	fn from_inner(inner: Self::Inner) -> Self {
		inner
	}
	fn into_inner(self) -> Self::Inner {
		self
	}
	fn with<R>(&self, f: impl FnOnce(&Self::Inner) -> R) -> R {
		f(&self)
	}
	fn inner_id(inner: Self::Inner) -> Self::Inner {
		inner
	}
}

impl<T: Linear> Default for Exp<T> {
	fn default() -> Self {
		Self(<T as Linear>::zero())
	}
}

impl<T: Linear> From<T> for Exp<T> {
	fn from(value: T) -> Exp<T> {
		Exp(value)
	}
}

impl<T> LinearIso<Exp<T>> for Exp<T>
where
	T: Linear,
{
	fn into_linear(value: Self) -> Exp<T> {
		value
	}
	fn from_linear(value: Exp<T>) -> Self {
		value
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

impl<T> Flux for Exp<T>
where
	T: Linear
{
	type Basis = Self;
	type Kind = Constant<Exp<T>>;
	fn basis(&self) -> Self::Basis {
		self.clone()
	}
	fn change(&self, basis: Self::Basis) -> Self::Kind {
		basis.into()
	}
}