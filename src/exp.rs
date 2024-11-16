//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use crate::{Flux, ToMoment, ToMomentMut};
use crate::kind::{FluxChange, FluxChangeUp, FluxKind};
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

impl<T> FluxChange for Exp<T>
where
	T: FluxChange
{
	type Basis = T::Basis;
	type Poly = Exp<T::Poly>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		Exp(self.0.into_poly(basis))
	}
}

impl<T> FluxChangeUp for Exp<T>
where
	T: FluxChangeUp
{
	type Up = Exp<T::Up>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Exp(self.0.up(basis))
	}
}

impl<T> FluxKind for Exp<T>
where
	T: FluxKind
{
	const DEGREE: usize = T::DEGREE;
	fn with_basis(basis: Self::Basis) -> Self {
		Self(T::with_basis(Basis::from_inner(basis.into_inner().ln())))
	}
	fn add_basis(self, basis: Self::Basis) -> Self {
		Self(self.0.add_basis(Basis::from_inner(basis.into_inner().ln())))
	}
	fn deriv(self) -> Self {
		Self(self.0.deriv())
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		Basis::from_inner(self.0.eval(time).into_inner().exp())
	}
}

impl<T> Flux for Exp<T>
where
	T: Flux
{
	type Basis = T::Basis;
	type Change = Exp<T::Change>;
	type Kind = Exp<T::Kind>;
	fn basis(&self) -> Self::Basis {
		self.0.basis()
	}
	fn change(&self, basis: Self::Basis) -> Self::Kind {
		Exp(self.0.change(basis))
	}
}

impl<T> ToMoment for Exp<T>
where
	T: ToMoment
{
	type Moment<'a> = Exp<T::Moment<'a>> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		Exp(self.0.to_moment(time))
	}
}

impl<T> ToMomentMut for Exp<T>
where
	T: ToMomentMut
{
	type MomentMut<'a> = Exp<T::MomentMut<'a>> where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		Exp(self.0.to_moment_mut(time))
	}
}