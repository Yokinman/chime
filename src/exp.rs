//! Exponential change over time.
//! 
//! Can be used to map a change over time from addition to multiplication. AKA
//! summation over time to multiplication over time.

use crate::change::{Change, ChangeUp};
use crate::linear::*;
use crate::poly::{Deriv, Poly, Translate};

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T>(pub T);

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

impl<T> Change for Exp<T>
where
	T: Change
{
	type Basis = T::Basis;
	type Poly = Exp<T::Poly>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		Exp(self.0.into_poly(basis.map_inner(Linear::ln)))
	}
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		Exp(self.0.scale(scalar))
	}
}

impl<T> ChangeUp<'*'> for Exp<T>
where
	T: ChangeUp<'+'>
{
	type Up = Exp<T::Up>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Exp(ChangeUp::<'+'>::up(self.0, basis.map_inner(Linear::ln)))
	}
}

impl<T> Poly for Exp<T>
where
	T: Poly
{
	const DEGREE: usize = T::DEGREE;
	type Basis = T::Basis;
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		self.0.eval(time).map_inner(Linear::exp)
	}
	fn zero() -> Self {
		Self(T::zero())
	}
}

impl<T> Deriv for Exp<T>
where
	T: Deriv
{
	type Deriv = Exp<T::Deriv>;
	fn deriv(self) -> Self::Deriv {
		Exp(Deriv::deriv(self.0))
	}
}

impl<T> Translate for Exp<T>
where
	T: Translate
{
	type Output = Exp<T::Output>;
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		Exp(self.0.translate(amount))
	}
}

/// ... natural logarithm
#[derive(Copy, Clone, Debug)]
pub struct Ln<T>(pub T);

impl<T> Poly for Ln<T>
where
	T: Poly
{
	const DEGREE: usize = T::DEGREE;
	
	type Basis = T::Basis;
	
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		self.0.eval(time).map_inner(Linear::ln)
	}
	
	fn zero() -> Self {
		Self(T::zero())
	}
}

impl<A, B, P> std::ops::Sub<B> for Ln<A>
where
	B: std::ops::Neg,
	Self: std::ops::Add<B::Output, Output=P>,
{
	type Output = P;
	fn sub(self, rhs: B) -> Self::Output {
		self + -rhs
	}
}

/// ...
pub trait ApplyExp {
	type Output;
	fn exp(self) -> Self::Output;
}

impl<T> ApplyExp for Ln<T> {
	type Output = T;
	fn exp(self) -> Self::Output {
		self.0
	}
}

impl<T> ApplyExp for T
where
	T: symb_poly::Exp
{
	type Output = T::Output;
	fn exp(self) -> Self::Output {
		T::exp(self)
	}
}

/// ...
pub trait ApplyLn {
	type Output;
	fn ln(self) -> Self::Output;
}

impl<T> ApplyLn for Exp<T> {
	type Output = T;
	fn ln(self) -> Self::Output {
		self.0
	}
}

impl<T> ApplyLn for T
where
	T: symb_poly::Ln
{
	type Output = T::Output;
	fn ln(self) -> Self::Output {
		T::ln(self)
	}
}