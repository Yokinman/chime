//! ...

use std::ops::{Add, Deref, DerefMut, Div, Mul, Neg, Sub};
use crate::exp::Exp;
use crate::change::{Change, ChangeUp, Sum};
use crate::change::sum::{AddPoly, Binomial, Monomial, MonomialAdd};
use crate::linear::{Basis, Linear, Vector};
use crate::poly::{Deriv, Poly, Translate};

/// ...
pub struct Nil<T>(std::marker::PhantomData<T>);

impl<T> Default for Nil<T> {
	fn default() -> Self {
		Self(std::marker::PhantomData)
	}
}

impl<T> Change for Nil<T>
where
	T: Basis
{
	type Basis = T;
	type Poly = Constant<T>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		Constant(basis)
	}
	fn scale(self, _scalar: <Self::Basis as Basis>::Inner) -> Self {
		self
	}
}

impl<T> ChangeUp<'+'> for Nil<T>
where
	T: Basis
{
	type Up = Sum<T, 1>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Sum([basis])
	}
}

impl<T> ChangeUp<'*'> for Nil<T>
where
	T: Basis
{
	type Up = Exp<Sum<T, 1>>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Exp(Sum([basis.map_inner(Linear::ln)]))
	}
}

impl<T, U> Add<U> for Nil<T>
where
	T: Basis,
	U: Change<Basis: Basis<Inner = T::Inner>>,
{
	type Output = U;
	fn add(self, rhs: U) -> Self::Output {
		rhs
	}
}

impl<T, U> Sub<U> for Nil<T>
where
	T: Basis,
	U: Change<Basis: Basis<Inner = T::Inner>>
		+ Neg<Output: Change<Basis: Basis<Inner = T::Inner>>>,
{
	type Output = <U as Neg>::Output;
	fn sub(self, rhs: U) -> Self::Output {
		-rhs
	}
}

impl<T, U> Mul<Exp<U>> for Nil<T>
where
	T: Basis,
	U: Change<Basis: Basis<Inner = T::Inner>>,
{
	type Output = Exp<U>;
	fn mul(self, rhs: Exp<U>) -> Self::Output {
		rhs
	}
}

impl<T, U> Div<Exp<U>> for Nil<T>
where
	T: Basis,
	U: Change<Basis: Basis<Inner = T::Inner>>
		+ Neg<Output: Change<Basis: Basis<Inner = T::Inner>>>,
{
	type Output = Exp<<U as Neg>::Output>;
	fn div(self, rhs: Exp<U>) -> Self::Output {
		Exp(-rhs.0)
	}
}

/// No change over time.
/// 
/// Equivalent "constant" flux kinds should implement both `Into<Constant<T>>`
/// and `From<Constant<T>>` (e.g. `SumPoly<T, 0>`).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Constant<T>(pub T);

impl<T> Constant<T> {
	pub fn map<U>(self, f: impl FnOnce(T) -> U) -> Constant<U> {
		Constant(f(self.0))
	}
}

impl<T> Deref for Constant<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}

impl<T> DerefMut for Constant<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<T: Basis> From<T> for Constant<T> {
	fn from(value: T) -> Self {
		Constant(value)
	}
}

impl<T: Basis> Poly for Constant<T> {
	const DEGREE: usize = 0;
	type Basis = T;
	fn eval(&self, _time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		self.0.clone()
	}
	fn zero() -> Self {
		Self(T::zero())
	}
}

impl<T> Deriv for Constant<T>
where
	T: Basis
{
	type Deriv = Self;
	fn deriv(self) -> Self::Deriv {
		Self(T::zero())
	}
}

impl<T> Translate for Constant<T>
where
	T: Basis
{
	type Output = Self;
	fn translate(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		self
	}
}

impl<T, const SIZE: usize> Vector<SIZE> for Constant<T>
where
	T: Vector<SIZE>
{
	type Output = Constant<T::Output>;
	fn index(&self, index: usize) -> Self::Output {
		Constant(self.0.index(index))
	}
}

impl<T: Basis> Neg for Constant<T> {
	type Output = Self;
	fn neg(self) -> Self::Output {
		Self(self.0.map_inner(|x| x.mul(Linear::from_f64(-1.))))
	}
}

impl<T> MonomialAdd for Constant<T>
where
	T: Basis
{
	type Output = Self;
	fn monom_add(self, rhs: Self) -> Self::Output {
		Constant(self.0.zip_map_inner(rhs.0, Linear::add))
	}
}

impl<T, A, P> Add<A> for Constant<T>
where
	Self: AddPoly<A, Output=P>
{
	type Output = P;
	fn add(self, rhs: A) -> Self::Output {
		self.add_poly(rhs)
	}
}

impl<T, A, P> Sub<A> for Constant<T>
where
	A: Neg,
	Self: Add<<A as Neg>::Output, Output=P>,
{
	type Output = P;
	fn sub(self, rhs: A) -> Self::Output {
		self + -rhs
	}
}

impl<T> Mul for Constant<T>
where
	T: Basis
{
	type Output = Self;
	fn mul(self, rhs: Self) -> Self::Output {
		self.map(|n| n.zip_map_inner(rhs.0, Linear::mul))
	}
}

impl<T> Div for Constant<T>
where
	T: Basis
{
	type Output = Self;
	fn div(self, rhs: Self) -> Self::Output {
		self.map(|n| n.zip_map_inner(rhs.0, Linear::div))
	}
}

impl<T, const D: usize> Mul<Monomial<T, D>> for Constant<T>
where
	T: Basis
{
	type Output = Monomial<T, D>;
	fn mul(self, rhs: Monomial<T, D>) -> Self::Output {
		Monomial(self.0.zip_map_inner(rhs.0, Linear::mul))
	}
}

impl<T, A, B, X, Y> Mul<Binomial<A, B>> for Constant<T>
where
	Self: Clone + Mul<A, Output=X> + Mul<B, Output=Y>
{
	type Output = Binomial<X, Y>;
	fn mul(self, rhs: Binomial<A, B>) -> Self::Output {
		Binomial {
			lhs: self.clone() * rhs.lhs,
			rhs: self * rhs.rhs,
		}
	}
}