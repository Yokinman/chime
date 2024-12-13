//! ...

use std::ops::{Add, Deref, DerefMut, Div, Mul, Neg, Sub};
use crate::exp::Exp;
use crate::change::{Change, ChangeUp, Sum};
use crate::linear::{Basis, Linear, Vector};
use crate::poly::{Deriv, Poly, PolyOffset};

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
	fn add_basis(mut self, basis: Self::Basis) -> Self {
		self.0 = self.0.zip_map_inner(basis, T::Inner::add);
		self
	}
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

impl<T> PolyOffset for Constant<T>
where
	T: Basis
{
	type Offset = Self;
	fn offset(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Offset {
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

impl<T> Mul<Constant<T>> for Constant<T>
where
	T: Mul<Output = T>
{
	type Output = Self;
	fn mul(mut self, rhs: Self) -> Self {
		self.0 = self.0 * rhs.0;
		self
	}
}

impl<A, B> Add<Constant<B>> for Constant<A>
where
	A: Add<B>,
{
	type Output = Constant<A::Output>;
	fn add(self, rhs: Constant<B>) -> Self::Output {
		Constant(self.0 + rhs.0)
	}
}

impl<A, B> Sub<Constant<B>> for Constant<A>
where
	A: Sub<B>,
{
	type Output = Constant<A::Output>;
	fn sub(self, rhs: Constant<B>) -> Self::Output {
		Constant(self.0 - rhs.0)
	}
}