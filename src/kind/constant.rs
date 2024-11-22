//! ...

use std::ops::{Add, Deref, DerefMut, Div, Mul, Neg, Sub};
use crate::{Flux, ToMoment, ToMomentMut};
use crate::exp::Exp;
use crate::kind::{FluxChange, FluxChangeUp, Poly};
use crate::linear::{Basis, Linear, Scalar, Vector};

/// ...
pub struct Nil<T>(std::marker::PhantomData<T>);

impl<T> Default for Nil<T> {
	fn default() -> Self {
		Self(std::marker::PhantomData)
	}
}

impl<T> FluxChange for Nil<T>
where
	T: Basis
{
	type Basis = T;
	type Poly = Constant<T>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		Constant(basis)
	}
	fn scale(self, _scalar: Scalar) -> Self {
		self
	}
}

impl<T> FluxChangeUp<'+'> for Nil<T>
where
	T: Basis
{
	type Up = crate::kind::sum::Sum<T, 1>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		crate::kind::sum::Sum([basis])
	}
}

impl<T> FluxChangeUp<'*'> for Nil<T>
where
	T: Basis
{
	type Up = Exp<crate::kind::sum::Sum<T, 1>>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Exp(crate::kind::sum::Sum([basis.map(Linear::ln)]))
	}
}

impl<T, U> Add<U> for Nil<T>
where
	T: Basis,
	U: FluxChange<Basis: Basis<Inner = T::Inner>>,
{
	type Output = U;
	fn add(self, rhs: U) -> Self::Output {
		rhs
	}
}

impl<T, U> Sub<U> for Nil<T>
where
	T: Basis,
	U: FluxChange<Basis: Basis<Inner = T::Inner>> + Neg<Output: FluxChange<Basis: Basis<Inner = T::Inner>>>,
{
	type Output = <U as Neg>::Output;
	fn sub(self, rhs: U) -> Self::Output {
		-rhs
	}
}

impl<T, U> Mul<Exp<U>> for Nil<T>
where
	T: Basis,
	U: FluxChange<Basis: Basis<Inner = T::Inner>>,
{
	type Output = Exp<U>;
	fn mul(self, rhs: Exp<U>) -> Self::Output {
		rhs
	}
}

impl<T, U> Div<Exp<U>> for Nil<T>
where
	T: Basis,
	U: FluxChange<Basis: Basis<Inner = T::Inner>> + Neg<Output: FluxChange<Basis: Basis<Inner = T::Inner>>>,
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

/// ...
pub struct ConstantIter<T>(T);

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

impl<T: Basis> Flux for Constant<T> {
	type Basis = T;
	type Change = Nil<T>;
	type Kind = Self;
	fn basis(&self) -> Self::Basis {
		self.0.clone()
	}
	fn change(&self) -> Self::Change {
		Nil::default()
	}
}

impl<T: Basis> ToMoment for Constant<T> {
	type Moment<'a> = Self;
	fn to_moment(&self, _time: Scalar) -> Self::Moment<'_> {
		self.clone()
	}
}

impl<T: Basis> ToMomentMut for Constant<T> {
	type MomentMut<'a> = &'a mut Self;
	fn to_moment_mut(&mut self, _time: Scalar) -> Self::MomentMut<'_> {
		self
	}
}

impl<T: Basis> Poly for Constant<T> {
	const DEGREE: usize = 0;
	type Basis = T;
	fn with_basis(value: Self::Basis) -> Self {
		Constant(value)
	}
	fn add_basis(mut self, basis: Self::Basis) -> Self {
		self.0 = self.0.zip_map(basis, T::Inner::add);
		self
	}
	fn deriv(self) -> Self {
		Self::zero()
	}
	fn eval(&self, _time: Scalar) -> Self::Basis {
		self.0.clone()
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

impl<T: Basis> Mul<Scalar> for Constant<T> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.0.map(|x| x.mul_scalar(rhs));
		self
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

impl<A, B> Add<B> for Constant<A>
where
	A: Basis,
	B: Poly<Basis = A>,
{
	type Output = B;
	fn add(self, rhs: B) -> Self::Output {
		rhs.add_basis(self.0)
	}
}

impl<A, B> Sub<B> for Constant<A>
where
	A: Basis,
	B: Poly<Basis = A> + Mul<Scalar, Output = B>,
{
	type Output = B;
	fn sub(self, rhs: B) -> Self::Output {
		(rhs * Scalar::from(-1.)).add_basis(self.0)
	}
}

impl<T> IntoIterator for Constant<T>
where
	T: IntoIterator,
{
	type Item = Constant<T::Item>;
	type IntoIter = ConstantIter<T::IntoIter>;
	fn into_iter(self) -> Self::IntoIter {
		ConstantIter(self.0.into_iter())
	}
}

impl<T> Iterator for ConstantIter<T>
where
	T: Iterator,
{
	type Item = Constant<T::Item>;
	fn next(&mut self) -> Option<Self::Item> {
		self.0.next().map(Constant)
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.0.size_hint()
	}
}