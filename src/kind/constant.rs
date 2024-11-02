//! ...

use std::ops::{Add, Deref, DerefMut, Mul, Sub};
use crate::{Change, Flux, ToMoment, ToMomentMut};
use crate::kind::{FluxAccum, FluxIntegral, FluxKind};
use crate::linear::{Basis, Linear, Scalar, Vector};

/// No change over time.
/// 
/// Equivalent "constant" flux kinds should implement both `Into<Constant<T>>`
/// and `From<Constant<T>>` (e.g. `Sum<T, 0>`).
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
	type Kind = Self;
	fn basis(&self) -> Self::Basis {
		self.0.clone()
	}
	fn change(&self, basis: Constant<Self::Basis>) -> Self::Kind {
		basis
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

impl<T: Basis> FluxKind for Constant<T> {
	const DEGREE: usize = 0;
	fn with_basis(value: Self::Basis) -> Self {
		Constant(value)
	}
	fn add_basis(mut self, value: Self::Basis) -> Self {
		self.0 = T::from_inner(self.0.into_inner().add(value.into_inner()));
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
		self.0 = T::from_inner(self.0.into_inner().mul_scalar(rhs));
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
	B: FluxKind<Basis = A>,
{
	type Output = B;
	fn add(self, rhs: B) -> Self::Output {
		rhs.add_basis(self.0)
	}
}

impl<A, B> Sub<B> for Constant<A>
where
	A: Basis,
	B: FluxKind<Basis = A> + Mul<Scalar, Output = B>,
{
	type Output = B;
	fn sub(self, rhs: B) -> Self::Output {
		(rhs * Scalar::from(-1.)).add_basis(self.0)
	}
}

impl<T, B> Add<Change<B>> for Constant<T>
where
	T: Basis,
	B: Flux<Kind: FluxIntegral>,
	Self: Add<<B::Kind as FluxIntegral>::Integ>,
{
	type Output = <Self as Add<<B::Kind as FluxIntegral>::Integ>>::Output;
	fn add(self, rhs: Change<B>) -> Self::Output {
		(FluxAccum(self) + rhs).0
	}
}

impl<T, B> Sub<Change<B>> for Constant<T>
where
	T: Basis,
	B: Flux<Kind: FluxIntegral>,
	Self: Sub<<B::Kind as FluxIntegral>::Integ>,
{
	type Output = <Self as Sub<<B::Kind as FluxIntegral>::Integ>>::Output;
	fn sub(self, rhs: Change<B>) -> Self::Output {
		(FluxAccum(self) - rhs).0
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