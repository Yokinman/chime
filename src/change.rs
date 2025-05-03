//! Defining a kind of change over time.

use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Sub};

use crate::linear::{Linear, Basis};
use crate::poly::Poly;

pub(crate) mod constant;
pub(crate) mod sum;
pub(crate) mod prod;
pub(crate) mod sumprod;

pub use constant::Nil;
pub use sum::{Sum2, Prod};

/// ...
pub trait Change {
	type Basis: Basis;
	type Poly: Poly<Basis = Self::Basis>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly;
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self;
}

mod _change_impls {
	use crate::linear::Basis;
	use super::Change;
	
	impl<T, const N: usize> Change for [T; N]
	where
		T: Change
	{
		type Basis = [T::Basis; N];
		type Poly = [T::Poly; N];
		fn into_poly(self, basis: Self::Basis) -> Self::Poly {
			let mut basis_iter = basis.into_iter();
			self.map(|x| x.into_poly(basis_iter.next().unwrap()))
		}
		fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
			self.map(|x| x.scale(scalar.clone()))
		}
	}
}

/// ...
pub trait ChangeUp<const OP: char>: Change {
	type Up: Change<Basis = Self::Basis>;
	fn up(self, basis: Self::Basis) -> Self::Up;
}

mod _flux_change_up_impls {
	use super::ChangeUp;
	
	impl<T> ChangeUp<'-'> for T
	where
		T: ChangeUp<'+'>
	{
		type Up = <T as ChangeUp<'+'>>::Up;
		fn up(self, basis: Self::Basis) -> Self::Up {
			ChangeUp::<'+'>::up(self, basis)
		}
	}
	
	impl<T> ChangeUp<'/'> for T
	where
		T: ChangeUp<'*'>
	{
		type Up = <T as ChangeUp<'*'>>::Up;
		fn up(self, basis: Self::Basis) -> Self::Up {
			ChangeUp::<'*'>::up(self, basis)
		}
	}
	
	impl<T, const N: usize> ChangeUp<'+'> for [T; N]
	where
		T: ChangeUp<'+'>
	{
		type Up = [T::Up; N];
		fn up(self, basis: Self::Basis) -> Self::Up {
			let mut basis_iter = basis.into_iter();
			self.map(|x| x.up(basis_iter.next().unwrap()))
		}
	}
	
	impl<T, const N: usize> ChangeUp<'*'> for [T; N]
	where
		T: ChangeUp<'*'>
	{
		type Up = [T::Up; N];
		fn up(self, basis: Self::Basis) -> Self::Up {
			let mut basis_iter = basis.into_iter();
			self.map(|x| x.up(basis_iter.next().unwrap()))
		}
	}
}

/// ...
pub trait ApplyChange<const OP: char, T>: Change {
	type Output: Change;
	fn apply_change(self, rhs: T) -> Self::Output;
}

impl<A, B> ApplyChange<'+', B> for A
where
	A: Change + Add<B, Output: Change>,
	B: Change,
{
	type Output = <A as Add<B>>::Output;
	fn apply_change(self, rhs: B) -> Self::Output {
		self + rhs
	}
}

impl<A, B> ApplyChange<'-', B> for A
where
	A: Change + Sub<B, Output: Change>,
	B: Change,
{
	type Output = <A as Sub<B>>::Output;
	fn apply_change(self, rhs: B) -> Self::Output {
		self - rhs
	}
}

impl<A, B> ApplyChange<'*', B> for A
where
	A: Change + Mul<B, Output: Change>,
	B: Change,
{
	type Output = <A as Mul<B>>::Output;
	fn apply_change(self, rhs: B) -> Self::Output {
		self * rhs
	}
}

impl<A, B> ApplyChange<'/', B> for A
where
	A: Change + Div<B, Output: Change>,
	B: Change,
{
	type Output = <A as Div<B>>::Output;
	fn apply_change(self, rhs: B) -> Self::Output {
		self / rhs
	}
}

/// ...
#[derive(Default)]
pub struct ChangeAccum<T>(T);

impl<T> ChangeAccum<T> {
	pub fn into_change(self) -> T {
		self.0
	}
}

impl<A, B> Add<B> for ChangeAccum<A>
where
	A: ApplyChange<'+', B>
{
	type Output = ChangeAccum<A::Output>;
	fn add(self, rhs: B) -> Self::Output {
		ChangeAccum(self.0.apply_change(rhs))
	}
}

impl<A, B> Sub<B> for ChangeAccum<A>
where
	A: ApplyChange<'-', B>
{
	type Output = ChangeAccum<A::Output>;
	fn sub(self, rhs: B) -> Self::Output {
		ChangeAccum(self.0.apply_change(rhs))
	}
}

impl<A, B> Mul<B> for ChangeAccum<A>
where
	A: ApplyChange<'*', B>
{
	type Output = ChangeAccum<A::Output>;
	fn mul(self, rhs: B) -> Self::Output {
		ChangeAccum(self.0.apply_change(rhs))
	}
}

impl<A, B> Div<B> for ChangeAccum<A>
where
	A: ApplyChange<'/', B>
{
	type Output = ChangeAccum<A::Output>;
	fn div(self, rhs: B) -> Self::Output {
		ChangeAccum(self.0.apply_change(rhs))
	}
}