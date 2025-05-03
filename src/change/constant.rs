//! ...

use std::ops::{Add, Div, Mul, Sub};
use symb_poly::Invar;
use crate::change::{Change, Sum2, Prod};
use crate::linear::{Basis, Linear};

/// ...
#[derive(Copy, Clone, Debug)]
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
	type Poly = Invar<Constant2<T>>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		Invar(Constant2(basis))
	}
	fn scale(self, _scalar: <Self::Basis as Basis>::Inner) -> Self {
		self
	}
}

impl<T, C> std::ops::Shr<C> for Nil<T> {
	type Output = C;
	fn shr(self, rhs: C) -> Self::Output {
		rhs
	}
}

impl<T, U> Add<crate::Rate<U>> for Nil<T>
where
	T: Basis,
	U: crate::Flux<Basis = T>,
	Sum2<Nil<T>, U::Change>: Change<Basis = T>,
{
	type Output = Sum2<Nil<T>, U::Change>;
	fn add(self, rhs: crate::Rate<U>) -> Self::Output {
		Sum2 {
			lhs: self,
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip()))
	}
}

impl<T, U> Sub<crate::Rate<U>> for Nil<T>
where
	T: Basis,
	U: crate::Flux<Basis = T>,
	Sum2<Nil<T>, U::Change>: Change<Basis = T>,
{
	type Output = Sum2<Nil<T>, U::Change>;
	fn sub(self, rhs: crate::Rate<U>) -> Self::Output {
		Sum2 {
			lhs: self,
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(-rhs.unit.as_secs_f64().recip()))
	}
}

impl<T, U> Mul<crate::Rate<U>> for Nil<T>
where
	T: Basis,
	U: crate::Flux<Basis = T>,
	Prod<Nil<T>, U::Change>: Change<Basis = T>,
{
	type Output = Prod<Nil<T>, U::Change>;
	fn mul(self, rhs: crate::Rate<U>) -> Self::Output {
		Prod {
			lhs: self,
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip()))
	}
}

impl<T, U> Div<crate::Rate<U>> for Nil<T>
where
	T: Basis,
	U: crate::Flux<Basis = T>,
	Prod<Nil<T>, U::Change>: Change<Basis = T>,
{
	type Output = Prod<Nil<T>, U::Change>;
	fn div(self, rhs: crate::Rate<U>) -> Self::Output {
		Prod {
			lhs: self,
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64()))
	}
}

/// ...
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Constant2<T>(pub T);

mod _constant_impls {
	use super::Constant2 as Constant;
	use num_traits::{Inv, Pow};
	use std::ops::{Add, Div, Mul, Neg, Sub};
	use symb_poly::{Exp, ExpInvar, FromInvar, InvarTerm, Ln, LnInvar, Num, NumKind};
	use crate::linear::{Basis, Linear};

	impl<T: Basis> Add for Constant<T> {
		type Output = Self;
		fn add(self, rhs: Constant<T>) -> Self::Output {
			Constant(self.0.zip_map_inner(rhs.0, Linear::add))
		}
	}

	impl<T: Basis> Sub for Constant<T> {
		type Output = Self;
		fn sub(self, rhs: Constant<T>) -> Self::Output {
			Constant(self.0.zip_map_inner(rhs.0, Linear::sub))
		}
	}

	impl<T: Basis> Mul for Constant<T> {
		type Output = Self;
		fn mul(self, rhs: Constant<T>) -> Self::Output {
			Constant(self.0.zip_map_inner(rhs.0, Linear::mul))
		}
	}

	impl<T: Basis> Div for Constant<T> {
		type Output = Self;
		fn div(self, rhs: Constant<T>) -> Self::Output {
			Constant(self.0.zip_map_inner(rhs.0, Linear::div))
		}
	}

	impl<T: Basis> Pow<Self> for Constant<T> {
		type Output = Self;
		fn pow(self, rhs: Constant<T>) -> Self::Output {
			Constant(self.0.zip_map_inner(rhs.0, Linear::pow))
		}
	}

	impl<T: Basis> Neg for Constant<T> {
		type Output = Self;
		fn neg(self) -> Self::Output {
			Constant(self.0.map_inner(|n| n.mul(Linear::from_f64(-1.))))
		}
	}

	impl<T: Basis> Inv for Constant<T> {
		type Output = Self;
		fn inv(self) -> Self::Output {
			Constant(self.0.map_inner(|n| n.pow(Linear::from_f64(-1.))))
		}
	}

	impl<T: Basis> Ln for Constant<T> {
		type Output = Self;
		fn ln(self) -> Self::Output {
			Constant(self.0.map_inner(Linear::ln))
		}
	}

	impl<T: Basis> Exp for Constant<T> {
		type Output = Self;
		fn exp(self) -> Self::Output {
			Constant(self.0.map_inner(Linear::exp))
		}
	}

	impl<T, N, D, R> FromInvar<Num<N, D, R>> for Constant<T>
	where
		T: Basis<Inner: From<Num<N, D, R>>>
	{
		type Output = Self;
		fn from_invar(value: Num<N, D, R>) -> Self::Output {
			Constant(T::from_inner(value.into()))
		}
	}

	impl<A, B, N> FromInvar<LnInvar<B, N>> for Constant<A>
	where
		Self: FromInvar<N>
			+ FromInvar<B, Output: Ln<Output: Pow<<Self as FromInvar<N>>::Output>>>
	{
		type Output = <<<Self as FromInvar<B>>::Output as Ln>::Output as Pow<<Self as FromInvar<N>>::Output>>::Output;
		fn from_invar(value: LnInvar<B, N>) -> Self::Output {
			Self::from_invar(value.0).ln().pow(Self::from_invar(value.1))
		}
	}

	impl<A, B> FromInvar<ExpInvar<B>> for Constant<A>
	where
		Self: FromInvar<B, Output: Exp>
	{
		type Output = <<Self as FromInvar<B>>::Output as Exp>::Output;
		fn from_invar(value: ExpInvar<B>) -> Self::Output {
			Self::from_invar(value.0).exp()
		}
	}

	impl<T> InvarTerm for Constant<T> {
		type Kind = NumKind;
	}
}