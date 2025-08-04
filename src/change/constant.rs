//! ...

/// ...
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Constant<T>(pub T);

mod _constant_impls {
	use super::Constant;
	use num_traits::{Inv, Pow};
	use std::ops::{Add, Div, Mul, Neg, Sub};
	use symb_poly::{CmpInvarInner, Exp, ExpInvar, FromInvar, InvarTerm, Ln, LnInvar, Num, NumKind};
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
	
	// ??? Temporary
	impl<T> FromInvar<Constant<T>> for Constant<T> {
		fn from_invar(value: Constant<T>) -> Self {
			value
		}
	}
	
	impl<T, N, D, R> FromInvar<Num<N, D, R>> for Constant<T>
	where
		T: Basis<Inner: From<Num<N, D, R>>>
	{
		fn from_invar(value: Num<N, D, R>) -> Self {
			Constant(T::from_inner(value.into()))
		}
	}

	impl<A, B, N> FromInvar<LnInvar<B, N>> for Constant<A>
	where
		Self: FromInvar<N> + FromInvar<B> + Ln<Output: Pow<Self, Output = Self>>
	{
		fn from_invar(value: LnInvar<B, N>) -> Self {
			Self::from_invar(value.0).ln().pow(Self::from_invar(value.1))
		}
	}

	impl<A, B> FromInvar<ExpInvar<B>> for Constant<A>
	where
		Self: FromInvar<B> + Exp<Output = Self>
	{
		fn from_invar(value: ExpInvar<B>) -> Self {
			Self::from_invar(value.0).exp()
		}
	}
	
	impl<T> InvarTerm for Constant<T>
	where
		T: Basis
	{
		type Kind = NumKind;
		fn is_indeterminant(&self) -> bool {
			self.0.is_indeterminant()
		}
	}
	
	// ??? Temporary
    impl<A, B, M, N> CmpInvarInner<Constant<B>, NumKind<M>, NumKind<N>> for Constant<A> {
        type Output = typenum::Less;
    }
}