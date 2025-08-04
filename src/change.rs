//! Defining a kind of change over time.

use crate::linear::Basis;
use crate::poly::Poly;

pub(crate) mod constant;
pub(crate) mod sum;
pub(crate) mod prod;
pub(crate) mod sumprod;

pub use sum::{Nil, Plus, Times};

/// A pattern of change.
pub trait Change<B: Basis> {
	type Poly: Poly<B>;
	fn into_poly(self, basis: B) -> Self::Poly;
	fn scale(self, scalar: B::Inner) -> Self;
}

mod _change_impls {
	use crate::change::sum::ChangeIntoPoly;
	use crate::linear::Basis;
	use crate::poly::Poly;
	use super::{Change, Nil, PopChange};
	use symb_poly::{Approach, Invar, Limit, Symbol, Var};
	
	impl<T> Change<T> for Nil
	where
		T: Basis + Copy,
		T::Inner: From<symb_poly::Num<typenum::Z0>> + From<symb_poly::Num<typenum::P1>> + From<symb_poly::Num<typenum::N1>>,
	{
		type Poly = Limit<Var<Symbol<typenum::U0>>, (Approach<Symbol<typenum::U0>, crate::change::constant::Constant<T>>,), Invar<crate::change::constant::Constant<T>>>;
		fn into_poly(self, basis: T) -> Self::Poly {
			Limit {
				expr: <Var<Symbol<typenum::U0>>>::default(),
				index: (Approach {
					symbol: <Symbol<typenum::U0>>::default(),
					value: crate::change::constant::Constant(basis.clone()),
				},),
				result: Invar(crate::change::constant::Constant(basis)),
			}
		}
		fn scale(self, _scalar: <T as Basis>::Inner) -> Self {
			self
		}
	}
	
	impl<A, B, T> Change<T> for (A, B)
	where
		A: Change<T> + PopChange,
		B: Change<T>,
		T: Basis,
		A::Inner: ChangeIntoPoly<A::Change, B, T, Output: Poly<T>>,
	{
		type Poly = <A::Inner as ChangeIntoPoly<A::Change, B, T>>::Output;
		fn into_poly(self, basis: T) -> Self::Poly {
			let (a, b) = self;
			let (a_inner, a_outer) = a.pop_change();
			a_inner.change_into_poly(a_outer, b, basis)
		}
		fn scale(self, scalar: <T as Basis>::Inner) -> Self {
			let (a, b) = self;
			(a.scale(scalar), b.scale(scalar))
		}
	}
	
	impl<T, const N: usize, B> Change<[B; N]> for [T; N]
	where
		T: Change<B>,
		B: Basis,
	{
		type Poly = [T::Poly; N];
		fn into_poly(self, basis: [B; N]) -> Self::Poly {
			let mut basis_iter = basis.into_iter();
			self.map(|x| x.into_poly(basis_iter.next().unwrap()))
		}
		fn scale(self, scalar: <[B; N] as Basis>::Inner) -> Self {
			self.map(|x| x.scale(scalar.clone()))
		}
	}
}

/// Combines two [`Change`] types into a sequence.
pub trait ApplyChange<C> {
	type Output;
	fn apply_change(self, change: C) -> Self::Output;
}

mod _impl_apply_change {
	use super::{ApplyChange, Nil, PopChange};
	
	impl<C> ApplyChange<C> for Nil {
		type Output = C;
		fn apply_change(self, change: C) -> Self::Output {
			change
		}
	}
	
	impl<A, B, C> ApplyChange<C> for (A, B)
	where
		B: ApplyChange<C, Output: PopChange>,
		A: ApplyChange<<B::Output as PopChange>::Inner>,
	{
		type Output = (A::Output, <B::Output as PopChange>::Change);
		fn apply_change(self, change: C) -> Self::Output {
			let (a, b) = self;
			let (b_inner, b_outer) = b.apply_change(change).pop_change();
			(a.apply_change(b_inner), b_outer)
		}
	}
}

/// Splits a [`Change`] into its basis change(s) and latter change.
pub trait PopChange {
	type Inner;
	type Change;
	fn pop_change(self) -> (Self::Inner, Self::Change);
}

mod _impl_pop_change {
	use super::{Nil, PopChange};
	
	impl PopChange for Nil {
		type Inner = Nil;
		type Change = Nil;
		fn pop_change(self) -> (Self::Inner, Self::Change) {
			(self, self)
		}
	}
	
	impl<A, B> PopChange for (A, B) {
		type Inner = A;
		type Change = B;
		fn pop_change(self) -> (Self::Inner, Self::Change) {
			self
		}
	}
}