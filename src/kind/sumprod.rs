//! ... https://www.desmos.com/calculator/ko5owr56jx

use std::ops::{Add, Sub};
use crate::{Flux, ToMoment, ToMomentMut};
use crate::kind::{FluxChange, Poly, Roots};
use crate::kind::constant::Constant;
use crate::linear::{Basis, Linear, Scalar};

/// ...
/// Represents the pattern
/// `a = (a + b)*c`
/// as
/// `f(x) = (a + bc/(c-1))*c^x - bc/(c-1)`
/// or maybe even
/// `f(x) = (a + b/ln(c))*c^x - b/ln(c)`
pub struct SumProd<T> {
	pub(crate) add_term: T,
	pub(crate) mul_term: T,
}

impl<T: Basis> FluxChange for SumProd<T> {
	type Basis = T;
	type Poly = SumProdPoly<T>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		let mul_term = self.mul_term;
		let add_term = self.add_term.map(|x| x
			.mul(mul_term.clone().into_inner())
			.div(mul_term.clone().into_inner().sub(Linear::from_f64(1.)))
			.add(basis.clone().into_inner())
		);
		SumProdPoly { basis, add_term, mul_term }
	}
	fn scale(self, scalar: Scalar) -> Self {
		Self {
			add_term: self.add_term.map(|x| x.mul_scalar(scalar)),
			mul_term: self.mul_term.map(|x| x.pow_scalar(scalar)),
		}
	}
}

/// ...
/// `b + a*(m^x - 1)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct SumProdPoly<T> {
	pub(crate) basis: T,
	pub(crate) add_term: T,
	pub(crate) mul_term: T,
}

impl<T> Flux for SumProdPoly<T>
where
	T: Basis
{
	type Basis = T;
	type Change = SumProd<T>;
	type Kind = Self;
	fn basis(&self) -> Self::Basis {
		self.basis.clone()
	}
	fn change(&self) -> Self::Change {
		// let mul_term = self.clone().mul_term;
		// let add_term = self.clone().add_term.map(|x| x
		// 	.mul(mul_term.clone().into_inner().sub(Linear::from_f64(1.)))
		// 	.div(mul_term.clone().into_inner())
		// );
		// SumProd { add_term, mul_term }
		todo!()
	}
}

impl<T> Poly for SumProdPoly<T>
where
	T: Basis
{
	const DEGREE: usize = usize::MAX;
	type Basis = T;
	fn with_basis(basis: Self::Basis) -> Self {
		Self {
			basis,
			add_term: T::zero(),
			mul_term: T::zero(),
		}
	}
	fn add_basis(mut self, basis: Self::Basis) -> Self {
		self.basis = self.basis.map(|x| x.add(basis.into_inner()));
		self
	}
	fn deriv(mut self) -> Self {
		self.add_term = self.add_term.map(|x| x.mul(self.mul_term.clone().into_inner().ln()));
		self.basis = self.add_term.clone();
		self
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		// !!! if add_term == inf && mul_term == 1 { return basis + add_term*x }
		self.basis.clone().map(|x| x
			.add(self.add_term.clone().into_inner()
				.mul(self.mul_term.clone().into_inner().pow_scalar(time)
					.sub(T::Inner::from_f64(1.)))))
	}
}

impl<T> ToMoment for SumProdPoly<T>
where
	T: Basis
{
	type Moment<'a> = Self where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		let mut x = self.clone();
		let _ = x.to_moment_mut(time);
		x
	}
}

impl<T> ToMomentMut for SumProdPoly<T>
where
	T: Basis
{
	type MomentMut<'a> = &'a mut Self where Self: 'a;
	fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
		self.basis = self.eval(time);
		self
	}
}

impl<T> Add<Constant<T>> for SumProdPoly<T>
where
	T: Basis
{
	type Output = Self;
	fn add(mut self, rhs: Constant<T>) -> Self::Output {
		self.basis = self.basis.map(|x| x.add(rhs.basis().into_inner()));
		self
	}
}

impl<T> Sub<Constant<T>> for SumProdPoly<T>
where
	T: Basis
{
	type Output = Self;
	fn sub(mut self, rhs: Constant<T>) -> Self::Output {
		self.basis = self.basis.map(|x| x.sub(rhs.basis().into_inner()));
		self
	}
}

impl Roots for SumProdPoly<f64> {
	type Output = [f64; 1];
	fn roots(self) -> <Self as Roots>::Output {
		// !!! if add_term == inf && mul_term == 1 { return [-self.basis / self.add_term] }
		[(1. - (self.basis / self.add_term)).log(self.mul_term)] 
	}
}

#[cfg(test)]
mod tests {
	use crate as chime;
	use crate::Flux;
	use crate::linear::Scalar;
	use crate::kind::{Poly, Roots};

	#[derive(Flux)]
	struct Test {
		#[basis]
		value: f64,
		#[change(add_per(chime::time::SEC))]
		add: f64,
		#[change(mul_per(chime::time::SEC))]
		mul: f64,
	}
	
	#[test]
	fn sumprod() {
		let a = Test { value: 1., add: 4., mul: 0.5 }.to_poly();
		assert_eq!(a.eval(Scalar::from(0.)), 1.);
		assert_eq!(a.eval(Scalar::from(1.)), 2.5);
		assert_eq!(a.eval(Scalar::from(2.)), 3.25);
		assert_eq!(a.eval(Scalar::from(3.)), 3.625);
		assert_eq!(a.eval(Scalar::from(f64::INFINITY)), 4.);
		assert_eq!(a.eval(Scalar::from(f64::NEG_INFINITY)), f64::NEG_INFINITY);
		let a = Test { value: 1., add: 4., mul: 2. }.to_poly();
		assert_eq!(a.eval(Scalar::from(0.)), 1.);
		assert_eq!(a.eval(Scalar::from(1.)), 10.);
		assert_eq!(a.eval(Scalar::from(2.)), 28.);
		assert_eq!(a.eval(Scalar::from(3.)), 64.);
		assert_eq!(a.eval(Scalar::from(f64::INFINITY)), f64::INFINITY);
		assert_eq!(a.eval(Scalar::from(f64::NEG_INFINITY)), -8.);
		
		assert_eq!(
			(a - chime::kind::constant::Constant(5.)).roots(),
			[0.5305147166987798],
		);
	}
}