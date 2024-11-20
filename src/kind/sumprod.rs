//! ... https://www.desmos.com/calculator/ko5owr56jx

use std::ops::{Add, Sub};
use crate::{ToMoment, ToMomentMut};
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
	type Poly = SumProdPoly<Constant<T>>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		let mul_term = self.mul_term;
		let add_term = self.add_term.map(|x| x
			.mul(mul_term.clone().into_inner())
			.div(mul_term.clone().into_inner().sub(Linear::from_f64(1.)))
			.add(basis.clone().into_inner())
		);
		let basis = Constant(basis);
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
pub struct SumProdPoly<T: Poly> {
	pub(crate) basis: T,
	pub(crate) add_term: T::Basis,
	pub(crate) mul_term: T::Basis,
}

impl<T: Poly> Poly for SumProdPoly<T> {
	const DEGREE: usize = usize::MAX;
	type Basis = T::Basis;
	fn with_basis(basis: Self::Basis) -> Self {
		Self {
			basis: T::with_basis(basis),
			add_term: Basis::zero(),
			mul_term: Basis::zero(),
		}
	}
	fn add_basis(mut self, basis: Self::Basis) -> Self {
		self.basis = self.basis.add_basis(basis);
		self
	}
	fn deriv(mut self) -> Self {
		self.add_term = self.add_term.map(|x| x.mul(self.mul_term.clone().into_inner().ln()));
		self.basis = self.basis.deriv().add_basis(self.add_term.clone());
		self
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		// !!! if add_term == inf && mul_term == 1 { return basis + add_term*x }
		self.basis.eval(time).map(|x| x
			.add(self.add_term.clone().into_inner()
				.mul(self.mul_term.clone().into_inner().pow_scalar(time)
					.sub(<T::Basis as Basis>::Inner::from_f64(1.)))))
	}
}

impl<T: Poly> ToMoment for SumProdPoly<T> {
	type Moment<'a> = Self where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		let mut x = self.clone();
		let _ = x.to_moment_mut(time);
		x
	}
}

impl<T: Poly> ToMomentMut for SumProdPoly<T> {
	type MomentMut<'a> = &'a mut Self where Self: 'a;
	fn to_moment_mut(&mut self, _time: Scalar) -> Self::MomentMut<'_> {
		// self.basis = self.eval(time);
		// self
		todo!()
	}
}

impl<T> Add<Constant<T::Basis>> for SumProdPoly<T>
where
	T: Poly + Add<Constant<T::Basis>, Output = T>
{
	type Output = Self;
	fn add(mut self, rhs: Constant<T::Basis>) -> Self::Output {
		self.basis = self.basis + rhs;
		self
	}
}

impl<T> Sub<Constant<T::Basis>> for SumProdPoly<T>
where
	T: Poly + Sub<Constant<T::Basis>, Output = T>
{
	type Output = Self;
	fn sub(mut self, rhs: Constant<T::Basis>) -> Self::Output {
		self.basis = self.basis - rhs;
		self
	}
}

impl Roots for SumProdPoly<Constant<f64>> {
	type Output = [f64; 1];
	fn roots(self) -> <Self as Roots>::Output {
		// !!! if add_term == inf && mul_term == 1 { return [-self.basis / self.add_term] }
		[(1. - (self.basis.0 / self.add_term)).log(self.mul_term)] 
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