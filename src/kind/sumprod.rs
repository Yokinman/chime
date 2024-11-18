//! ... https://www.desmos.com/calculator/ko5owr56jx

use crate::{Flux, ToMoment, ToMomentMut};
use crate::kind::{FluxChange, Poly};
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
		SumProdPoly {
			basis,
			add_term: self.add_term,
			mul_term: self.mul_term,
		}
	}
	fn scale(self, scalar: Scalar) -> Self {
		Self {
			add_term: self.add_term.map(|x| x.mul_scalar(scalar)),
			mul_term: self.mul_term.map(|x| x.pow_scalar(scalar)),
		}
	}
}

/// ...
#[derive(Copy, Clone, Debug)]
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
		SumProd {
			add_term: self.add_term.clone(),
			mul_term: self.mul_term.clone(),
		}
	}
}

impl<T> Poly for SumProdPoly<T>
where
	T: Basis
{
	const DEGREE: usize = usize::MAX;
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
	fn deriv(self) -> Self {
		// let m = self.mul_term.clone().into_inner().ln() / ???;
		// self.basis = self.basis.map(|x| x.mul(m.clone()));
		// self.add_term = T::zero();
		// self
		todo!()
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		let m = self.mul_term.clone().into_inner();
		// if m == Linear::from_f64(1.) {
		// 	return a + bx
		// }
		let a = self.add_term.clone().into_inner()
			.mul(m.clone().div(m.clone().sub(Linear::from_f64(1.))));
		self.basis.clone().map(|x| x
			.add(a.clone())
			.mul(m.pow_scalar(time))
			.sub(a))
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

#[cfg(test)]
mod tests {
	use crate as chime;
	use crate::Flux;
	use crate::linear::Scalar;
	use crate::kind::Poly;
	
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
		let a = Test { value: 1., add: 4., mul: 0.5 }.to_kind();
		assert_eq!(a.eval(Scalar::from(0.)), 1.);
		assert_eq!(a.eval(Scalar::from(1.)), 2.5);
		assert_eq!(a.eval(Scalar::from(2.)), 3.25);
		assert_eq!(a.eval(Scalar::from(3.)), 3.625);
		assert_eq!(a.eval(Scalar::from(f64::INFINITY)), 4.);
		assert_eq!(a.eval(Scalar::from(f64::NEG_INFINITY)), f64::NEG_INFINITY);
		let a = Test { value: 1., add: 4., mul: 2. }.to_kind();
		assert_eq!(a.eval(Scalar::from(0.)), 1.);
		assert_eq!(a.eval(Scalar::from(1.)), 10.);
		assert_eq!(a.eval(Scalar::from(2.)), 28.);
		assert_eq!(a.eval(Scalar::from(3.)), 64.);
		assert_eq!(a.eval(Scalar::from(f64::INFINITY)), f64::INFINITY);
		assert_eq!(a.eval(Scalar::from(f64::NEG_INFINITY)), -8.);
	}
}