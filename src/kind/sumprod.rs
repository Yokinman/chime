//! ... https://www.desmos.com/calculator/ko5owr56jx

use std::ops::{Add, Sub};
use crate::{ToMoment, ToMomentMut};
use crate::kind::{Change, FluxChangeUp, Poly, Roots};
use crate::kind::constant::Constant;
use crate::kind::sum::SumPoly;
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

impl<T: Basis> Change for SumProd<T> {
	type Basis = T;
	type Poly = SumProdPoly<Constant<T>>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		let mul_term = self.mul_term;
		let add_term = T::each_map_inner(
			[basis.clone(), self.add_term, mul_term.clone()],
			|[a, b, c]| {
				if c == Linear::from_f64(1.) {
					return b // a + bx
				}
				a.add(b.mul(c.clone()).div(c.sub(Linear::from_f64(1.))))
			}
		);
		SumProdPoly {
			basis: Constant(basis),
			add_term,
			mul_term,
		}
	}
	fn scale(self, scalar: Scalar) -> Self {
		Self {
			add_term: self.add_term.map_inner(|x| x.mul_scalar(scalar)),
			mul_term: self.mul_term.map_inner(|x| x.pow_scalar(scalar)),
		}
	}
}

impl<T: Basis> FluxChangeUp<'+'> for SumProd<T> {
	type Up = SumProd2<T>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		SumProd2 {
			basis,
			change: self,
		}
	}
}

/// ...
pub struct SumProd2<T> {
	basis: T,
	change: SumProd<T>,
}

impl<T: Basis> Change for SumProd2<T> {
	type Basis = T;
	type Poly = SumProdPoly<SumPoly<T, 1>>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		let deriv = self.change.into_poly(self.basis);
		let sum_term = T::each_map_inner(
			[deriv.basis.0, deriv.add_term.clone(), deriv.mul_term.clone()],
			|[a, b, c]| {
				if c == Linear::from_f64(1.) {
					return a
				}
				a.sub(b)
			}
		);
		let add_term = deriv.add_term.zip_map_inner(deriv.mul_term.clone(), |b, c| {
			if c == Linear::from_f64(1.) {
				return b.mul_scalar(Scalar::from(0.5))
			}
			b.div(c.ln())
		});
		let mul_term = deriv.mul_term;
		// This is specifically the integral of `SumProd` to avoid confusion
		// when doing predictions on a flux value and its change over time.
		SumProdPoly {
			basis: SumPoly::new(basis, [sum_term]),
			add_term,
			mul_term,
		}
	}
	fn scale(mut self, scalar: Scalar) -> Self {
		self.basis = self.basis.map_inner(|x| x.mul_scalar(scalar));
		self.change = self.change.scale(scalar);
		self
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
		self.add_term = self.add_term.zip_map_inner(self.mul_term.clone(), |a, b| a.mul(b.ln()));
		self.basis = self.basis.deriv().add_basis(self.add_term.clone());
		self
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		Basis::each_map_inner(
			[self.basis.eval(time), self.add_term.clone(), self.mul_term.clone()],
			|[a, b, c]| {
				if c == Linear::from_f64(1.) {
					// Treat `*1` multipliers as `x^t` terms.
					return a.add(b.mul_scalar(match T::DEGREE {
						0 => time,
						1 => time*time,
						_ => unimplemented!("this is a temporary hacky way to do this")
						// !!! https://www.desmos.com/calculator/pgflaxsrxs
						// - Give `SumProdPoly` a `const DEGREE: usize` parameter.
						// - Use the `T` parameter's `FluxChangeUp::<'+'>::Up` type.
						// - Make the `SumProd` change impls return 2-variant enums.
					}))
				}
				a.add(b.mul(c.pow_scalar(time).sub(Linear::from_f64(1.))))
			}
		)
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
		if self.mul_term == 1. {
			return [-self.basis.0 / self.add_term];
		}
		[(1. - (self.basis.0 / self.add_term)).log(self.mul_term)] 
	}
}

impl Roots for SumProdPoly<SumPoly<f64, 1>> {
	type Output = [f64; 2];
	fn roots(self) -> <Self as Roots>::Output {
		if self.mul_term == 1. {
			return SumPoly::new(self.basis.0, [self.basis.1[0], self.add_term])
				.roots()
		}
		// Convert to `a + b*x + c^x` and solve using product log.
		let a = (self.basis.0 / self.add_term) - 1.;
		let b = self.basis.1[0] / self.add_term;
		let c = self.mul_term;
		let c_ln = c.ln();
		let z = c_ln / (b * c.powf(a / b));
		[
			-(lambert_w::lambert_wm1(z) / c_ln) - (a / b),
			-(lambert_w::lambert_w0(z)  / c_ln) - (a / b),
		]
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
	
	#[derive(Flux)]
	struct Test2 {
		#[basis]
		value: f64,
		#[change(add_per(chime::time::SEC))]
		add: Test,
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
		
		for root in (a - chime::kind::constant::Constant(5.)).roots() {
			let val = a.eval(Scalar::from(root));
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let b = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 2. }
		}.to_poly();
		assert_eq!(b.eval(Scalar::from(0.)), 1.);
		assert_eq!(b.eval(Scalar::from(1.)), 5.984255368000671);
		assert_eq!(b.eval(Scalar::from(2.)), 23.95276610400201);
		assert_eq!(b.eval(Scalar::from(3.)), 67.8897875760047);
		assert_eq!(b.eval(Scalar::from(-1.)), 2.5078723159996645);
		assert_eq!(b.eval(Scalar::from(-2.)), 7.261808473999498);
		
		for root in (b - chime::kind::constant::Constant(5.)).roots() {
			let val = b.eval(Scalar::from(root));
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let c = Test { value: 1., add: 4., mul: 1. }.to_poly();
		assert_eq!(c.eval(Scalar::from(0.)), 1.);
		assert_eq!(c.eval(Scalar::from(1.)), 5.);
		assert_eq!(c.eval(Scalar::from(2.)), 9.);
		assert_eq!(c.eval(Scalar::from(3.)), 13.);
		
		for root in (c - chime::kind::constant::Constant(5.)).roots() {
			let val = c.eval(Scalar::from(root));
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let d = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 1. }
		}.to_poly();
		assert_eq!(d.eval(Scalar::from(0.)), 1.);
		assert_eq!(d.eval(Scalar::from(1.)), 4.);
		assert_eq!(d.eval(Scalar::from(2.)), 11.);
		assert_eq!(d.eval(Scalar::from(3.)), 22.);
		
		for root in (d - chime::kind::constant::Constant(5.)).roots() {
			let val = d.eval(Scalar::from(root));
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
	}
}