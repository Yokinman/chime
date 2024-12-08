//! ... https://www.desmos.com/calculator/ko5owr56jx

use std::ops::{Add, Sub};
use crate::change::{Change, ChangeUp};
use crate::change::constant::Constant;
use crate::change::sum::SumPoly;
use crate::linear::{Basis, Linear};
use crate::poly::{Deriv, Poly, PolyOffset, Roots};

/// The pattern of alternating addition and multiplication, `a = (a + b) * c`.
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
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		Self {
			add_term: self.add_term.map_inner(|x| x.mul(scalar.clone())),
			mul_term: self.mul_term.map_inner(|x| x.pow(scalar.clone())),
		}
	}
}

impl<T: Basis> ChangeUp<'+'> for SumProd<T> {
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
				return b.div(Linear::from_f64(2.))
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
	fn scale(mut self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		self.basis = self.basis.map_inner(|x| x.mul(scalar.clone()));
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
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		Basis::each_map_inner(
			[self.basis.eval(time), self.add_term.clone(), self.mul_term.clone()],
			|[a, b, c]| {
				if c == Linear::from_f64(1.) {
					// Treat `*1` multipliers as `x^t` terms.
					return a.add(b.mul(match T::DEGREE {
						0 => time,
						1 => time.mul(time),
						_ => unimplemented!("this is a temporary hacky way to do this")
						// !!! https://www.desmos.com/calculator/pgflaxsrxs
						// - Give `SumProdPoly` a `const DEGREE: usize` parameter.
						// - Use the `T` parameter's `ChangeUp::<'+'>::Up` type.
						// - Make the `SumProd` change impls return 2-variant enums.
					}))
				}
				a.add(b.mul(c.pow(time).sub(Linear::from_f64(1.))))
			}
		)
	}
}

impl<T> Deriv for SumProdPoly<T>
where
	T: Deriv
{
	type Deriv = SumProdPoly<T::Deriv>;
	fn deriv(self) -> Self::Deriv {
		let add_term = self.add_term.zip_map_inner(self.mul_term.clone(), |a, b| a.mul(b.ln()));
		SumProdPoly {
			add_term: add_term.clone(),
			mul_term: self.mul_term,
			basis: Deriv::deriv(self.basis).add_basis(add_term),
		}
	}
}

impl<T> PolyOffset for SumProdPoly<T>
where
	T: PolyOffset,
	Self: Add<T::Offset, Output: Poly<Basis = Self::Basis>>,
{
	type Offset = <Self as Add<T::Offset>>::Output;
	fn offset(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Offset {
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
	use crate::poly::{Poly, Roots};

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
		assert_eq!(a.eval(0.), 1.);
		assert_eq!(a.eval(1.), 2.5);
		assert_eq!(a.eval(2.), 3.25);
		assert_eq!(a.eval(3.), 3.625);
		assert_eq!(a.eval(f64::INFINITY), 4.);
		assert_eq!(a.eval(f64::NEG_INFINITY), f64::NEG_INFINITY);
		let a = Test { value: 1., add: 4., mul: 2. }.to_poly();
		assert_eq!(a.eval(0.), 1.);
		assert_eq!(a.eval(1.), 10.);
		assert_eq!(a.eval(2.), 28.);
		assert_eq!(a.eval(3.), 64.);
		assert_eq!(a.eval(f64::INFINITY), f64::INFINITY);
		assert_eq!(a.eval(f64::NEG_INFINITY), -8.);
		
		for root in (a - chime::change::constant::Constant(5.)).roots() {
			let val = a.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let b = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 2. }
		}.to_poly();
		assert_eq!(b.eval(0.), 1.);
		assert_eq!(b.eval(1.), 5.984255368000671);
		assert_eq!(b.eval(2.), 23.95276610400201);
		assert_eq!(b.eval(3.), 67.8897875760047);
		assert_eq!(b.eval(-1.), 2.5078723159996645);
		assert_eq!(b.eval(-2.), 7.261808473999498);
		
		for root in (b - chime::change::constant::Constant(5.)).roots() {
			let val = b.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let c = Test { value: 1., add: 4., mul: 1. }.to_poly();
		assert_eq!(c.eval(0.), 1.);
		assert_eq!(c.eval(1.), 5.);
		assert_eq!(c.eval(2.), 9.);
		assert_eq!(c.eval(3.), 13.);
		
		for root in (c - chime::change::constant::Constant(5.)).roots() {
			let val = c.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
		
		let d = Test2 {
			value: 1.,
			add: Test { value: 1., add: 4., mul: 1. }
		}.to_poly();
		assert_eq!(d.eval(0.), 1.);
		assert_eq!(d.eval(1.), 4.);
		assert_eq!(d.eval(2.), 11.);
		assert_eq!(d.eval(3.), 22.);
		
		for root in (d - chime::change::constant::Constant(5.)).roots() {
			let val = d.eval(root);
			assert!((val - 5.).abs() < 1e-12, "{:?}", val);
		}
	}
}