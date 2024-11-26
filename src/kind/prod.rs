//! ...

use crate::exp::Exp;
use crate::kind::FluxChange;
use crate::kind::sum::{Sum, SumPoly};
use crate::kind::sumprod::SumProd;
use crate::linear::{Basis, Linear};

/// Represents the pattern of repeated multiplication, `a = a * b`.
pub type Prod<T, const DEGREE: usize> = Exp<Sum<T, DEGREE>>;

/// ...
pub type ProdPoly<T, const DEGREE: usize> = Exp<SumPoly<T, DEGREE>>;

impl<T: FluxChange + std::ops::Neg> std::ops::Neg for Exp<T> {
	type Output = Exp<<T as std::ops::Neg>::Output>;
	fn neg(self) -> Self::Output {
		Exp(-self.0)
	}
}

impl<T: Basis> std::ops::Mul<Prod<T, 1>> for Sum<T, 1> {
	type Output = SumProd<T>;
	fn mul(self, rhs: Prod<T, 1>) -> Self::Output {
		let Sum([add_term]) = self;
		let Exp(Sum([mul_term])) = rhs;
		SumProd {
			add_term,
			mul_term: mul_term.map(Linear::exp),
		}
	}
}

impl<T: Basis> std::ops::Div<Prod<T, 1>> for Sum<T, 1> {
	type Output = SumProd<T>;
	fn div(self, rhs: Prod<T, 1>) -> Self::Output {
		let Sum([add_term]) = self;
		let Exp(Sum([mul_term])) = rhs;
		SumProd {
			add_term,
			mul_term: mul_term.map(|x| T::Inner::from_f64(1.).div(x.exp())),
		}
	}
}

impl<A, B> std::ops::Mul<Exp<B>> for Exp<A>
where
	A: FluxChange + std::ops::Add<B>,
	B: FluxChange,
{
	type Output = Exp<A::Output>;
	fn mul(self, rhs: Exp<B>) -> Self::Output {
		Exp(self.0 + rhs.0)
	}
}

impl<A, B> std::ops::Div<Exp<B>> for Exp<A>
where
	A: FluxChange + std::ops::Sub<B>,
	B: FluxChange,
{
	type Output = Exp<A::Output>;
	fn div(self, rhs: Exp<B>) -> Self::Output {
		Exp(self.0 - rhs.0)
	}
}

#[cfg(test)]
mod _test {
	use crate as chime;
	use crate::Flux;
	use crate::kind::{FluxChange, Poly};
	use crate::linear::Scalar;
	use super::SumPoly;
	
	#[derive(Flux)]
	pub struct Test {
		#[basis]
		value: f64,
		#[change(mul_per(crate::time::SEC))]
		mul: f64,
	}
	
	#[test]
	fn prod() {
		let a = SumPoly::new(5., [1.5]);
		assert_eq!(a.eval(Scalar::from(0.)), 5.);
		assert_eq!(a.eval(Scalar::from(1.)), 6.5);
		assert_eq!(a.eval(Scalar::from(2.)), 8.);
		let b = Test { value: 1., mul: 2. }.to_poly();
		assert_eq!(b.eval(Scalar::from(0.)), 1.);
		assert_eq!(b.eval(Scalar::from(1.)), 2.);
		assert_eq!(b.eval(Scalar::from(2.)), 4.);
		let c = (a.change() * b.change()).into_poly(a.basis());
		assert_eq!(c.eval(Scalar::from(0.)), 5.);
		assert_eq!(c.eval(Scalar::from(1.)), 13.);
		assert_eq!(c.eval(Scalar::from(2.)), 29.);
		let d = (a.change() / b.change()).into_poly(a.basis());
		assert_eq!(d.eval(Scalar::from(0.)), 5.);
		assert_eq!(d.eval(Scalar::from(1.)), 3.25);
		assert_eq!(d.eval(Scalar::from(2.)), 2.375);
	}
}