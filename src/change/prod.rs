//! ...

use crate::exp::Exp;
use crate::change::{Change, Sum, SumProd};
use crate::linear::{Basis, Linear};
use crate::poly::SumPoly;

/// Represents the pattern of repeated multiplication, `a = a * b`.
pub type Prod<T, const DEGREE: usize> = Exp<Sum<T, DEGREE>>;

/// ...
pub type ProdPoly<T, const DEGREE: usize> = Exp<SumPoly<T, DEGREE>>;

impl<T: Change + std::ops::Neg> std::ops::Neg for Exp<T> {
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
			mul_term: mul_term.map_inner(Linear::exp),
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
			mul_term: mul_term.map_inner(|x| T::Inner::from_f64(1.).div(x.exp())),
		}
	}
}

impl<A, B> std::ops::Mul<Exp<B>> for Exp<A>
where
	A: Change + std::ops::Add<B>,
	B: Change,
{
	type Output = Exp<A::Output>;
	fn mul(self, rhs: Exp<B>) -> Self::Output {
		Exp(self.0 + rhs.0)
	}
}

impl<A, B> std::ops::Div<Exp<B>> for Exp<A>
where
	A: Change + std::ops::Sub<B>,
	B: Change,
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
	use crate::change::Change;
	use crate::poly::Poly;
	
	#[derive(Flux)]
	pub struct Test {
		#[basis]
		value: f64,
		#[change(mul_per(crate::time::SEC))]
		mul: f64,
	}
	
	#[test]
	fn prod() {
		let a_change = crate::change::sum::Sum([1.5]);
		let a_basis = 5.;
		let a = a_change.into_poly(a_basis);
		assert_eq!(a.eval(0.), 5.);
		assert_eq!(a.eval(1.), 6.5);
		assert_eq!(a.eval(2.), 8.);
		let b = Test { value: 1., mul: 2. };
		let b_change = b.change();
		let b = b.to_poly();
		assert_eq!(b.eval(0.), 1.);
		assert_eq!(b.eval(1.), 2.);
		assert_eq!(b.eval(2.), 4.);
		let c = (a_change * b_change).into_poly(a_basis);
		assert_eq!(c.eval(0.), 5.);
		assert_eq!(c.eval(1.), 13.);
		assert_eq!(c.eval(2.), 29.);
		let d = (a_change / b_change).into_poly(a_basis);
		assert_eq!(d.eval(0.), 5.);
		assert_eq!(d.eval(1.), 3.25);
		assert_eq!(d.eval(2.), 2.375);
	}
}