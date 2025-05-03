//! ...

use crate::exp::Exp;
use crate::change::Change;

impl<T: Change + std::ops::Neg> std::ops::Neg for Exp<T> {
	type Output = Exp<<T as std::ops::Neg>::Output>;
	fn neg(self) -> Self::Output {
		Exp(-self.0)
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
		let a_change = crate::change::sum::Sum2 {
			lhs: crate::change::constant::Nil::default(),
			rhs: crate::change::constant::Nil::default(),
			rhs_basis: 1.5,
		};
		let a_basis = 5.;
		let a = a_change.into_poly(a_basis);
		assert_eq!(a.eval(0.), 5.);
		assert_eq!(a.eval(1.), 6.5);
		assert_eq!(a.eval(2.), 8.);
		let b = Test { value: 1., mul: 2. };
		let mut b_change = b.change();
		let b = b.to_poly();
		assert_eq!(b.eval(0.), 1.);
		assert_eq!(b.eval(1.), 2.);
		assert_eq!(b.eval(2.), 4.);
		let c = (a_change >> b_change).into_poly(a_basis);
		assert_eq!(c.eval(0.), 5.);
		assert_eq!(c.eval(1.), 14.32808512266689);
		assert_eq!(c.eval(2.), 32.98425536800067);
		b_change.rhs_basis = 1. / b_change.rhs_basis;
		let d = (a_change >> b_change).into_poly(a_basis);
		assert_eq!(d.eval(0.), 5.);
		assert_eq!(d.eval(1.), 3.0410106403333614);
		assert_eq!(d.eval(2.), 2.0615159605000417);
	}
}