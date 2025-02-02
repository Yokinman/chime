//! Summation over time.

use std::cmp::Ordering;
use std::ops::{Add, Index, IndexMut, Mul, Neg, Sub};
use crate::{change::*, linear::*};
use crate::change::constant::Constant;
use crate::exp::{ApplyExp, ApplyLn, Ln};
use crate::poly::{Deriv, Poly, Roots, Translate};

/// The pattern of repeated addition, `a = a + b`.
#[derive(Copy, Clone, Debug)]
pub struct Sum<T, const DEGREE: usize>(pub(crate) [T; DEGREE]);

impl<T, const D: usize> Change for Sum<T, D>
where
	T: Basis
{
	type Basis = T;
	type Poly = SumPoly<T, D>;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		let mut i = 0.;
		let mut n = 1.;
		let terms = self.0.map(|term| {
			i += 1.;
			n *= i;
			term.map_inner(|x| x.mul(Linear::from_f64(1. / n)))
		});
		SumPoly::new(basis, terms)
	}
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		Self(self.0.map(|term| term.map_inner(|x| x.mul(scalar.clone()))))
	}
}

impl<T> ChangeUp<'+'> for Sum<T, 0>
where
	T: Basis
{
	type Up = Sum<T, 1>;
	fn up(self, basis: Self::Basis) -> Self::Up {
		Sum([basis])
	}
}

impl<T, const D: usize> Neg for Sum<T, D>
where
	T: Basis
{
	type Output = Self;
	fn neg(self) -> Self::Output {
		Self(self.0.map(|term| term.map_inner(|x| x.mul(Linear::from_f64(-1.)))))
	}
}

/// ...
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Monomial<T, const DEGREE: usize>(pub(crate) T);

impl<T, const D: usize> Poly for Monomial<T, D>
where
	T: Basis
{
	const DEGREE: usize = D;
	type Basis = T;
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		self.0.clone()
			.map_inner(|n| n.mul(time.pow(Linear::from_f64(D as f64))))
	}
	fn zero() -> Self {
		Self(T::zero())
	}
}

impl<T, const D: usize> Deriv for Monomial<T, D>
where
	T: Basis,
	Self: MonomialDown<Down: Deriv<Basis = T>>,
{
	type Deriv = <Self as MonomialDown>::Down;
	fn deriv(mut self) -> Self::Deriv {
		self.0 = self.0.map_inner(|x| x.mul(Linear::from_f64(D as f64)));
		self.downgrade()
	}
}

impl<T, const D: usize> Translate for Monomial<T, D>
where
	T: Basis,
	Constant<T>: MonomialTranslate<D, Basis = T>,
{
	type Output = <Constant<T> as MonomialTranslate<D>>::Output;
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		Constant(self.0).monom_translate(amount)
	}
}

impl<T, const D: usize> MonomialAdd for Monomial<T, D>
where
	T: Basis
{
	type Output = Self;
	fn monom_add(self, rhs: Self) -> Self::Output {
		Self(self.0.zip_map_inner(rhs.0, Linear::add))
	}
}

impl<T, const D: usize, A, P> Add<A> for Monomial<T, D>
where
	Self: AddPoly<A, Output=P>
{
	type Output = P;
	fn add(self, rhs: A) -> Self::Output {
		self.add_poly(rhs)
	}
}

impl<T, const D: usize, A, P> Sub<A> for Monomial<T, D>
where
	A: Neg,
	Self: Add<A::Output, Output=P>,
{
	type Output = P;
	fn sub(self, rhs: A) -> Self::Output {
		self + -rhs
	}
}

impl<T, const D: usize> Neg for Monomial<T, D>
where
	T: Basis
{
	type Output = Self;
	fn neg(self) -> Self::Output {
		Self(self.0.map_inner(|x| x.mul(Linear::from_f64(-1.))))
	}
}

impl<T, const D: usize, A, P> Mul<A> for Monomial<T, D>
where
	Self: ApplyLn<Output: Add<A::Output, Output: ApplyExp<Output=P>>>,
	A: ApplyLn,
{
	type Output = P;
	fn mul(self, rhs: A) -> Self::Output {
		(self.ln() + rhs.ln()).exp()
	}
}

impl<T, const D: usize> ApplyLn for Monomial<T, D> {
	type Output = Ln<Self>;
	fn ln(self) -> Self::Output {
		Ln(self)
	}
}

impl<const D: usize> Roots for Monomial<f64, D> {
	type Output = [f64; D];
	fn roots(self) -> <Self as Roots>::Output {
		[0. / self.0; D]
	}
}

impl<const D: usize> Roots for Binomial<Constant<f64>, Monomial<f64, D>> {
	type Output = [f64; D];
	fn roots(self) -> <Self as Roots>::Output {
		let Binomial { lhs: Constant(a), rhs: Monomial(b) } = self;
		
		if D == 1 {
			return [-a / b; D]
		}
		
		if a == 0. {
			return Monomial::<_, D>(b).roots()
		}
		
		let mut roots = [f64::NAN; D];
		roots[0] = (-a / b).powf((D as f64).recip());
		if D % 2 == 0 {
			roots[1] = -roots[0];
		}
		roots
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Monomial<f64, 2>>> {
	type Output = [f64; 2];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(b), rhs: Monomial(c) },
		} = self;
		
		let mut n = -b / c;
		
		 // Precision Breakdown (Discriminant Quotient):
		let x = -a / b;
		if x.is_nan() || n.is_nan() || (x / n).abs() < 1e-8 {
			return if x.is_nan() {
				[n, n]
			} else {
				[x, n-x]
			}
		}
		
		 // General Quadratic:
		n /= 2.;
		let x = n.mul_add(n, -a/c).sqrt();
		[n-x, n+x]
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Monomial<f64, 3>>> {
	type Output = [f64; 3];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(b), rhs: Monomial(d) }
		} = self;
		
		 // Pseudo-linear:
		if b == 0. {
			return Binomial { lhs: Constant(a), rhs: Monomial::<_, 3>(d) }.roots()
		}
		
		 // Depressed Cubic:
		let p = -a / (2. * d);
		let q = -b / (3. * d);
		let r = p / q;
		let discriminant = if q < 0. {
			1.
		} else {
			let disc = f64::ln((r*r) / q.abs());
			if disc.abs() < 1e-15 {
				0.
			} else {
				disc
			}
		};
		match discriminant.partial_cmp(&0.) {
			 // 3 Real Roots:
			Some(Ordering::Less) => {
				let sqrt_q = q.sqrt();
				let angle = f64::acos(r / sqrt_q);
				debug_assert!(!angle.is_nan(), "{:?}", self);
				use std::f64::consts::TAU;
				[
					2. * sqrt_q * f64::cos((angle         ) / 3.),
					2. * sqrt_q * f64::cos((angle +    TAU) / 3.),
					2. * sqrt_q * f64::cos((angle + 2.*TAU) / 3.),
				] 
			},
			
			 // 1 Real Root:
			Some(Ordering::Greater) => {
				let n = q.mul_add(-q*q, p*p).sqrt();
				debug_assert!(!n.is_nan());
				[(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN]
			},
			
			 // 2 Real Roots:
			_ => [-r, -r, 2.*r]
		}
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 2>, Monomial<f64, 3>>> {
	type Output = [f64; 3];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(c), rhs: Monomial(d) }
		} = self;
		Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(0.),
				rhs: Binomial { lhs: Monomial::<_, 2>(c), rhs: Monomial::<_, 3>(d) }
			}
		}.roots()
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Binomial<Monomial<f64, 2>, Monomial<f64, 3>>>> {
	type Output = [f64; 3];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial(b),
				rhs: Binomial { lhs: Monomial(c), rhs: Monomial(d) }
			}
		} = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x   * c) / b).abs() < 1e-5 && // ??? Adjust as needed
			((x*x * d) / b).abs() < 1e-16
		) {
			let [y, z] = Binomial {
				lhs: Constant(b),
				rhs: Binomial { lhs: Monomial::<_, 1>(c), rhs: Monomial::<_, 2>(d) },
			}.roots();
			return if x.is_nan() {
				[y, y, z]
			} else {
				[x, y-x, z]
			}
		}
		
		let mut n = -c / d;
		
		 // Weak Leading Term:
		if n.is_nan() || (
			(b / (n   * c)).abs() < 1e-9 && // ??? Adjust as needed
			(a / (n*n * c)).abs() < 1e-13
		) {
			let [x, y] = Binomial {
				lhs: Constant(a),
				rhs: Binomial { lhs: Monomial::<_, 1>(b), rhs: Monomial::<_, 2>(c) }
			}.roots();
			return [x, y, n]
		}
		
		 // Depressed Cubic:
		if c == 0. {
			return Binomial {
				lhs: Constant(a),
				rhs: Binomial { lhs: Monomial::<_, 1>(b), rhs: Monomial::<_, 3>(d) }
			}.roots()
		}
		
		 // General Cubic:
		n /= 3.;
		let depressed_cubic = Binomial {
			lhs: Constant(n.mul_add(n.mul_add(c * (2./3.), b), a)),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(n.mul_add(c, b)),
				rhs: Monomial::<_, 3>(d),
			}
		};
		depressed_cubic.roots().map(|x| x + n)
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Monomial<f64, 4>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(b), rhs: Monomial(e) }
		} = self;
		Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(b),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(0.),
					rhs: Binomial { lhs: Monomial::<_, 3>(0.), rhs: Monomial::<_, 4>(e) }
				}
			}
		}.roots()
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 2>, Monomial<f64, 4>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(c), rhs: Monomial(e) }
		} = self;
		
		 // Pseudo-linear:
		if c == 0. {
			return Binomial { lhs: Constant(a), rhs: Monomial::<_, 4>(e) }.roots()
		}
		
		 // Biquadratic:
		let [x, y] = Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial::<_, 1>(c), rhs: Monomial::<_, 2>(e) }
		}.roots();
		let (x, y) = (x.sqrt(), y.sqrt());
		return [-x, x, -y, y];
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 3>, Monomial<f64, 4>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial { lhs: Monomial(d), rhs: Monomial(e) }
		} = self;
		Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(0.),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(0.),
					rhs: Binomial { lhs: Monomial::<_, 3>(d), rhs: Monomial::<_, 4>(e) }
				}
			}
		}.roots()
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Binomial<Monomial<f64, 2>, Monomial<f64, 4>>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial(b),
				rhs: Binomial { lhs: Monomial(c), rhs: Monomial(e) }
			},
		} = self;
		
		 // Biquadratic:
		if b == 0. {
			return Binomial {
				lhs: Constant(a),
				rhs: Binomial { lhs: Monomial::<_, 2>(c), rhs: Monomial::<_, 4>(e) }
			}.roots()
		}
		
		 // Depressed Quartic:
		let p = a / e;
		let q = b / (2. * e);
		let r = c / (2. * e);
		let resolvent_cubic = Binomial {
			lhs: Constant(-(q * q) / 2.),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(r.mul_add(r, -p)),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(2. * r),
					rhs: Monomial::<_, 3>(1.),
				}
			}
		};
		let mut m = resolvent_cubic.roots().into_iter()
			.find(|&r| r > 0. && r.is_finite())
			.unwrap_or_else(|| panic!(
				"expected a positive root from: {:?}",
				resolvent_cubic
			));
		let y = resolvent_cubic.eval(m)
			/ m.mul_add(m.mul_add(3., 4.*r), resolvent_cubic.rhs.lhs.0);
		if m > y && y.is_finite() {
			m -= y; // Newton-Raphson step
		}
		let sqrt_2m = (2. * m).sqrt();
		let [x, y] = Binomial {
			lhs: Constant((q / sqrt_2m) + r + m),
			rhs: Binomial { lhs: Monomial::<_, 1>(-sqrt_2m), rhs: Monomial::<_, 2>(1.) }
		}.roots();
		let [z, w] = Binomial {
			lhs: Constant(-(q / sqrt_2m) + r + m),
			rhs: Binomial { lhs: Monomial::<_, 1>(sqrt_2m), rhs: Monomial::<_, 2>(1.) }
		}.roots();
		[x, y, z, w]
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Binomial<Monomial<f64, 3>, Monomial<f64, 4>>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial(b),
				rhs: Binomial { lhs: Monomial(d), rhs: Monomial(e) }
			}
		} = self;
		Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(b),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(0.),
					rhs: Binomial { lhs: Monomial::<_, 3>(d), rhs: Monomial::<_, 4>(e) }
				}
			}
		}.roots()
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 2>, Binomial<Monomial<f64, 3>, Monomial<f64, 4>>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial(c),
				rhs: Binomial { lhs: Monomial(d), rhs: Monomial(e) }
			}
		} = self;
		Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(0.),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(c),
					rhs: Binomial { lhs: Monomial::<_, 3>(d), rhs: Monomial::<_, 4>(e) }
				}
			}
		}.roots()
	}
}

impl Roots for Binomial<Constant<f64>, Binomial<Monomial<f64, 1>, Binomial<Monomial<f64, 2>, Binomial<Monomial<f64, 3>, Monomial<f64, 4>>>>> {
	type Output = [f64; 4];
	fn roots(self) -> Self::Output {
		let Binomial {
			lhs: Constant(a),
			rhs: Binomial {
				lhs: Monomial(b),
				rhs: Binomial {
					lhs: Monomial(c),
					rhs: Binomial { lhs: Monomial(d), rhs: Monomial(e) }
				}
			}
		} = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x     * c) / b).abs() < 1e-4  && // ??? Adjust as needed
			((x*x   * d) / b).abs() < 1e-12 &&
			((x*x*x * e) / b).abs() < 1e-20
		) {
			let [y, z, w] = Binomial {
				lhs: Constant(b),
				rhs: Binomial {
					lhs: Monomial::<_, 1>(c),
					rhs: Binomial { lhs: Monomial::<_, 2>(d), rhs: Monomial::<_, 3>(e) }
				}
			}.roots();
			return if x.is_nan() {
				[y, y, z, w]
			} else {
				[x, y-x, z, w]
			}
		}
		
		let mut n = -d / e;
		
		 // Weak Leading Term:
		if n.is_nan() || (
		    (c / (n     * d)).abs() < 1e-7  && // ??? Adjust as needed
			(b / (n*n   * d)).abs() < 1e-10 &&
			(a / (n*n*n * d)).abs() < 1e-15
		) {
			let [x, y, z] = Binomial {
				lhs: Constant(a),
				rhs: Binomial {
					lhs: Monomial::<_, 1>(b),
					rhs: Binomial { lhs: Monomial::<_, 2>(c), rhs: Monomial::<_, 3>(d) }
				}
			}.roots();
			return [x, y, z, n]
		}
		
		 // Depressed Quartic:
		if d == 0. {
			return Binomial {
				lhs: Constant(a),
				rhs: Binomial {
					lhs: Monomial::<_, 1>(b),
					rhs: Binomial { lhs: Monomial::<_, 2>(c), rhs: Monomial::<_, 4>(e) }
				}
			}.roots()
		}
		
		 // General Quartic:
		n /= 4.;
		let depressed_quartic = Binomial {
			lhs: Constant(n.mul_add(n.mul_add(n.mul_add(d * (3./4.), c), b), a)),
			rhs: Binomial {
				lhs: Monomial::<_, 1>(n.mul_add(n.mul_add(d, c) * 2., b)),
				rhs: Binomial {
					lhs: Monomial::<_, 2>(n.mul_add(d * (3./2.), c)),
					rhs: Monomial::<_, 4>(e),
				}
			}
		};
		let [x, y, z, w] = depressed_quartic.roots();
		[x+n, y+n, z+n, w+n]
	}
}

impl<T, const A: usize, const B: usize, const N: usize> Roots for Binomial<Monomial<f64, A>, T>
where
	Self: Poly + BinomialDown<Down: Roots<Output = [f64; N]>>,
	T: Roots<Output = [f64; B]>,
{
	type Output = [f64; B];
	fn roots(self) -> <Self as Roots>::Output {
		let mut roots = [0.; B];
		roots[1..=N].copy_from_slice(&self.binom_downgrade().roots());
		roots
	}
}

/// ...
pub trait AddPoly<T> {
	type Output;
	fn add_poly(self, rhs: T) -> Self::Output;
}

impl<A, B, O, P> AddPoly<B> for A
where
	A: MonomialCmp<B, Order=O> + MonomialAdd<B, O, Output=P>,
{
	type Output = P;
	fn add_poly(self, rhs: B) -> Self::Output {
		self.monom_add(rhs)
	}
}

/// ...
pub trait MonomialAdd<Rhs=Self, OrdMarker=order::Same> {
	type Output;
	fn monom_add(self, rhs: Rhs) -> Self::Output;
}

impl<A, B> MonomialAdd<B, order::Below> for A {
	type Output = Binomial<A, B>;
	fn monom_add(self, rhs: B) -> Self::Output {
		Binomial {
			lhs: self,
			rhs,
		}
	}
}

impl<A, B> MonomialAdd<B, order::Above> for A {
	type Output = Binomial<B, A>;
	fn monom_add(self, rhs: B) -> Self::Output {
		Binomial {
			lhs: rhs,
			rhs: self,
		}
	}
}

impl<A, B, T> MonomialAdd<Binomial<A, B>, order::ToLhs<order::Below>> for T {
	type Output = Binomial<T, Binomial<A, B>>;
	fn monom_add(self, rhs: Binomial<A, B>) -> Self::Output {
		Binomial {
			lhs: self,
			rhs,
		}
	}
}

impl<A, B, T> MonomialAdd<Binomial<A, B>, order::ToLhs<order::Same>> for T
where
	T: MonomialAdd<A, order::Same>
{
	type Output = Binomial<T::Output, B>;
	fn monom_add(self, rhs: Binomial<A, B>) -> Self::Output {
		Binomial {
			lhs: self.monom_add(rhs.lhs),
			rhs: rhs.rhs,
		}
	}
}

impl<A, B, T> MonomialAdd<Binomial<A, B>, order::ToLhs<order::Above>> for T
where
	T: MonomialCmp<B> + MonomialAdd<B, T::Order>,
{
	type Output = Binomial<A, T::Output>;
	fn monom_add(self, rhs: Binomial<A, B>) -> Self::Output {
		Binomial {
			lhs: rhs.lhs,
			rhs: self.monom_add(rhs.rhs),
		}
	}
}

impl<A, B, T, O> MonomialAdd<T, order::FromLhs<O>> for Binomial<A, B>
where
	A: MonomialAdd<T, O>,
	B: MonomialCmp<A::Output> + MonomialAdd<A::Output, B::Order>,
{
	type Output = B::Output;
	fn monom_add(self, rhs: T) -> Self::Output {
		self.rhs.monom_add(self.lhs.monom_add(rhs))
	}
}

impl<T, U, V, const A: usize, const B: usize> MonomialAdd<Ln<Monomial<T, B>>, order::Same> for Ln<Monomial<T, A>>
where
	Monomial<T, A>: MonomialDown<Down = V>,
	Monomial<T, B>: MonomialUp<Up = U>,
	Ln<V>: Add<Ln<U>>,
{
	type Output = <Ln<V> as Add<Ln<U>>>::Output;
	fn monom_add(self, rhs: Ln<Monomial<T, B>>) -> Self::Output {
		Ln(self.0.downgrade()) + Ln(rhs.0.upgrade())
	}
}

impl<T, const D: usize> MonomialAdd<Ln<Monomial<T, D>>, order::Same> for Ln<Constant<T>>
where
	T: Basis
{
	type Output = Ln<Monomial<T, D>>;
	fn monom_add(self, rhs: Ln<Monomial<T, D>>) -> Self::Output {
		Ln(Monomial(self.0.0.zip_map_inner(rhs.0.0, Linear::mul)))
	}
}

impl<T, const D: usize> MonomialAdd<Ln<Constant<T>>, order::Same> for Ln<Monomial<T, D>>
where
	T: Basis
{
	type Output = Ln<Monomial<T, D>>;
	fn monom_add(self, rhs: Ln<Constant<T>>) -> Self::Output {
		Ln(Monomial(self.0.0.zip_map_inner(rhs.0.0, Linear::mul)))
	}
}

impl<T, const D: usize, A, B, P> MonomialAdd<Ln<Binomial<A, B>>, order::Same> for Ln<Monomial<T, D>>
where
	Binomial<A, B>: Mul<Monomial<T, D>, Output = P>,
{
	type Output = Ln<P>;
	fn monom_add(self, rhs: Ln<Binomial<A, B>>) -> Self::Output {
		Ln(rhs.0 * self.0)
	}
}

/// ...
pub trait MonomialMul<Rhs=Self> {
	type Output;
	fn monom_mul(self, rhs: Rhs) -> Self::Output;
}

impl<T, U, V, const A: usize, const B: usize> MonomialMul<Monomial<T, B>> for Monomial<T, A>
where
	Monomial<T, A>: MonomialDown<Down = V>,
	Monomial<T, B>: MonomialUp<Up = U>,
	V: MonomialMul<U>,
{
	type Output = V::Output;
	fn monom_mul(self, rhs: Monomial<T, B>) -> Self::Output {
		self.downgrade().monom_mul(rhs.upgrade())
	}
}

impl<T, const D: usize> MonomialMul<Monomial<T, D>> for Constant<T>
where
	T: Basis
{
	type Output = Monomial<T, D>;
	fn monom_mul(self, rhs: Monomial<T, D>) -> Self::Output {
		Monomial(self.0.zip_map_inner(rhs.0, Linear::mul))
	}
}

impl<T, const D: usize> MonomialMul<Constant<T>> for Monomial<T, D>
where
	T: Basis
{
	type Output = Monomial<T, D>;
	fn monom_mul(self, rhs: Constant<T>) -> Self::Output {
		Monomial(self.0.zip_map_inner(rhs.0, Linear::mul))
	}
}

/// ...
pub trait MonomialTranslate<const N: usize, OrdMarker = order::Below>: Poly {
	type Output: Poly<Basis = Self::Basis>;
	fn monom_translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output;
}

impl<T, O, P, const N: usize> MonomialTranslate<N> for Constant<T>
where
	T: Basis,
	P: Poly,
	Monomial<T, 1>: MonomialCmp<Monomial<T, N>, Order=O> + MonomialTranslate<N, O, Output=P, Basis=T>,
	Self: MonomialOrder<P, Basis = T>,
{
	type Output = Binomial<Self, P>;
	fn monom_translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		let fac = amount.mul(Linear::from_f64(-1.))
			.pow(Linear::from_f64(N as f64));
		Binomial {
			lhs: Self(self.0.clone().map_inner(|n| n.mul(fac))),
			rhs: self.upgrade().monom_translate(amount),
		}
	}
}

impl<T, const N: usize> MonomialTranslate<N, order::Same> for Monomial<T, N>
where
	T: Basis
{
	type Output = Self;
	fn monom_translate(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		self
	}
}

impl<T, O, P, const D: usize, const N: usize> MonomialTranslate<N> for Monomial<T, D>
where
	T: Basis,
	P: Poly,
	Self: MonomialUp<Up: MonomialCmp<Monomial<T, N>, Order=O> + MonomialTranslate<N, O, Output=P, Basis=T>>
		+ MonomialOrder<P, Basis = T>,
{
	type Output = Binomial<Self, P>;
	fn monom_translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		// binomial coefficient
		const fn binom(n: usize, k: usize) -> usize {
		    if n < k {
		        return 0
		    }
		    let mut numer = 1;
		    let mut denom = 1;
		    let mut n = n;
		    let mut k = k;
		    while k != 0 {
		        numer *= n;
		        denom *= k;
		        n -= 1;
		        k -= 1;
		    }
		    numer / denom
		}
		
		let fac = amount.mul(Linear::from_f64(-1.))
			.pow(Linear::from_f64((N - D) as f64))
			.mul(Linear::from_f64(binom(N, D) as f64));
		Binomial {
			lhs: Monomial(self.0.clone().map_inner(|n| n.mul(fac))),
			rhs: self.upgrade().monom_translate(amount),
		}
	}
}

/// ...
pub trait MonomialUp {
	type Up;
	fn upgrade(self) -> Self::Up;
}

impl<T> MonomialUp for Constant<T> {
	type Up = Monomial<T, 1>;
	fn upgrade(self) -> Self::Up {
		Monomial(self.0)
	}
}

/// ...
pub trait MonomialDown {
	type Down;
	fn downgrade(self) -> Self::Down;
}

impl<T> MonomialDown for Monomial<T, 1> {
	type Down = Constant<T>;
	fn downgrade(self) -> Self::Down {
		Constant(self.0)
	}
}

/// ...
pub trait BinomialDown {
	type Down;
	fn binom_downgrade(self) -> Self::Down;
}

impl<T> BinomialDown for T
where
	T: MonomialDown
{
	type Down = T::Down;
	fn binom_downgrade(self) -> Self::Down {
		self.downgrade()
	}
}

impl<A, B> BinomialDown for Binomial<A, B>
where
	A: MonomialDown,
	B: BinomialDown,
{
	type Down = Binomial<A::Down, B::Down>;
	fn binom_downgrade(self) -> Self::Down {
		Binomial {
			lhs: self.lhs.downgrade(),
			rhs: self.rhs.binom_downgrade(),
		}
	}
}

/// ...
mod order {
	pub struct Above;
	pub struct Below;
	pub struct Same;
	pub struct ToLhs<T>(T);
	pub struct FromLhs<T>(T);
}

/// ...
pub trait MonomialCmp<T> {
	type Order;
}

impl<A, B, O> MonomialCmp<B> for A
where
	A: MonomialDown,
	B: MonomialDown,
	A::Down: MonomialCmp<B::Down, Order = O>,
{
	type Order = O;
}

impl<T, B> MonomialCmp<Constant<B>> for T
where
	T: MonomialDown
{
	type Order = order::Above;
}

impl<T, B> MonomialCmp<T> for Constant<B>
where
	T: MonomialDown
{
	type Order = order::Below;
}

impl<B> MonomialCmp<Constant<B>> for Constant<B> {
	type Order = order::Same;
}

impl<A, B, T, O> MonomialCmp<Binomial<A, B>> for T
where
	A: MonomialDown,
	T: MonomialDown,
	T::Down: MonomialCmp<A::Down, Order=O>,
{
	type Order = order::ToLhs<O>;
}

impl<A, B, T> MonomialCmp<Binomial<Constant<A>, B>> for T
where
	T: MonomialDown,
{
	type Order = order::ToLhs<order::Above>;
}

impl<A, B, T> MonomialCmp<Binomial<A, B>> for Constant<T>
where
	A: MonomialDown,
{
	type Order = order::ToLhs<order::Below>;
}

impl<A, B, T> MonomialCmp<Binomial<Constant<A>, B>> for Constant<T> {
	type Order = order::ToLhs<order::Same>;
}

impl<A, B, T> MonomialCmp<T> for Binomial<A, B>
where
	A: MonomialCmp<T>,
{
	type Order = order::FromLhs<A::Order>;
}

impl<T, U> MonomialCmp<Ln<U>> for Ln<T> {
	type Order = order::Same;
}

/// ...
pub trait MonomialOrder<T>: Poly {
	fn binom_eval(self, rhs: T, time: <Self::Basis as Basis>::Inner) -> Self;
}

impl<T, const D: usize> MonomialOrder<Monomial<T, D>> for Constant<T>
where
	T: Basis,
{
	// a + bx^B
	fn binom_eval(
		self,
		rhs: Monomial<T, D>,
		time: <Self::Basis as Basis>::Inner,
	) -> Self {
		Self(self.0.zip_map_inner(rhs.eval(time), Linear::add))
	}
}

impl<T, const A: usize, const B: usize> MonomialOrder<Monomial<T, B>> for Monomial<T, A>
where
	T: Basis,
	Self: MonomialCmp<Monomial<T, B>, Order = order::Below>,
{
	// ax^A + bx^B -> (a + bx^{B-A})x^A
	fn binom_eval(
		self,
		rhs: Monomial<T, B>,
		time: <Self::Basis as Basis>::Inner,
	) -> Self {
		let fac = time.pow(Linear::from_f64((B - A) as f64));
		Self(self.0.zip_map_inner(
			rhs.0,
			|a, b| a.add(b.mul(fac))
		))
	}
}

impl<T, A, B> MonomialOrder<Binomial<A, B>> for T
where
	A: MonomialOrder<B>,
	B: Poly,
	T: MonomialOrder<A, Basis = A::Basis>,
{
	// n + (ax^A + bx^B) -> n + (a + bx^{B-A})x^A
	fn binom_eval(
		self,
		rhs: Binomial<A, B>,
		time: <Self::Basis as Basis>::Inner,
	) -> Self {
		self.binom_eval(rhs.lhs.binom_eval(rhs.rhs, time), time)
	}
}

/// Summation over time.
/// 
/// `Binomial<Monomial<T, 2>, Binomial<Monomial<T, 1>, Constant<T>>>`: `a + bx + cx^2`
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Binomial<A, B> {
	pub(crate) lhs: A,
	pub(crate) rhs: B,
}

impl<A, B> Poly for Binomial<A, B>
where
	A: MonomialOrder<B>,
	B: Poly,
{
	const DEGREE: usize = B::DEGREE;
	type Basis = A::Basis;
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		self.lhs.clone()
			.binom_eval(self.rhs.clone(), time)
			.eval(time)
	}
	fn zero() -> Self {
		Self {
			lhs: A::zero(),
			rhs: B::zero(),
		}
	}
}

impl<T, const D: usize> Deriv for Binomial<Monomial<T::Basis, D>, T>
where
	T: Deriv,
	Monomial<T::Basis, D>: MonomialOrder<T, Basis = T::Basis> + Deriv,
	Binomial<<Monomial<T::Basis, D> as Deriv>::Deriv, T::Deriv>: Deriv<Basis = T::Basis>,
{
	type Deriv = Binomial<<Monomial<T::Basis, D> as Deriv>::Deriv, T::Deriv>;
	fn deriv(self) -> Self::Deriv {
		Binomial {
			lhs: self.lhs.deriv(),
			rhs: self.rhs.deriv(),
		}
	}
}

impl<T> Deriv for Binomial<Constant<T::Basis>, T>
where
	T: Deriv,
	Constant<T::Basis>: MonomialOrder<T, Basis = T::Basis>,
{
	type Deriv = T::Deriv;
	fn deriv(self) -> Self::Deriv {
		self.rhs.deriv()
	}
}

impl<T, A, B, O> Translate for Binomial<A, B>
where
	Self: Poly<Basis=T>,
	T: Basis,
	A: Translate<Basis = T>,
	B: Translate<Basis = T>,
	A::Output: Add<B::Output, Output=O>,
	O: Poly<Basis = T>,
{
	type Output = O;
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		self.lhs.translate(amount) + self.rhs.translate(amount)
	}
}

impl<A, B, C, P> Add<C> for Binomial<A, B>
where
	Self: AddPoly<C, Output=P>
{
	type Output = P;
	fn add(self, rhs: C) -> Self::Output {
		self.add_poly(rhs)
	}
}

impl<A, B, C, P> Sub<C> for Binomial<A, B>
where
	C: Neg,
	Self: Add<C::Output, Output=P>,
{
	type Output = P;
	fn sub(self, rhs: C) -> Self::Output {
		self + -rhs
	}
}

impl<A, B> Neg for Binomial<A, B>
where
	A: Neg,
	B: Neg,
{
	type Output = Binomial<A::Output, B::Output>;
	fn neg(self) -> Self::Output {
		Binomial {
			lhs: -self.lhs,
			rhs: -self.rhs,
		}
	}
}

impl<A, B, C, P> Mul<C> for Binomial<A, B>
where
	A: Mul<C, Output: Add<B::Output, Output=P>>,
	B: Mul<C>,
	C: Clone,
{
	type Output = P;
	fn mul(self, rhs: C) -> Self::Output {
		(self.lhs * rhs.clone()) + (self.rhs * rhs)
	}
}

impl<A, B> ApplyLn for Binomial<A, B> {
	type Output = Ln<Self>;
	fn ln(self) -> Self::Output {
		Ln(self)
	}
}

#[test]
fn binomial_temp() {
	let a = SumPoly::new(5., [10., 0.9, 0., 4.]);
	let b = (Constant(5.) + Monomial::<_, 1>(5.) + Monomial::<_, 1>(5.)) + (Monomial::<_, 4>(4.) - Monomial::<_, 2>(-0.9));
	for i in 0..100 {
		let t = i as f64;
		println!("> {:?}", (a.eval(t), b.eval(t)));
		assert_eq!(a.eval(t), b.eval(t));
		assert_eq!(a.deriv().eval(t), b.deriv().eval(t));
	}
	
	assert_eq!(
		(Monomial::<_, 1>(10.) + Constant(-2.) + Monomial::<_, 4>(4.) + Monomial::<_, 3>(-32.)).roots(),
		SumPoly(-2., [10., 0., -32., 4.]).roots()
	);
	assert!(a.roots().iter().zip(b.roots())
		.all(|(x, y)| x.total_cmp(&y) == Ordering::Equal));
	
	let Binomial {
		lhs: Constant(52.599999999999994),
		rhs: Binomial {
			lhs: Monomial::<_, 1>(-121.6),
			rhs: Binomial {
				lhs: Monomial::<_, 2>(96.9),
				rhs: Binomial {
					lhs: Monomial::<_, 3>(-32.),
					rhs: Monomial::<_, 4>(4.),
				}
			},
		}
	} = b.translate(2.) else {
		panic!("> {:?}", b.translate(2.));
	};
	
	assert_eq!(
		(Monomial::<_, 5>(4.) + Monomial::<_, 3>(5.)) * (Monomial::<_, 2>(2.) + Monomial::<_, 3>(4.)),
		Binomial {
			lhs: Monomial::<f64, 5>(10.),
			rhs: Binomial {
				lhs: Monomial::<_, 6>(20.),
				rhs: Binomial {
					lhs: Monomial::<_, 7>(8.),
					rhs: Monomial::<_, 8>(16.),
				}
			},
		}
	);
}

/// Summation over time.
/// 
/// - `SumPoly<0>`: `a`
/// - `SumPoly<1>`: `a + bx`
/// - `SumPoly<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct SumPoly<T, const DEGREE: usize>(pub(crate) T, pub(crate) [T; DEGREE]);

impl<T, const D: usize> SumPoly<T, D> {
	pub fn new(basis: T, terms: [T; D]) -> Self {
		Self(basis, terms)
	}
	
	pub fn map<U>(self, mut f: impl FnMut(T) -> U) -> SumPoly<U, D> {
		SumPoly(f(self.0), self.1.map(f))
	}
	
	pub fn iter(&self) -> impl Iterator<Item = &T> {
		std::iter::once(&self.0).chain(&self.1)
	}
	
	pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
		std::iter::once(&mut self.0).chain(&mut self.1)
	}
	
	#[allow(clippy::should_implement_trait)]
	pub fn into_iter(self) -> impl Iterator<Item = T> {
		// !!! impl IntoIterator
		std::iter::once(self.0).chain(self.1)
	}
}

impl<T: Basis, const D: usize> From<T> for SumPoly<T, D> {
	fn from(value: T) -> Self {
		Self(value, std::array::from_fn(|_| T::zero()))
	}
}

impl<T, const DEG: usize> Index<usize> for SumPoly<T, DEG> {
	type Output = T;
	fn index(&self, index: usize) -> &Self::Output {
		if index == 0 {
			&self.0
		} else {
			&self.1[index - 1]
		}
	}
}

impl<T, const D: usize> IndexMut<usize> for SumPoly<T, D> {
	fn index_mut(&mut self, index: usize) -> &mut Self::Output {
		if index == 0 {
			&mut self.0
		} else {
			&mut self.1[index - 1]
		}
	}
}

impl<T, const D: usize> PartialEq for SumPoly<T, D>
where
	T: PartialEq,
{
	fn eq(&self, other: &Self) -> bool {
		if self.0 != other.0 {
			return false
		}
		for i in 0..D {
			if self.1[i] != other.1[i] {
				return false
			}
		}
		true
	}
}

impl<T: Basis, const D: usize> Neg for SumPoly<T, D> {
	type Output = Self;
	fn neg(mut self) -> Self::Output {
		self.0 = self.0.map_inner(|x| x.mul(Linear::from_f64(-1.)));
		self.1 = self.1.map(|term| term.map_inner(|x| x.mul(Linear::from_f64(-1.))));
		self
	}
}

impl<T: Basis, const D: usize> Poly for SumPoly<T, D> {
	const DEGREE: usize = D;
	
	type Basis = T;
	
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
		if time == Linear::zero() {
			return self.0.clone()
		}
		let mut value = T::zero();
		for degree in 1..=D {
			value = self.1[D - degree].clone()
				.zip_map_inner(value, |a, b| a.add(b.mul(time)))
		}
		self.0.clone()
			.zip_map_inner(value, |a, b| a.add(b.mul(time)))
	}

	fn zero() -> Self {
		Self(T::zero(), std::array::from_fn(|_| T::zero()))
	}
}

impl<T: Basis> Deriv for SumPoly<T, 0> {
	type Deriv = Self;
	fn deriv(self) -> Self::Deriv {
		Self(T::zero(), [])
	}
}

impl<T: Basis> Translate for SumPoly<T, 0> {
	type Output = Self;
	fn translate(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		self
	}
}

impl<T, const DEGREE: usize, const SIZE: usize> Vector<SIZE> for SumPoly<T, DEGREE>
where
	T: Vector<SIZE>
{
	type Output = SumPoly<T::Output, DEGREE>;
	fn index(&self, index: usize) -> Self::Output {
		SumPoly(
			self.0.index(index),
			std::array::from_fn(|i| self.1[i].index(index))
		)
	}
}

#[allow(clippy::from_over_into)]
impl<T: Basis> Into<Constant<T>> for SumPoly<T, 0> {
	fn into(self) -> Constant<T> {
		Constant::from(self.0)
	}
}

impl<T: Basis> From<Constant<T>> for SumPoly<T, 0> {
	fn from(Constant(value): Constant<T>) -> Self {
		Self::from(value)
	}
}

impl<A, B, const D: usize> Add<B> for SumPoly<A, D>
where
	A: Basis,
	B: Poly + Into<Constant<A>>,
{
	type Output = Self;
	fn add(mut self, rhs: B) -> Self {
		self.0 = self.0.zip_map_inner(rhs.into().0, Linear::add);
		self
	}
}

/// Degree sequential ordering.
macro_rules! impl_deg_order {
	(1  1  $($num:tt)*) => { impl_deg_order!(2  $($num)*); };
	(2  2  $($num:tt)*) => { impl_deg_order!(4  $($num)*); };
	(4  4  $($num:tt)*) => { impl_deg_order!(8  $($num)*); };
	// (8  8  $($num:tt)*) => { impl_deg_order!(16 $($num)*); }; // !!! This is really slow
	// (16 16 $($num:tt)*) => { impl_deg_order!(32 $($num)*); };
	// (32 32 $($num:tt)*) => { impl_deg_order!(64 $($num)*); };
	(8) => {/* break */};
	($($num:tt)+) => {
		impl<T> MonomialUp for Monomial<T, { $($num +)+ 0 }> {
			type Up = Monomial<T, { $($num +)+ 1 }>;
			fn upgrade(self) -> Self::Up {
				Monomial(self.0)
			}
		}
		
		impl<T> MonomialDown for Monomial<T, { $($num +)+ 1 }> {
			type Down = Monomial<T, { $($num +)+ 0 }>;
			fn downgrade(self) -> Self::Down {
				Monomial(self.0)
			}
		}
		
		impl<T: Basis> Deriv for SumPoly<T, { $($num +)+ 0 }> {
			type Deriv = SumPoly<T, { $($num +)+ 0 - 1 }>;
			fn deriv(self) -> Self::Deriv {
				SumPoly(
					self.1[0].clone(),
					std::array::from_fn(|i| {
						let term = self.1[i + 1].clone();
						term.map_inner(|x| x.mul(Linear::from_f64((i + 2) as f64)))
					})
				)
			}
		}
		impl<T: Basis> Translate for SumPoly<T, { $($num +)+ 0 }> {
			type Output = Self;
			fn translate(mut self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
				if amount == Linear::zero() {
					return self;
				}
				let mut deriv = self.clone();
				self.0 = deriv.eval(amount);
				const D: usize = { $($num +)+ 0 };
				for degree in 1..D {
					deriv = {
						std::mem::swap(&mut deriv.0, &mut deriv.1[0]);
						deriv.1.rotate_left(1);
						deriv.1[D-1] = T::zero();
						let mut d = 1.;
						deriv.1 = deriv.1.map(|term| {
							d += 1.;
							term.map_inner(|x| x.mul(Linear::from_f64(d)))
						});
						deriv
					};
					deriv = deriv.map(|term| term.map_inner(|x| {
						x.mul(Linear::from_f64(1. / (degree as f64)))
					}));
					self.1[degree-1] = deriv.eval(amount);
				}
				self
			}
		}
		impl<T: Basis> ChangeUp<'+'> for Sum<T, { $($num +)+ 0 }> {
			type Up = Sum<T, { $($num +)+ 1 }>;
			fn up(self, basis: Self::Basis) -> Self::Up {
				Sum(unsafe {
					// SAFETY: This seems to work. Maybe I'll validate this
					// later, lol.
					let mut terms = std::mem::MaybeUninit::uninit();
					let ptr = terms.as_mut_ptr() as *mut T;
					ptr.write(basis);
					let ptr = ptr.add(1) as *mut [T; { $($num +)+ 0 }];
					ptr.write(self.0);
					terms.assume_init()
				})
			}
		}
		impl<T: Basis> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, 0> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				rhs
			}
		}
		impl<T: Basis> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.0.into_iter();
				let mut b = rhs.0.into_iter();
				Sum(std::array::from_fn(|_| unsafe {
					// SAFETY: Sizes of all input & output arrays are equal.
					a.next().unwrap_unchecked()
						.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
				}))
			}
		}
		impl<T: Basis> Add<SumPoly<T, { $($num +)+ 0 }>> for SumPoly<T, 0> {
			type Output = SumPoly<T, { $($num +)+ 0 }>;
			fn add(self, rhs: SumPoly<T, { $($num +)+ 0 }>) -> Self::Output {
				rhs
			}
		}
		impl<T: Basis> Add<SumPoly<T, { $($num +)+ 0 }>> for SumPoly<T, { $($num +)+ 0 }> {
			type Output = SumPoly<T, { $($num +)+ 0 }>;
			fn add(self, rhs: SumPoly<T, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				SumPoly(
					self.0.zip_map_inner(rhs.0, Linear::add),
					std::array::from_fn(|_| unsafe {
						// SAFETY: Sizes of all input & output arrays are equal.
						a.next().unwrap_unchecked()
							.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
					}),
				)
			}
		}
		impl<T: Basis> Mul for SumPoly<T, { $($num +)+ 0 }> { // Squaring
			type Output = SumPoly<T, { 2 * ($($num +)+ 0) }>;
			fn mul(self, rhs: Self) -> Self::Output {
				const SIZE: usize = $($num +)+ 0;
				let SumPoly(a_value, a_terms) = self;
				let SumPoly(b_value, b_terms) = rhs;
				let terms = std::array::from_fn(|i|
					if i < SIZE {
						let mut term = T::each_map_inner(
							[
								a_terms[i].clone(),
								b_terms[i].clone(),
								b_value.clone(),
								a_value.clone(),
							],
							|[a, b, a_val, b_val]| {
								a.mul(b_val).add(b.mul(a_val))
							}
						);
						for j in 0..i {
							term = T::each_map_inner(
								[term, a_terms[j].clone(), b_terms[i-j-1].clone()],
								|[x, a, b]| x.add(a.mul(b))
							);
						}
						term
					} else {
						let mut term = T::zero();
						for j in (i - SIZE)..SIZE {
							term = T::each_map_inner(
								[term, a_terms[j].clone(), b_terms[i-j-1].clone()],
								|[x, a, b]| x.add(a.mul(b))
							);
						}
						term
					}
				);
				SumPoly(a_value.zip_map_inner(b_value, |a, b| a.mul(b)), terms)
			}
		}
		impl_deg_add!({ $($num +)+ 0 }, 1 $($num)+);
		impl_deg_order!(1 $($num)+);
	};
}
macro_rules! impl_deg_add {
	($a:tt, 1  1  $($num:tt)*) => { impl_deg_add!($a, 2  $($num)*); };
	($a:tt, 2  2  $($num:tt)*) => { impl_deg_add!($a, 4  $($num)*); };
	($a:tt, 4  4  $($num:tt)*) => { impl_deg_add!($a, 8  $($num)*); };
	// ($a:tt, 8  8  $($num:tt)*) => { impl_deg_add!($a, 16 $($num)*); };
	// ($a:tt, 16 16 $($num:tt)*) => { impl_deg_add!($a, 32 $($num)*); };
	// ($a:tt, 32 32 $($num:tt)*) => { impl_deg_add!($a, 64 $($num)*); };
	($a:tt, 8) => {/* break */};
	($a:tt, $($num:tt)+) => {
		impl<T: Basis> Add<Sum<T, $a>> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, $a>) -> Self::Output {
				let mut a = self.0.into_iter();
				let mut b = rhs.0.into_iter();
				Sum(std::array::from_fn(|i| unsafe {
					// SAFETY: `a` is the same size as the output array, and
					// `b` is the size of `$a` (bad naming).
					if i < $a {
						a.next().unwrap_unchecked()
							.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
					} else {
						a.next().unwrap_unchecked()
					}
				}))
			}
		}
		impl<T: Basis> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, $a> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.0.into_iter();
				let mut b = rhs.0.into_iter();
				Sum(std::array::from_fn(|i| unsafe {
					// SAFETY: `b` is the same size as the output array, and
					// `a` is the size of `$a`.
					if i < $a {
						a.next().unwrap_unchecked()
							.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
					} else {
						b.next().unwrap_unchecked()
					}
				}))
			}
		}
		impl<T: Basis> Add<SumPoly<T, $a>> for SumPoly<T, { $($num +)+ 0 }> {
			type Output = SumPoly<T, { $($num +)+ 0 }>;
			fn add(self, rhs: SumPoly<T, $a>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				SumPoly(
					self.0.zip_map_inner(rhs.0, Linear::add),
					std::array::from_fn(|i| unsafe {
						// SAFETY: `a` is the same size as the output array, and
						// `b` is the size of `$a` (bad naming).
						if i < $a {
							a.next().unwrap_unchecked()
								.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
						} else {
							a.next().unwrap_unchecked()
						}
					}),
				)
			}
		}
		impl<T: Basis> Add<SumPoly<T, { $($num +)+ 0 }>> for SumPoly<T, $a> {
			type Output = SumPoly<T, { $($num +)+ 0 }>;
			fn add(self, rhs: SumPoly<T, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				SumPoly::new(
					self.0.zip_map_inner(rhs.0, Linear::add),
					std::array::from_fn(|i| unsafe {
						// SAFETY: `b` is the same size as the output array, and
						// `a` is the size of `$a`.
						if i < $a {
							a.next().unwrap_unchecked()
								.zip_map_inner(b.next().unwrap_unchecked(), Linear::add)
						} else {
							b.next().unwrap_unchecked()
						}
					}),
				)
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);

impl<K, T, const N: usize> Sub<K> for Sum<T, N>
where
	K: Change + Neg<Output: Change<Basis: Basis<Inner = T::Inner>>>,
	T: Basis,
	Self: Add<<K as Neg>::Output>,
{
	type Output = <Self as Add<<K as Neg>::Output>>::Output;
	fn sub(self, rhs: K) -> Self::Output {
		self + -rhs
	}
}

impl<K, T, const N: usize> Sub<K> for SumPoly<T, N>
where
	K: Neg,
	Self: Add<K::Output>,
{
	type Output = <Self as Add<K::Output>>::Output;
	fn sub(self, rhs: K) -> Self::Output {
		self + -rhs
	}
}

impl Roots for SumPoly<f64, 0> {
	type Output = [f64; 0];
	fn roots(self) -> <Self as Roots>::Output {
		[]
	}
}

impl Roots for SumPoly<f64, 1> {
	type Output = [f64; 1];
	fn roots(self) -> <Self as Roots>::Output {
		let SumPoly(a, [b]) = self;
		[-a / b]
	}
}

impl Roots for SumPoly<f64, 2> {
	type Output = [f64; 2];
	fn roots(self) -> <Self as Roots>::Output {
		let SumPoly(a, [b, c]) = self;
		let mut n = -b / c;
		
		 // Precision Breakdown (Discriminant Quotient):
		let x = -a / b;
		if x.is_nan() || n.is_nan() || (x / n).abs() < 1e-8 {
			return if x.is_nan() {
				[n, n]
			} else {
				[x, n-x]
			}
		}
		
		 // General Quadratic:
		n /= 2.;
		let x = n.mul_add(n, -a/c).sqrt();
		[n-x, n+x]
	}
}

impl Roots for SumPoly<f64, 3> {
	type Output = [f64; 3];
	fn roots(self) -> <Self as Roots>::Output {
		let SumPoly(a, [b, c, d]) = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x   * c) / b).abs() < 1e-5 && // ??? Adjust as needed
			((x*x * d) / b).abs() < 1e-16
		) {
			let [y, z] = SumPoly(b, [c, d]).roots();
			return if x.is_nan() {
				[y, y, z]
			} else {
				[x, y-x, z]
			}
		}
		
		let mut n = -c / d;
		
		 // Weak Leading Term:
		if n.is_nan() || (
			(b / (n   * c)).abs() < 1e-9 && // ??? Adjust as needed
			(a / (n*n * c)).abs() < 1e-13
		) {
			let [x, y] = SumPoly(a, [b, c]).roots();
			return [x, y, n]
		}
		
		if c == 0. {
			 // Pseudo-linear:
			if b == 0. {
				return [(-a / d).cbrt(); 3]
			}
			
			 // Depressed Cubic:
			let p = -a / (2. * d);
			let q = -b / (3. * d);
			let r = p / q;
			let discriminant = if q < 0. {
				1.
			} else {
				let disc = f64::ln((r*r) / q.abs());
				if disc.abs() < 1e-15 {
					0.
				} else {
					disc
				}
			};
			return match discriminant.partial_cmp(&0.) {
				 // 3 Real Roots:
				Some(Ordering::Less) => {
					let sqrt_q = q.sqrt();
					let angle = f64::acos(r / sqrt_q);
					debug_assert!(!angle.is_nan(), "{:?}", self);
					use std::f64::consts::TAU;
					[
						2. * sqrt_q * f64::cos((angle         ) / 3.),
						2. * sqrt_q * f64::cos((angle +    TAU) / 3.),
						2. * sqrt_q * f64::cos((angle + 2.*TAU) / 3.),
					] 
				},
				
				 // 1 Real Root:
				Some(Ordering::Greater) => {
					let n = q.mul_add(-q*q, p*p).sqrt();
					debug_assert!(!n.is_nan());
					[(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN]
				},
				
				 // 2 Real Roots:
				_ => [-r, -r, 2.*r]
			}
		}
		
		 // General Cubic:
		n /= 3.;
		let depressed_cubic = SumPoly(
			n.mul_add(n.mul_add(c * (2./3.), b), a),
			[
				n.mul_add(c, b),
				0.,
				d,
			]
		);
		let [x, y, z] = depressed_cubic.roots();
		[x+n, y+n, z+n]
	}
}

impl Roots for SumPoly<f64, 4> {
	type Output = [f64; 4];
	fn roots(self) -> <Self as Roots>::Output {
		let SumPoly(a, [b, c, d, e]) = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x     * c) / b).abs() < 1e-4  && // ??? Adjust as needed
			((x*x   * d) / b).abs() < 1e-12 &&
			((x*x*x * e) / b).abs() < 1e-20
		) {
			let [y, z, w] = SumPoly(b, [c, d, e]).roots();
			return if x.is_nan() {
				[y, y, z, w]
			} else {
				[x, y-x, z, w]
			}
		}
		
		let mut n = -d / e;
		
		 // Weak Leading Term:
		if n.is_nan() || (
		    (c / (n     * d)).abs() < 1e-7  && // ??? Adjust as needed
			(b / (n*n   * d)).abs() < 1e-10 &&
			(a / (n*n*n * d)).abs() < 1e-15
		) {
			let [x, y, z] = SumPoly(a, [b, c, d]).roots();
			return [x, y, z, n]
		}
		
		if d == 0. {
			if b == 0. {
				if c == 0. {
					 // Pseudo-linear:
					let r = (-a / e).sqrt().sqrt();
					return [r, r, -r, -r]
				}
				
				 // Biquadratic:
				let [x, y] = SumPoly(a, [c, e]).roots();
				let (x, y) = (x.sqrt(), y.sqrt());
				return [-x, x, -y, y];
			}
			
			 // Depressed Quartic:
			let p = a / e;
			let q = b / (2. * e);
			let r = c / (2. * e);
			let resolvent_cubic = SumPoly(
				-(q * q) / 2.,
				[
					r.mul_add(r, -p),
					2. * r,
					1.,
				]
			);
			let mut m = resolvent_cubic.roots().into_iter()
				.find(|&r| r > 0. && r.is_finite())
				.unwrap_or_else(|| panic!(
					"expected a positive root from: {:?}",
					resolvent_cubic
				));
			let y = resolvent_cubic.eval(m)
				/ m.mul_add(m.mul_add(3., 4.*r), resolvent_cubic.1[0]);
			if m > y && y.is_finite() {
				m -= y; // Newton-Raphson step
			}
			let sqrt_2m = (2. * m).sqrt();
			let [x, y] = SumPoly( (q / sqrt_2m) + r + m, [-sqrt_2m, 1.]).roots();
			let [z, w] = SumPoly(-(q / sqrt_2m) + r + m, [ sqrt_2m, 1.]).roots();
			return [x, y, z, w]
		}
		
		 // General Quartic:
		n /= 4.;
		let depressed_quartic = SumPoly(
			n.mul_add(n.mul_add(n.mul_add(d * (3./4.), c), b), a),
			[
				n.mul_add(n.mul_add(d, c) * 2., b),
				n.mul_add(d * (3./2.), c),
				0.,
				e,
			]
		);
		let [x, y, z, w] = depressed_quartic.roots();
		[x+n, y+n, z+n, w+n]
	}
}

impl<A, B, const D: usize> Roots for SumPoly<Iso<A, B>, D>
where
	A: Basis,
	B: LinearIso<A>,
	SumPoly<A, D>: Roots,
{
	type Output = <SumPoly<A, D> as Roots>::Output;
	fn roots(self) -> <Self as Roots>::Output {
		let sum = SumPoly::<A, D>::new(
			self.0.into_inner(),
			self.1.map(|x| x.into_inner())
		);
		<SumPoly<A, D> as Roots>::roots(sum)
	}
}

#[cfg(feature = "glam")]
impl<const D: usize> Roots for SumPoly<glam::DVec2, D>
where
	SumPoly<f64, D>: Poly<Basis=f64> + Roots,
	<SumPoly<f64, D> as Roots>::Output: IntoIterator<Item=f64>,
{
	type Output = [f64; D];
	fn roots(self) -> <Self as Roots>::Output {
		let mut root_list = [f64::NAN; D];
		let mut root_count = 0;
		
		let mut x_poly = SumPoly::zero();
		let mut y_poly = SumPoly::zero();
		for (index, coeff) in self.into_iter().enumerate() {
			x_poly[index] = coeff.x;
			y_poly[index] = coeff.y;
		}
		
		let mut x_roots = x_poly.roots().into_iter().collect::<Vec<_>>();
		let mut y_roots = y_poly.roots().into_iter().collect::<Vec<_>>();
		x_roots.sort_unstable_by(f64::total_cmp);
		y_roots.sort_unstable_by(f64::total_cmp);
		
		let mut x_iter = x_roots.into_iter();
		let mut y_iter = y_roots.into_iter();
		let mut x = x_iter.next();
		let mut y = y_iter.next();
		while let (Some(a), Some(b)) = (x.as_ref(), y.as_ref()) {
			match a.total_cmp(b) {
				Ordering::Less    => x = x_iter.next(),
				Ordering::Greater => y = y_iter.next(),
				Ordering::Equal   => {
					root_list[root_count] = *a;
					root_count += 1;
					x = x_iter.next();
					y = y_iter.next();
					// !!! What if the next value(s) are equal to the current?
					// Should only one of these be advanced?
				}
			}
		}
		
		root_list
	}
}

#[cfg(test)]
mod tests {
	use crate::time::*;
	use crate::temporal::Temporal;
	use crate::pred::Prediction;
	use crate::*;
	use super::*;
	
	fn real_roots<K>(poly: Temporal<K>) -> impl Iterator<Item=f64>
	where
		K: Roots,
		<K as Roots>::Output: IntoIterator<Item=f64>,
	{
		poly.inner.roots().into_iter()
			.filter(|r| r.is_finite())
	}
	
	#[test]
	fn add() {
		let a = SumPoly(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = SumPoly(7.1, [5.9, 3.1]);
		assert_eq!(a + b, SumPoly(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
		assert_eq!(b + a, SumPoly(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
	}
	
	#[test]
	fn sub() {
		let a = SumPoly(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = SumPoly(7.1, [5.9, 3.1]);
		assert_eq!(a - b, SumPoly(
			-5.6,
			[-3.4000000000000004, 0.10000000000000009, 4.5, 5.7]
		));
		assert_eq!(b - a, SumPoly(
			5.6,
			[3.4000000000000004, -0.10000000000000009, -4.5, -5.7]
		));
	}
	
	#[test]
	fn scalar_mul() {
		let a = SumPoly(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = SumPoly(7.1, [5.9, 3.1]);
		assert_eq!(a.map(|x| x * 1.5), SumPoly(2.25, [3.75, 4.800000000000001, 6.75, 8.55]));
		assert_eq!(b.map(|x| x * 1.5), SumPoly(10.649999999999999, [8.850000000000001, 4.65]));
	}
	
	fn assert_roots<K: Poly>(p: K, expected_roots: &[f64])
	where
		K: Roots,
		<K as Roots>::Output: IntoIterator<Item=f64>,
	{
		let mut r = real_roots(Temporal::from(p)).into_iter()
			.collect::<Vec<f64>>();
		let mut expected_roots = expected_roots.into_iter()
			.copied()
			.collect::<Vec<f64>>();
		r.sort_unstable_by(f64::total_cmp);
		expected_roots.sort_unstable_by(f64::total_cmp);
		assert_eq!(
			r.len(), expected_roots.len(),
			"{:?} vs {:?}",
			r, expected_roots
		);
		for (x, y) in r.into_iter().zip(expected_roots.into_iter()) {
			assert!(
				(x - y).abs() < 1e-9 * 10_f64.powf(1. + y.abs().log10().ceil().max(0.)),
				"{:?} vs {:?}",
				x, y
			);
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(SumPoly(2., []), &[]);
		assert_roots(SumPoly(0., []), &[]);
	}
	
	#[test]
	fn linear() {
		assert_roots(SumPoly(20., [-4.]), &[5.]);
		assert_roots(SumPoly(20., [4.]), &[-5.]);
		assert_roots(SumPoly(-0., [4./3.]), &[0.]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(SumPoly(20., [4., 7.]), &[]);
		assert_roots(
			SumPoly(20., [4., -7.]),
			&[-10./7., 2.]
		);
		let r = real_roots(Temporal::from(SumPoly(-40./3., [-2./3., 17./100.])))
			.into_iter().collect::<Vec<_>>();
		assert_roots(
			SumPoly(40./3., [2./3., -17./100.]),
			r.as_slice()
		);
		assert_roots(
			SumPoly(0., [4./6., -17./100.]),
			&[0., 200./51.]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			SumPoly(6., [-2077.5, -17000./77., 6712./70.]),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		assert_roots(
			SumPoly(-3.6350710512584225e-33, [
				0.10240000000000019,
				-0.5455003017809628,
				1.
			]),
			&[3.54987407349455e-32]
		);
		assert_eq!(
			real_roots(Temporal::from(SumPoly(
				-236263115684.8131,
				[-9476965815.566229, -95034754.784949, 1.0]
			))).count(),
			3
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> SumPoly<f64, 3> {
			SumPoly(a, [
				b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s),
				c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + d/(2.*3.*s*s*s),
				d/(2.*3.*s*s*s),
			])
		}
		assert_roots(
			sum_poly(1000000000., 1., 1., 1., -1.),
			&[4591405718.8520711602007]
		);
		assert_roots(
			sum_poly(1000000000.*60., 1., 1., 1., -1.),
			&[275484343229.41354766068234]
		);
		assert_roots(
			sum_poly(1000000000.*60.*60., 1., 1., 1., -1.),
			&[16529060593863.10213769318]
		);
		assert_roots(
			sum_poly(1000000000.*60.*60., 1000., 300., 10., -4.),
			&[
				-55432951711728.79099027553,
				-13201315814560.382043758431976512,
				95634267526286.173034033963943,
			]
		);
	}
	
	#[test]
	fn quartic() {
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> SumPoly<f64, 4> {
			SumPoly(a, [
				b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s),
				c/(2.*s*s) + d/(2.*3.*s*s*s) + (2.*d)/(2.*3.*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s),
				d/(2.*3.*s*s*s) + e/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s),
				e/(2.*3.*4.*s*s*s*s),
			])
		}
		assert_roots(
			sum_poly(1000000000.*60.*60., 7.5, -6.4, -10., 7.8, 2.7),
			&[
				-51671989204288.265625,
				-5830778969957.388671875,
				2846556446843.169921875,
				13056211727396.474609375,
			]
		);
		let t = 1000000000.*60.*60.;
		assert_roots(
			SumPoly(
				7500.,
				[0., -1000./(t*t), 0., 27./(t*t*t*t)]
			),
			&[
				-2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
				-2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				 2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				 2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
			]
		);
		assert_roots(
			SumPoly(
				4.906800313619897e-6,
				[
					15145.164260286278,
					23999.999769993723,
					-236643.191,
					250000.0
				]
			),
			&[
				-0.19230902079054427,
				-3.23984644667874e-10,
				0.4732863823239847,
				0.6655954027905443
			]
		);
		assert_roots(
			SumPoly(-35.99999999999999, [8.046918796812349e-6, 2999.99988981303, -54772.25541522836, 250000.0]),
			&[-0.067702232, 0.177246743]
		);
		assert_eq!(Temporal::new(SumPoly(-181.99999999999994, [-7.289202428347473e-6, -500.0]), 1209618449*NANOSEC).eval(1209618450*NANOSEC), -181.99999999999994);
		assert_roots(
			SumPoly(9.094947017729282e-13, [-1.7967326421342023e-5, -8983.663173028655, 997.5710159206409, 250000.0]),
			&[-0.191570016, -1.11113e-8, 9.11132e-9, 0.187579734]
		);
	}
	
	#[cfg(feature = "glam")]
	#[test]
	fn vec() {
		use crate as chime;
		
		#[derive(PartialEq, Flux, ToMoment, ToMomentMut)]
		struct Pos {
			#[basis]
			value: Iso<glam::DVec2, glam::IVec2>,
			#[change(add_per(SEC))]
			spd: Spd,
		}
		
		#[derive(PartialEq, Flux, ToMoment, ToMomentMut)]
		struct Spd {
			#[basis]
			value: Iso<glam::DVec2, glam::IVec2>
		}
		
		let a_pos = Temporal::from(Pos {
			value: glam::IVec2::new(10, 10).into(),
			spd: Spd {
				value: glam::IVec2::new(6, 4).into()
			}
		});
		
		let b_pos = Temporal::from(Pos {
			value: glam::IVec2::new(14, 18).into(),
			spd: Spd {
				value: glam::IVec2::new(4, 0).into()
			}
		});
		
		// println!("{:?}", a_pos.poly() - b_pos.poly());
		
		assert_eq!(
			Vec::from_iter(a_pos.to_poly().when_eq(b_pos.to_poly())
				.into_ranges(Time::ZERO)
				.inclusive()),
			[(2*SEC - 83333333*NANOSEC, 2*SEC + 83333333*NANOSEC)]
		);
	}
	
	#[test]
	fn precise() {
		// https://www.desmos.com/calculator/1z97cqlopx
		
		 // Basic Check:
		let a = Temporal::new(SumPoly::new(-193.99999999999997, [4.481238217799146e-6, -500.]), SEC);
		let b = Temporal::new(Constant::from(-194.), SEC);
		assert_eq!(
			Vec::from_iter(a.when_eq(b)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(SEC-5*NANOSEC, SEC-3*NANOSEC),
				(SEC+12*NANOSEC, SEC+14*NANOSEC)
			]
		);
		
		 // Distance Check:
		let a = Temporal::new([
			SumPoly::new(0.0036784761334161292, [1.1687626970174242e-7, 0.]),
			SumPoly::new(-182.00000057575835, [-7.537214753878195e-7, -500.])
		], SEC);
		let b = Temporal::new([
			Constant::from(-3.8808053943969956e-5),
			Constant::from(-193.99999999999997)
		], SEC);
		let dis = Temporal::new(Constant::from(12.), SEC);
		assert_eq!(
			Vec::from_iter(a.when_dis_eq(b, dis)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(780910981*NANOSEC, 780910981*NANOSEC),
				(SEC-13*NANOSEC, SEC-10*NANOSEC),
				(SEC+2*NANOSEC, SEC+8*NANOSEC),
				(1219089016*NANOSEC, 1219089016*NANOSEC),
			]
		);
	}
}