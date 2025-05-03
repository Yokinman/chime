//! Summation over time.

use std::ops::{Add, Mul, Sub};
use symb_poly::{SymbolizeInner, Unsymbolize};
use crate::{change::*, Flux, Rate};
use crate::poly::{Poly, Roots};

mod _impl_roots {
	use crate::constant::Constant2;
	use crate::poly::Poly;
	use std::cmp::Ordering;
	use num_traits::Pow;
	use super::Roots;
	use symb_poly::{Func, Invar, Num, Power, Prod, Sum, Var};
	use typenum::{P2, P3, P4};

	impl Roots for Invar<Constant2<f64>> {
		type Output = [f64; 0];
		fn roots(self) -> <Self as Roots>::Output {
			[0.; 0]
		}
	}

	impl Roots for Func<Prod, (Invar<Constant2<f64>>, Var)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), _) = self.into_args();
			[0. / a; 1]
		}
	}
	
	impl Roots for Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), _) = self.into_args();
			[0. / a; 2]
		}
	}
	
	impl Roots for Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)> {
		type Output = [f64; 3];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), _) = self.into_args();
			[0. / a; 3]
		}
	}
	
	impl Roots for Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)> {
		type Output = [f64; 4];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), _) = self.into_args();
			[0. / a; 4]
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Prod, (Invar<Constant2<f64>>, Var)>
	)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_prod) = self.into_args();
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			[-a / b; 1]
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_prod) = self.into_args();
			
			if a == 0. {
				return b_prod.roots()
			}
			
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 2];
			roots[0] = (-a / b).powf((2 as f64).recip());
			roots[1] = -roots[0];
			roots
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_prod) = self.into_args();
			
			if a == 0. {
				return b_prod.roots()
			}
			
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 3];
			roots[0] = (-a / b).powf((3 as f64).recip());
			roots
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_prod) = self.into_args();
			
			if a == 0. {
				return b_prod.roots()
			}
			
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 4];
			roots[0] = (-a / b).powf((4 as f64).recip());
			roots[1] = -roots[0];
			roots
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>
		)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, c_prod) = b_sum.into_args();
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let (Invar(Constant2(c)), _) = c_prod.into_args();
			
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
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, d_prod) = b_sum.into_args();
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let (Invar(Constant2(d)), d_var) = d_prod.into_args();
			
			 // Pseudo-linear:
			if b == 0. {
				return (Invar(Constant2(a)) + Invar(Constant2(d))*d_var).roots()
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
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			(self + (Invar(Constant2(0.)) * <Var>::default())).roots()
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (c_prod, d_prod) = c_sum.into_args();
			let (Invar(Constant2(b)), b_var) = b_prod.into_args();
			let (Invar(Constant2(c)), c_var) = c_prod.into_args();
			let (Invar(Constant2(d)), d_var) = d_prod.into_args();
			
			 // Weak Constant:
			let x = -a / b;
			if x.is_nan() || (
				((x   * c) / b).abs() < 1e-5 && // ??? Adjust as needed
				((x*x * d) / b).abs() < 1e-16
			) {
				let [y, z] = (Invar(Constant2(b)) + Invar(Constant2(c))*b_var + Invar(Constant2(d))*c_var).roots();
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
				let [x, y] = (Invar(Constant2(a)) + Invar(Constant2(b))*b_var + Invar(Constant2(c))*c_var).roots();
				return [x, y, n]
			}
			
			 // Depressed Cubic:
			if c == 0. {
				return (Invar(Constant2(a)) + Invar(Constant2(b))*b_var + Invar(Constant2(d))*d_var).roots()
			}
			
			 // General Cubic:
			n /= 3.;
			let depressed_cubic = Invar(Constant2(n.mul_add(n.mul_add(c * (2./3.), b), a)))
				+ Invar(Constant2(n.mul_add(c, b)))*b_var
				+ Invar(Constant2(d))*d_var;
			depressed_cubic.roots().map(|x| x + n)
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self
				+ Invar(Constant2(0.))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant2(0.))*<Var>::default().pow(Invar(Num::<P3>::default())))
				.roots()
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), c_sum) = self.into_args();
			let (c_prod, e_prod) = c_sum.into_args();
			let (Invar(Constant2(c)), _) = c_prod.into_args();
			
			 // Pseudo-linear:
			if c == 0. {
				return (Invar(Constant2(a)) + e_prod).roots()
			}
			
			let (Invar(Constant2(e)), _) = e_prod.into_args();
			
			 // Biquadratic:
			let [x, y] = (Invar(Constant2(a))
				+ Invar(Constant2(c))*<Var>::default()
				+ Invar(Constant2(e))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			let (x, y) = (x.sqrt(), y.sqrt());
			return [-x, x, -y, y];
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self
				+ Invar(Constant2(0.))
				+ Invar(Constant2(0.))*<Var>::default().pow(Invar(Num::<P2>::default())))
				.roots()
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			
			 // Biquadratic:
			if b == 0. {
				return (Invar(Constant2(a)) + c_sum).roots()
			}
			
			let (c_prod, e_prod) = c_sum.into_args();
			let (Invar(Constant2(c)), _) = c_prod.into_args();
			let (Invar(Constant2(e)), _) = e_prod.into_args();
			
			 // Depressed Quartic:
			let p = a / e;
			let q = b / (2. * e);
			let r = c / (2. * e);
			let resolvent_b = r.mul_add(r, -p);
			let resolvent_cubic = Invar(Constant2(-(q * q) / 2.))
				+ Invar(Constant2(resolvent_b))*<Var>::default()
				+ Invar(Constant2(2. * r))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant2(1.))*<Var>::default().pow(Invar(Num::<P3>::default()));
			let mut m = resolvent_cubic.roots().into_iter()
				.find(|&r| r > 0. && r.is_finite())
				.unwrap_or_else(|| panic!(
					"expected a positive root from: {:?}",
					resolvent_cubic
				));
			let y = resolvent_cubic.eval(m)
				/ m.mul_add(m.mul_add(3., 4.*r), resolvent_b);
			if m > y && y.is_finite() {
				m -= y; // Newton-Raphson step
			}
			let sqrt_2m = (2. * m).sqrt();
			let [x, y] = (Invar(Constant2((q / sqrt_2m) + r + m))
				+ Invar(Constant2(-sqrt_2m))*<Var>::default()
				+ Invar(Constant2(1.))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			let [z, w] = (Invar(Constant2(-(q / sqrt_2m) + r + m))
				+ Invar(Constant2(sqrt_2m))*<Var>::default()
				+ Invar(Constant2(1.))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			[x, y, z, w]
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self + Invar(Constant2(0.))*<Var>::default().pow(Invar(Num::<P2>::default()))).roots()
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self + Invar(Constant2(0.))*<Var>::default()).roots()
		}
	}
	
	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Sum, (
					Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
					Func<Prod, (Invar<Constant2<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
				)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (c_prod, d_sum) = c_sum.into_args();
			let (d_prod, e_prod) = d_sum.into_args();
			let (Invar(Constant2(b)), b_var) = b_prod.into_args();
			let (Invar(Constant2(c)), c_var) = c_prod.into_args();
			let (Invar(Constant2(d)), d_var) = d_prod.into_args();
			let (Invar(Constant2(e)), e_var) = e_prod.into_args();
			
			 // Weak Constant:
			let x = -a / b;
			if x.is_nan() || (
				((x     * c) / b).abs() < 1e-4  && // ??? Adjust as needed
				((x*x   * d) / b).abs() < 1e-12 &&
				((x*x*x * e) / b).abs() < 1e-20
			) {
				let [y, z, w] = (Invar(Constant2(b)) + Invar(Constant2(c))*b_var + Invar(Constant2(d))*c_var + Invar(Constant2(e))*d_var).roots();
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
				let [x, y, z] = (Invar(Constant2(a)) + Invar(Constant2(b))*b_var + Invar(Constant2(c))*c_var + Invar(Constant2(d))*d_var).roots();
				return [x, y, z, n]
			}
			
			 // Depressed Quartic:
			if d == 0. {
				return (Invar(Constant2(a)) + Invar(Constant2(b))*b_var + Invar(Constant2(c))*c_var + Invar(Constant2(e))*e_var).roots()
			}
			
			 // General Quartic:
			n /= 4.;
			let depressed_quartic = Invar(Constant2(n.mul_add(n.mul_add(n.mul_add(d * (3./4.), c), b), a)))
				+ Invar(Constant2(n.mul_add(n.mul_add(d, c) * 2., b)))*<Var>::default()
				+ Invar(Constant2(n.mul_add(d * (3./2.), c)))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant2(e))*<Var>::default().pow(Invar(Num::<P4>::default()));
			let [x, y, z, w] = depressed_quartic.roots();
			[x+n, y+n, z+n, w+n]
		}
	}
	
	impl<A, B, const M: usize, const N: usize> Roots for Func<Sum, (Func<Prod, (Invar<Constant2<f64>>, A)>, B)>
	where
		Self: Poly + std::ops::Div<Var, Output: Roots<Output = [f64; M]>>,
		B: Roots<Output = [f64; N]>,
	{
		type Output = [f64; N];
		fn roots(self) -> <Self as Roots>::Output {
			let mut roots = [0.; N];
			roots[1..=M].copy_from_slice(&(self / Var::default()).roots());
			roots
		}
	}

	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Prod, (
			Invar<Constant2<f64>>,
			Func<Power, (Invar<Constant2<f64>>, Var)>
		)>
	)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_prod) = self.into_args();
			let (Invar(Constant2(b)), c_power) = b_prod.into_args();
			let (Invar(Constant2(c)), _) = c_power.into_args();
			if c == 1. {
				return [1. - a/b]
			}
			[(-a/b).log(c)]
		}
	}

	impl Roots for Func<Sum, (
		Invar<Constant2<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant2<f64>>, Var)>,
			Func<Prod, (
				Invar<Constant2<f64>>,
				Func<Power, (Invar<Constant2<f64>>, Var)>
			)>
		)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots>::Output {
			let (Invar(Constant2(a)), b_sum) = self.into_args();
			let (b_prod, c_prod) = b_sum.into_args();
			let (Invar(Constant2(b)), _) = b_prod.into_args();
			let (Invar(Constant2(c)), d_power) = c_prod.into_args();
			let (Invar(Constant2(d)), _) = d_power.into_args();

			if d == 1. {
				todo!()
			}
			// Convert to `a + b*x + c^x` and solve using product log.
			let x = a / c;
			let y = b / c;
			let d_ln = d.ln();
			let z = d_ln / (y * d.powf(x / y));
			[
				-(lambert_w::lambert_wm1(z) / d_ln) - (x / y),
				-(lambert_w::lambert_w0(z)  / d_ln) - (x / y),
			]
		}
	}
}

/// The pattern of repeated addition, `a = a + b`.
/// 
/// In essence, the second parameter is added to the first parameter, e.g.
/// 
/// - `Sum<Prod<A, B>, C>` ~ `a = a*b + c`
/// - `Sum<A, Prod<B, C>>` ~ `a = a + b*c`
#[derive(Copy, Clone, Debug, Default)]
pub struct Sum2<A: Change, B> {
	pub(crate) lhs: A,
	pub(crate) rhs: B,
	pub(crate) rhs_basis: A::Basis,
}

impl<A, B> Change for Sum2<A, B>
where
	A: AddChange<B::Poly>,
	B: Change<Basis = <A as Change>::Basis>,
{
	type Basis = A::Basis;
	type Poly = A::Output;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		self.lhs.add_change(basis, self.rhs.into_poly(self.rhs_basis))
	}
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		Self {
			lhs: self.lhs.scale(scalar),
			rhs: self.rhs.scale(scalar),
			rhs_basis: self.rhs_basis.map_inner(|n| n.mul(scalar)),
		}
	}
}

impl<A, B, C> Add<Rate<C>> for Sum2<A, B>
where
	A: Change,
	Self: std::ops::Shr<Sum2<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Sum2<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = <Self as std::ops::Shr<Sum2<Nil<C::Basis>, C::Change>>>::Output;
	fn add(self, rhs: Rate<C>) -> Self::Output {
		self >> Sum2 {
			lhs: Nil::default(),
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip()))
	}
}

impl<A, B, C> Sub<Rate<C>> for Sum2<A, B>
where
	A: Change,
	Self: std::ops::Shr<Sum2<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Sum2<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = <Self as std::ops::Shr<Sum2<Nil<C::Basis>, C::Change>>>::Output;
	fn sub(self, rhs: Rate<C>) -> Self::Output {
		self >> Sum2 {
			lhs: Nil::default(),
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(-rhs.unit.as_secs_f64().recip()))
	}
}

impl<A, B, C> Mul<Rate<C>> for Sum2<A, B>
where
	A: Change,
	Self: std::ops::Shr<Prod<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Prod<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = <Self as std::ops::Shr<Prod<Nil<C::Basis>, C::Change>>>::Output;
	fn mul(self, rhs: Rate<C>) -> Self::Output {
		self >> Prod {
			lhs: Nil::default(),
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip()))
	}
}

impl<A, B, C> Div<Rate<C>> for Sum2<A, B>
where
	A: Change,
	Self: std::ops::Shr<Prod<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Prod<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = <Self as std::ops::Shr<Prod<Nil<C::Basis>, C::Change>>>::Output;
	fn div(self, rhs: Rate<C>) -> Self::Output {
		self >> Prod {
			lhs: Nil::default(),
			rhs: rhs.amount.change(),
			rhs_basis: rhs.amount.basis(),
		}.scale(Linear::from_f64(rhs.unit.as_secs_f64()/*.recip()*/))
	}
}

impl<A, B, T> std::ops::Shr<Nil<T>> for Sum2<A, B>
where
	A: Change
{
	type Output = Self;
	fn shr(self, _rhs: Nil<T>) -> Self::Output {
		self
	}
}

impl<A, B, X, Y> std::ops::Shr<Sum2<X, Y>> for Sum2<A, B>
where
	A: Change,
	X: Change<Basis = A::Basis>,
	A: std::ops::Shr<X, Output: Change<Basis = A::Basis>>,
	B: std::ops::Shr<Y>,
{
	type Output = Sum2<A::Output, B::Output>;
	fn shr(self, rhs: Sum2<X, Y>) -> Self::Output {
		Sum2 {
			lhs: self.lhs >> rhs.lhs,
			rhs: self.rhs >> rhs.rhs,
			rhs_basis: self.rhs_basis.zip_map_inner(rhs.rhs_basis, Linear::add),
		}
	}
}

impl<A, B, X, Y> std::ops::Shr<Prod<X, Y>> for Sum2<A, B>
where
	A: Change,
	X: Change<Basis = A::Basis>,
	Self: std::ops::Shr<X, Output: Change<Basis = A::Basis>>,
{
	type Output = Prod<<Self as std::ops::Shr<X>>::Output, Y>;
	fn shr(self, rhs: Prod<X, Y>) -> Self::Output {
		Prod {
			lhs: self >> rhs.lhs,
			rhs: rhs.rhs,
			rhs_basis: rhs.rhs_basis,
		}
	}
}

/// The pattern of repeated multiplication, `a = a * b`.
/// 
/// In essence, the second parameter is multiplied into the first parameter, e.g.
/// 
/// - `Prod<Sum<A, B>, C>` ~ `a = (a+b) * c`
/// - `Prod<A, Sum<B, C>>` ~ `a = a * (b+c)`
#[derive(Copy, Clone, Debug, Default)]
pub struct Prod<A: Change, B> {
	pub(crate) lhs: A,
	pub(crate) rhs: B,
	pub(crate) rhs_basis: A::Basis,
}

impl<A, B> Change for Prod<A, B>
where
	A: MulChange<B::Poly>,
	B: Change<Basis = <A as Change>::Basis>,
{
	type Basis = A::Basis;
	type Poly = A::Output;
	fn into_poly(self, basis: Self::Basis) -> Self::Poly {
		self.lhs.mul_change(basis, self.rhs.into_poly(self.rhs_basis))
	}
	fn scale(self, scalar: <Self::Basis as Basis>::Inner) -> Self {
		Self {
			lhs: self.lhs.scale(scalar),
			rhs: self.rhs.scale(scalar),
			rhs_basis: self.rhs_basis.map_inner(|n| n.mul(scalar)),
		}
	}
}

impl<A, B, C> Add<Rate<C>> for Prod<A, B>
where
	A: Change,
	B: std::ops::Shr<Sum2<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Sum2<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = Prod<A, B::Output>;
	fn add(self, rhs: Rate<C>) -> Self::Output {
		Prod {
			lhs: self.lhs,
			rhs: self.rhs >> Sum2 {
				lhs: Nil::default(),
				rhs: rhs.amount.change(),
				rhs_basis: rhs.amount.basis(),
			}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip())),
			rhs_basis: self.rhs_basis,
		}
	}
}

impl<A, B, C> Mul<Rate<C>> for Prod<A, B>
where
	A: Change,
	B: std::ops::Shr<Prod<Nil<C::Basis>, C::Change>>,
	C: Flux<Basis = A::Basis>,
	Prod<Nil<C::Basis>, C::Change>: Change<Basis = A::Basis>,
{
	type Output = Prod<A, B::Output>;
	fn mul(self, rhs: Rate<C>) -> Self::Output {
		Prod {
			lhs: self.lhs,
			rhs: self.rhs >> Prod {
				lhs: Nil::default(),
				rhs: rhs.amount.change(),
				rhs_basis: rhs.amount.basis(),
			}.scale(Linear::from_f64(rhs.unit.as_secs_f64().recip())),
			rhs_basis: self.rhs_basis,
		}
	}
}

impl<A, B, T> std::ops::Shr<Nil<T>> for Prod<A, B>
where
	A: Change
{
	type Output = Self;
	fn shr(self, _rhs: Nil<T>) -> Self::Output {
		self
	}
}

impl<A, B, X, Y> std::ops::Shr<Prod<X, Y>> for Prod<A, B>
where
	A: Change,
	X: Change<Basis = A::Basis>,
	A: std::ops::Shr<X, Output: Change<Basis = A::Basis>>,
	B: std::ops::Shr<Y>,
{
	type Output = Prod<A::Output, B::Output>;
	fn shr(self, rhs: Prod<X, Y>) -> Self::Output {
		Prod {
			lhs: self.lhs >> rhs.lhs,
			rhs: self.rhs >> rhs.rhs,
			rhs_basis: self.rhs_basis.zip_map_inner(rhs.rhs_basis, Linear::mul),
		}
	}
}

impl<A, B, X, Y> std::ops::Shr<Sum2<X, Y>> for Prod<A, B>
where
	A: Change,
	X: Change<Basis = A::Basis>,
	Self: std::ops::Shr<X, Output: Change<Basis = A::Basis>>,
{
	type Output = Sum2<<Self as std::ops::Shr<X>>::Output, Y>;
	fn shr(self, rhs: Sum2<X, Y>) -> Self::Output {
		Sum2 {
			lhs: self >> rhs.lhs,
			rhs: rhs.rhs,
			rhs_basis: rhs.rhs_basis,
		}
	}
}

/// Helper trait for implementing [`Change`] for [`Sum2`].
pub trait AddChange<C>: Change {
	type Output: Poly<Basis = Self::Basis>;
	fn add_change(self, basis: Self::Basis, change: C) -> Self::Output;
}

// Nil<T> + Integ!(C)
impl<C, O, P, T> AddChange<C> for Nil<T>
where
	Self: Change<Basis=T, Poly: Add<O, Output=P>>,
	C: Integ<Output=O>,
	P: Poly<Basis=T>,
	T: Basis,
{
	type Output = P;
	fn add_change(self, basis: Self::Basis, change: C) -> Self::Output {
		self.into_poly(basis) + change.integ()
	}
}

// (A + Integ!(B)) + Integ!(C)
impl<A, B, C, O, P, T> AddChange<C> for Sum2<A, B>
where
	Self: Change<Basis=T, Poly: Add<P, Output=O>>,
	A: Change<Basis=T>,
	B: Change<Basis=T>,
	C: Integ<Output=P>,
	O: Poly<Basis=T>,
	T: Basis,
{
	type Output = O;
	fn add_change(self, basis: Self::Basis, change: C) -> Self::Output {
		self.into_poly(basis) + change.integ()
	}
}

// (A * ProdInteg!(B)) + Integ!(C) -> (A + Integ!(C / ProdInteg!(B))) * ProdInteg!(B)
impl<A, B, C, O, P, Q, T> AddChange<C> for Prod<A, B>
where
	Self: Change<Basis=T>,
	A: Change<Basis=T, Poly: SymbolizeInner<<B::Poly as SymbolizeInner<C::NextSymbol>>::NextSymbol, Output: Add<Q, Output: Mul<P, Output=O>>>>,
	B: Change<Basis=T, Poly: SymbolizeInner<C::NextSymbol, Output: Div<P, Output: Integ<Output=Q>>>>,
	C: SymbolizeInner<typenum::U0, Output: ProdInteg<Output=P>>,
	O: Unsymbolize<Output: Poly<Basis=T>>,
	P: Clone,
	T: Basis,
{
	type Output = O::Output;
	fn add_change(self, basis: Self::Basis, change: C) -> Self::Output {
		let (c_poly, c_symbol) = change.symbolize_inner(Default::default());
		let p = c_poly.prod_integ();
		let (b_poly, b_symbol) = self.rhs.into_poly(self.rhs_basis).symbolize_inner(c_symbol);
		let q = (b_poly / p.clone()).integ();
		let (a_poly, _) = self.lhs.into_poly(basis).symbolize_inner(b_symbol);
		((a_poly + q) * p).unsymbolize()
	}
}

/// Helper trait for implementing [`Change`] for [`Prod`].
pub trait MulChange<C>: Change {
	type Output: Poly<Basis = Self::Basis>;
	fn mul_change(self, basis: Self::Basis, change: C) -> Self::Output;
}

// Nil<T> * ProdInteg!(C)
impl<C, O, P, T> MulChange<C> for Nil<T>
where
	Self: Change<Basis=T, Poly: Mul<P, Output=O>>,
	C: ProdInteg<Output=P>,
	O: Poly<Basis=T>,
	T: Basis,
{
	type Output = O;
	fn mul_change(self, basis: Self::Basis, change: C) -> Self::Output {
		self.into_poly(basis) * change.prod_integ()
	}
}

// (A + Integ!(B)) * ProdInteg!(C) -> (A + Integ!(B * C/ProdInteg!(C))) * ProdInteg!(C)
impl<A, B, C, O, P, Q, T> MulChange<C> for Sum2<A, B>
where
	Self: Change<Basis=T>,
	A: Change<Basis=T, Poly: SymbolizeInner<<B::Poly as SymbolizeInner<C::NextSymbol>>::NextSymbol, Output: Add<Q, Output: Mul<P, Output=O>>>>,
	B: Change<Basis=T, Poly: SymbolizeInner<C::NextSymbol, Output: Mul<C::Output, Output: Div<P, Output: Integ<Output=Q>>>>>,
	C: SymbolizeInner<typenum::U0, Output: ProdInteg<Output=P> + Clone>,
	O: Unsymbolize<Output: Poly<Basis=T>>,
	P: Clone,
	T: Basis,
{
	type Output = O::Output;
	fn mul_change(self, basis: Self::Basis, change: C) -> Self::Output {
		let (c_poly, c_symbol) = change.symbolize_inner(Default::default());
		let p = c_poly.clone().prod_integ();
		let (b_poly, b_symbol) = self.rhs.into_poly(self.rhs_basis).symbolize_inner(c_symbol);
		let q = (b_poly * c_poly / p.clone()).integ();
		let (a_poly, _) = self.lhs.into_poly(basis).symbolize_inner(b_symbol);
		((a_poly + q) * p)
			.unsymbolize()
	}
}

// (A * ProdInteg!(B)) * ProdInteg!(C)
impl<A, B, C, O, P, T> MulChange<C> for Prod<A, B>
where
	Self: Change<Basis=T, Poly: Mul<P, Output=O>>,
	A: Change<Basis=T>,
	B: Change<Basis=T>,
	C: ProdInteg<Output=P>,
	O: Poly<Basis=T>,
	T: Basis,
{
	type Output = O;
	fn mul_change(self, basis: Self::Basis, change: C) -> Self::Output {
		self.into_poly(basis) * change.prod_integ()
	}
}

/// Integration, aka continuous summation.
pub trait Integ {
	type Output;
	fn integ(self) -> Self::Output;
}

impl<T> Integ for symb_poly::Invar<T>
where
	Self: symb_poly::IntegExpr<symb_poly::Variable>
{
	type Output = <Self as symb_poly::IntegExpr<symb_poly::Variable>>::Output;
	fn integ(self) -> Self::Output {
		symb_poly::IntegExpr::integ_expr(self, Default::default())
	}
}

impl<T> Integ for symb_poly::Var<T>
where
	Self: symb_poly::IntegExpr<symb_poly::Variable>
{
	type Output = <Self as symb_poly::IntegExpr<symb_poly::Variable>>::Output;
	fn integ(self) -> Self::Output {
		symb_poly::IntegExpr::integ_expr(self, Default::default())
	}
}

impl<F, A> Integ for symb_poly::Func<F, A>
where
	Self: symb_poly::IntegExpr<symb_poly::Variable>
{
	type Output = <Self as symb_poly::IntegExpr<symb_poly::Variable>>::Output;
	fn integ(self) -> Self::Output {
		symb_poly::IntegExpr::integ_expr(self, Default::default())
	}
}

/// Product integration, aka continuous product.
/// 
/// Blanket implemented for all types that implement [`ApplyLn`] -> [`Integ`] ->
/// [`ApplyExp`].
pub trait ProdInteg {
	type Output;
	fn prod_integ(self) -> Self::Output;
}

impl<T> ProdInteg for symb_poly::Invar<T>
where
	Self: symb_poly::Symbolize<Output: symb_poly::Ln<Output: symb_poly::IntegExpr<symb_poly::Variable, Output: symb_poly::Exp<Output: symb_poly::Unsymbolize>>>>
{
	type Output = <<<<<Self as symb_poly::Symbolize>::Output as symb_poly::Ln>::Output as symb_poly::IntegExpr<symb_poly::Variable>>::Output as symb_poly::Exp>::Output as symb_poly::Unsymbolize>::Output;
	fn prod_integ(self) -> Self::Output {
		symb_poly::Unsymbolize::unsymbolize(symb_poly::Exp::exp(symb_poly::IntegExpr::integ_expr(symb_poly::Ln::ln(symb_poly::Symbolize::symbolize(self)), Default::default())))
	}
}

impl<T> ProdInteg for symb_poly::Var<T>
where
	Self: symb_poly::Symbolize<Output: symb_poly::Ln<Output: symb_poly::IntegExpr<symb_poly::Variable, Output: symb_poly::Exp<Output: symb_poly::Unsymbolize>>>>
{
	type Output = <<<<<Self as symb_poly::Symbolize>::Output as symb_poly::Ln>::Output as symb_poly::IntegExpr<symb_poly::Variable>>::Output as symb_poly::Exp>::Output as symb_poly::Unsymbolize>::Output;
	fn prod_integ(self) -> Self::Output {
		symb_poly::Unsymbolize::unsymbolize(symb_poly::Exp::exp(symb_poly::IntegExpr::integ_expr(symb_poly::Ln::ln(symb_poly::Symbolize::symbolize(self)), Default::default())))
	}
}

impl<F, A> ProdInteg for symb_poly::Func<F, A>
where
	Self: symb_poly::Symbolize<Output: symb_poly::Ln<Output: symb_poly::IntegExpr<symb_poly::Variable, Output: symb_poly::Exp<Output: symb_poly::Unsymbolize>>>>
{
	type Output = <<<<<Self as symb_poly::Symbolize>::Output as symb_poly::Ln>::Output as symb_poly::IntegExpr<symb_poly::Variable>>::Output as symb_poly::Exp>::Output as symb_poly::Unsymbolize>::Output;
	fn prod_integ(self) -> Self::Output {
		symb_poly::Unsymbolize::unsymbolize(symb_poly::Exp::exp(symb_poly::IntegExpr::integ_expr(symb_poly::Ln::ln(symb_poly::Symbolize::symbolize(self)), Default::default())))
	}
}

// #[cfg(feature = "glam")]
// impl<const D: usize> Roots for SumPoly<glam::DVec2, D>
// where
// 	SumPoly<f64, D>: Poly<Basis=f64> + Roots,
// 	<SumPoly<f64, D> as Roots>::Output: IntoIterator<Item=f64>,
// {
// 	type Output = [f64; D];
// 	fn roots(self) -> <Self as Roots>::Output {
// 		let mut root_list = [f64::NAN; D];
// 		let mut root_count = 0;
//
// 		let mut x_poly = SumPoly::zero();
// 		let mut y_poly = SumPoly::zero();
// 		for (index, coeff) in self.into_iter().enumerate() {
// 			x_poly[index] = coeff.x;
// 			y_poly[index] = coeff.y;
// 		}
//
// 		let mut x_roots = x_poly.roots().into_iter().collect::<Vec<_>>();
// 		let mut y_roots = y_poly.roots().into_iter().collect::<Vec<_>>();
// 		x_roots.sort_unstable_by(f64::total_cmp);
// 		y_roots.sort_unstable_by(f64::total_cmp);
//
// 		let mut x_iter = x_roots.into_iter();
// 		let mut y_iter = y_roots.into_iter();
// 		let mut x = x_iter.next();
// 		let mut y = y_iter.next();
// 		while let (Some(a), Some(b)) = (x.as_ref(), y.as_ref()) {
// 			match a.total_cmp(b) {
// 				Ordering::Less    => x = x_iter.next(),
// 				Ordering::Greater => y = y_iter.next(),
// 				Ordering::Equal   => {
// 					root_list[root_count] = *a;
// 					root_count += 1;
// 					x = x_iter.next();
// 					y = y_iter.next();
// 					// !!! What if the next value(s) are equal to the current?
// 					// Should only one of these be advanced?
// 				}
// 			}
// 		}
//
// 		root_list
// 	}
// }

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
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(a + b, poly!(|x| 8.6 + 8.4*x + 6.300000000000001*x^2 + 4.5*x^3 + 5.7*x^4));
		assert_eq!(b + a, poly!(|x| 8.6 + 8.4*x + 6.300000000000001*x^2 + 4.5*x^3 + 5.7*x^4));
	}
	
	#[test]
	fn sub() {
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(
			a - b,
			poly!(|x| -5.6 - 3.4000000000000004*x + 0.10000000000000009*x^2 + 4.5*x^3 + 5.7*x^4)
		);
		assert_eq!(
			b - a,
			poly!(|x| 5.6 + 3.4000000000000004*x - 0.10000000000000009*x^2 - 4.5*x^3 - 5.7*x^4)
		);
	}
	
	#[test]
	fn scalar_mul() {
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(poly!(|_| 1.5 * a), poly!(|x| 2.25 + 3.75*x + 4.800000000000001*x^2 + 6.75*x^3 + 8.55*x^4));
		assert_eq!(poly!(|_| 1.5 * b), poly!(|x| 10.649999999999999 + 8.850000000000001*x + 4.65*x^2));
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
	
	// #[test]
	// fn constant() {
	// 	assert_roots(poly!(|_| 2.), &[]);
	// 	assert_roots(poly!(|_| 0.), &[]);
	// }
	
	#[test]
	fn linear() {
		assert_roots(poly!(|x| 20. + -4.*x), &[5.]);
		assert_roots(poly!(|x| 20. + 4.*x), &[-5.]);
		assert_roots(poly!(|x| 4./3.*x), &[0.]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(poly!(|x| 20. + 4.*x + 7.*x^2), &[]);
		assert_roots(
			poly!(|x| 20. + 4.*x - 7.*x^2),
			&[-10./7., 2.]
		);
		let r = real_roots(Temporal::from(poly!(|x| -40./3. - 2./3.*x + 17./100.*x^2)))
			.into_iter().collect::<Vec<_>>();
		assert_roots(
			poly!(|x| 40./3. + 2./3.*x - 17./100.*x^2),
			r.as_slice()
		);
		assert_roots(
			poly!(|x| 4./6.*x - 17./100.*x^2),
			&[0., 200./51.]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			poly!(|x| 6. - 2077.5*x - 17000./77.*x^2 + 6712./70.*x^3),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		assert_roots(
			poly!(|x| -3.6350710512584225e-33
				+ 0.10240000000000019*x
				+ -0.5455003017809628*x^2
				+ 1.*x^3),
			&[3.54987407349455e-32]
		);
		assert_eq!(
			real_roots(Temporal::from(poly!(|x| -236263115684.8131
				- 9476965815.566229*x
				- 95034754.784949*x^2
				+ 1.0*x^3
			))).count(),
			3
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> impl Roots<Output: IntoIterator<Item=f64>> {
			poly!(|x| a
				+ (b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s))*x
				+ (c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + d/(2.*3.*s*s*s))*x^2
				+ (d/(2.*3.*s*s*s))*x^3)
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
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> impl Roots<Output: IntoIterator<Item=f64>> {
			poly!(|x| a
				+ (b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s))*x
				+ (c/(2.*s*s) + d/(2.*3.*s*s*s) + (2.*d)/(2.*3.*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s))*x^2
				+ (d/(2.*3.*s*s*s) + e/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s))*x^3
				+ (e/(2.*3.*4.*s*s*s*s))*x^4)
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
			poly!(|x| 7500. - 1000./(t*t)*x^2 + 27./(t*t*t*t)*x^4),
			&[
				-2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
				-2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				 2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				 2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
			]
		);
		assert_roots(
			poly!(|x| 4.906800313619897e-6
				+ 15145.164260286278*x
				+ 23999.999769993723*x^2
				- 236643.191*x^3
				+ 250000.0*x^4
			),
			&[
				-0.19230902079054427,
				-3.23984644667874e-10,
				0.4732863823239847,
				0.6655954027905443
			]
		);
		assert_roots(
			poly!(|x| -35.99999999999999 + 8.046918796812349e-6*x + 2999.99988981303*x^2 - 54772.25541522836*x^3 + 250000.0*x^4),
			&[-0.067702232, 0.177246743]
		);
		assert_eq!(Temporal::new(poly!(|x| -181.99999999999994 - 7.289202428347473e-6*x - 500.*x^2), 1209618449*NANOSEC).eval(1209618450*NANOSEC), -181.99999999999994);
		assert_roots(
			poly!(|x| 9.094947017729282e-13 - 1.7967326421342023e-5*x - 8983.663173028655*x^2 + 997.5710159206409*x^3 + 250000.*x^4),
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

		todo!("re-add Roots for glam::* impls")
		// assert_eq!(
		// 	Vec::from_iter(a_pos.to_poly().when_eq(b_pos.to_poly())
		// 		.into_ranges(Time::ZERO)
		// 		.inclusive()),
		// 	[(2*SEC - 83333333*NANOSEC, 2*SEC + 83333333*NANOSEC)]
		// );
	}
	
	#[test]
	fn precise() {
		// https://www.desmos.com/calculator/1z97cqlopx
		
		 // Basic Check:
		let a = Temporal::new(poly!(|x| -193.99999999999997 + 4.481238217799146e-6*x - 500.*x^2), SEC);
		let b = Temporal::new(poly!(|_| -194.), SEC);
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
			poly!(|x| 0.0036784761334161292 + 1.1687626970174242e-7*x + 0.*x^2),
			poly!(|x| -182.00000057575835 - 7.537214753878195e-7*x - 500.*x^2)
		], SEC);
		let b = Temporal::new([
			poly!(|_| -3.8808053943969956e-5),
			poly!(|_| -193.99999999999997)
		], SEC);
		let dis = Temporal::new(poly!(|_| 12.), SEC);
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