//! ...

use std::cmp::Ordering;
use time::{Time, TimeUnit};
use crate::change::{LinearValue, Scalar};

pub type Degree = usize;

/// `a + (b + (c + (d + ..) (t+2)/3) (t+1)/2) t`
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum FluxPolynomial<T: LinearValue> {
	Constant  ([T; 1]),
	Linear    ([T; 2]),
	Quadratic ([T; 3]),
	Cubic     ([T; 4]),
	Quartic   ([T; 5]),
}

impl<T: LinearValue> FluxPolynomial<T> {
	pub fn degree(&self) -> Degree {
		match self {
			Self::Constant(..)  => 0,
			Self::Linear(..)    => 1,
			Self::Quadratic(..) => 2,
			Self::Cubic(..)     => 3,
			Self::Quartic(..)   => 4,
		}
	}
	
	// !!! Derivative is:
	// `a + (b + (c + (d + ..) (t+2)/3) (t+1)/2) t`
	// => (b + c/2 + d/3 + ..) + ((c + d/2 + ..) + (d + ..) (t+1)/2) t
	
	pub fn term(&self, term_index: Degree) -> T {
		match self {
			Self::Constant (list) => list[term_index],
			Self::Linear   (list) => list[term_index],
			Self::Quadratic(list) => list[term_index],
			Self::Cubic    (list) => list[term_index],
			Self::Quartic  (list) => list[term_index],
		}
	}
	
	pub fn add_term(mut self, term_index: Degree, term: T) -> Result<Self, Self> {
		match term_index.cmp(&(self.degree() + 1)) {
			Ordering::Less => {
				match self {
					Self::Constant (ref mut list) => list[term_index] = list[term_index] + term,
					Self::Linear   (ref mut list) => list[term_index] = list[term_index] + term,
					Self::Quadratic(ref mut list) => list[term_index] = list[term_index] + term,
					Self::Cubic    (ref mut list) => list[term_index] = list[term_index] + term,
					Self::Quartic  (ref mut list) => list[term_index] = list[term_index] + term,
				}
				Ok(self)
			},
			Ordering::Equal => match self {
				Self::Constant ([a         ]) => Ok(Self::Linear   ([a,          term])),
				Self::Linear   ([a, b      ]) => Ok(Self::Quadratic([a, b,       term])),
				Self::Quadratic([a, b, c   ]) => Ok(Self::Cubic    ([a, b, c,    term])),
				Self::Cubic    ([a, b, c, d]) => Ok(Self::Quartic  ([a, b, c, d, term])),
				Self::Quartic(..) => Err(self),
			},
			Ordering::Greater => Err(self),
		}
	}
	
	pub fn at(&self, time: Time) -> T {
		let t = (time >> TimeUnit::Nanosecs) as f64;
		let b_time = 1;
		let c_time = 1;
		let d_time = 1;
		let e_time = 1;
		match *self {
			Self::Constant ([a        ]) => a,
			Self::Linear   ([a,b      ]) => a + b*Scalar(t),
			Self::Quadratic([a,b,c    ]) => a + (b + c*Scalar((t+1.0)/2.0))*Scalar(t),
			Self::Cubic    ([a,b,c,d  ]) => a + (b + (c + d*Scalar((t+2.0)/3.0))*Scalar((t+1.0)/2.0))*Scalar(t),
			Self::Quartic  ([a,b,c,d,e]) => a + (b + (c + (d + e*Scalar((t+3.0)/4.0))*Scalar((t+2.0)/3.0))*Scalar((t+1.0)/2.0))*Scalar(t),
		}
	}
}

impl FluxPolynomial<f64> {
	pub fn roots(self, offset: Time) -> Vec<Time> {
		//! Returns all real roots in the range `0..=u64::MAX` as a sorted list.
		
		let b_time = 1.0;
		let c_time = 1.0;
		let d_time = 1.0;
		let e_time = 1.0;
		
		let polynomial = match self {
			Self::Constant(a) => Polynomial::Constant(a),
			Self::Linear([a, b]) => Polynomial::Linear([
				a,
				b/b_time,
			]),
			Self::Quadratic([a, b, c]) => {
				Polynomial::Quadratic([
					a,
					(b + (c/2.0)/c_time)/b_time,
					(c/2.0)/c_time,
				])
			},
			Self::Cubic([a, b, c, d]) => {
				Polynomial::Cubic([
					a,
					(b     +    (c/2.0 + (d/3.0)/d_time)/c_time)/b_time,
					(c/2.0 + ((d+d/2.0)/3.0)/d_time)/c_time,
					(d/2.0/3.0)/d_time,
				])
			},
			Self::Quartic([a, b, c, d, e]) => {
				Polynomial::Quartic([
					a,
					(b         +           (c/2.0      +          (d/3.0 + (e/4.0)/e_time)/d_time)/c_time)/b_time,
					(c/2.0     +        ((d+d/2.0)/3.0 + ((e+e/2.0+e/3.0)/4.0)/e_time)/d_time)/c_time,
					(d/2.0/3.0 + ((e/2.0+(e+e/2.0)/3.0)/4.0)/e_time)/d_time,
					(e/2.0/3.0/4.0)/e_time,
				])
			},
		};
		
		let mut roots: Vec<Time> = polynomial.roots()
			.iter()
			.filter_map(|&r| {
				let r = r? - ((offset >> TimeUnit::Nanosecs) as f64);
				debug_assert!(!r.is_nan());
				if r < 0.0 || r > (u64::MAX as f64) {
					None
				} else {
					Some((r as u64)*TimeUnit::Nanosecs)
				}
			})
			.collect();
		
		roots.sort_unstable();
		roots
	}
}

/// `a + b x + c x^2 + d x^3 + ..`
#[derive(Debug)]
pub enum Polynomial {
	Constant ([f64; 1]),
	Linear   ([f64; 2]),
	Quadratic([f64; 3]),
	Cubic    ([f64; 4]),
	Quartic  ([f64; 5]),
}

impl Polynomial {
	pub fn roots(self) -> Vec<Option<f64>> {
		const ZERO: f64 = 0.0;
		const ONE:  f64 = 1.0;
		
		fn sqrt(num: f64) -> Option<f64> {
			if num < 0.0 {
				None
			} else {
				Some(num.sqrt())
			}
		}
		fn cbrt(num: f64) -> f64 {
			num.cbrt()
		}
		
		match self {
			Self::Constant(_) => vec![],
			
			 // Downgrade:
			Self::Linear   ([a, b         ]) if b == ZERO => [vec![None], Self::Constant ([a         ]).roots()].concat(),
			Self::Quadratic([a, b, c      ]) if c == ZERO => [vec![None], Self::Linear   ([a, b      ]).roots()].concat(),
			Self::Cubic    ([a, b, c, d   ]) if d == ZERO => [vec![None], Self::Quadratic([a, b, c   ]).roots()].concat(),
			Self::Quartic  ([a, b, c, d, e]) if e == ZERO => [vec![None], Self::Cubic    ([a, b, c, d]).roots()].concat(),
			
			 // Shift Over:
			Self::Linear   ([a, b         ]) if a == ZERO => [vec![Some(ZERO)], Self::Constant ([b         ]).roots()].concat(),
			Self::Quadratic([a, b, c      ]) if a == ZERO => [vec![Some(ZERO)], Self::Linear   ([b, c      ]).roots()].concat(),
			Self::Cubic    ([a, b, c, d   ]) if a == ZERO => [vec![Some(ZERO)], Self::Quadratic([b, c, d   ]).roots()].concat(),
			Self::Quartic  ([a, b, c, d, e]) if a == ZERO => [vec![Some(ZERO)], Self::Cubic    ([b, c, d, e]).roots()].concat(),
			
			 // Linear:
			Self::Linear([a, b]) => vec![Some(-a / b)],
			
			 // Pseudo-linear:
			Self::Quadratic([a, b, c]) if b == ZERO => {
				let r = sqrt(-a / c);
				vec![r, r.map(|r| -r)]
			},
			Self::Cubic([a, b, c, d]) if b == ZERO && c == ZERO => {
				let r = Some(cbrt(-a / d));
				vec![r, r, r]
			},
			Self::Quartic([a, b, c, d, e]) if b == ZERO && c == ZERO && d == ZERO => {
				let r = sqrt(-a / e).map(|r| sqrt(r).unwrap());
				vec![
					r, r.map(|r| -r),
					r, r.map(|r| -r),
				]
			},
			
			 // Quadratic:
			Self::Quadratic([a, b, c]) => {
				let n = -b / (2.0 * c);
				let n1 = sqrt(n*n - a/c);
				vec![n1.map(|n1| n + n1), n1.map(|n1| n - n1)]
			},
			
			 // Depressed Cubic:
			Self::Cubic([a, b, c, d]) if c == ZERO => {
				let p = -a / (2.0 * d);
				let q = -b / (3.0 * d);
				let discriminant = p*p - q*q*q;
				match discriminant.partial_cmp(&0.0) {
					 // 3 Real Roots:
					Some(Ordering::Less) => {
						let sqrt_q = sqrt(q).unwrap();
						let angle = f64::acos(p / (q * sqrt_q)) / 3.0;
						use std::f64::consts::TAU;
						vec![
							Some(2.0 * sqrt_q * f64::cos(angle + TAU*0.0/3.0)),
							Some(2.0 * sqrt_q * f64::cos(angle + TAU*1.0/3.0)),
							Some(2.0 * sqrt_q * f64::cos(angle + TAU*2.0/3.0)),
						]
					},
					
					 // 1 Real Root:
					Some(Ordering::Greater) => {
						let n = sqrt(discriminant).unwrap();
						vec![Some(cbrt(p + n) + cbrt(p - n)), None, None]
					},
					
					_ => unreachable!()
				}
			},
			
			 // Cubic:
			Self::Cubic([a, b, c, d]) => {
				let n = -c / (3.0 * d);
				let depressed_cubic = Self::Cubic([
					a + (n * (b + (n * c * 2.0)/3.0)),
					b + (n * c),
					ZERO,
					d,
				]);
				depressed_cubic.roots()
					.iter()
					.map(|x| x.map(|x| x + n))
					.collect()
			},
			
			 // Biquadratic:
			Self::Quartic([a, b, c, d, e]) if b == ZERO && d == ZERO => {
				let roots: Vec<Option<f64>> = Self::Quadratic([a, c, e])
					.roots()
					.iter().map(|&r| sqrt(r?))
					.collect();
				
				vec![
					roots[0], roots[0].map(|x| -x),
					roots[1], roots[1].map(|x| -x)
				]
			},
			
			 // Depressed Quartic:
			Self::Quartic([a, b, c, d, e]) if d == ZERO => {
				let p = a / e;
				let q = b / (2.0 * e);
				let r = c / (2.0 * e);
				let resolvent_cubic = Self::Cubic([
					-(q * q) / 2.0,
					(r * r) - p,
					2.0 * r,
					ONE
				]);
				match *resolvent_cubic.roots().as_slice() {
					[Some(m), _, _] |
					[_, Some(m), _] |
					[_, _, Some(m)] => {
						let sqrt_2m = sqrt(2.0 * m).unwrap();
						let quad_a = Self::Quadratic([ (q / sqrt_2m) + r + m, -sqrt_2m, ONE]);
						let quad_b = Self::Quadratic([-(q / sqrt_2m) + r + m,  sqrt_2m, ONE]);
						[quad_a.roots(), quad_b.roots()].concat()
					},
					_ => unreachable!()
				}
			},
			
			 // Quartic:
			Self::Quartic([a, b, c, d, e]) => {
				let n = -d / (4.0 * e);
				let depressed_quartic = Self::Quartic([
					a + (n * (b + (n * (c + (n * d * 3.0)/4.0)))),
					b + (n * (c + (n * d)) * 2.0),
					c + (n * d * 3.0)/2.0,
					ZERO,
					e
				]);
				depressed_quartic.roots()
					.iter()
					.map(|x| x.map(|x| x + n))
					.collect()
			},
		}
	}
}

#[cfg(test)]
mod root_tests {
	use super::*;
	
	fn assert_roots(p: Polynomial, mut expected_roots: Vec<Option<f64>>) {
		let mut r = p.roots();
		r.sort_unstable_by(|x, y| x.map(|x| x as u64).cmp(&y.map(|y| y as u64)));
		expected_roots.sort_unstable_by(|x, y| x.map(|x| x as u64).cmp(&y.map(|y| y as u64)));
		for i in 0..r.len() {
			if r[i].is_some() && expected_roots[i].is_some() {
				// println!("{:?} <> {:?}", r[i], expected_roots[i]);
				assert!(f64::abs(r[i].unwrap() - expected_roots[i].unwrap()) < 0.1);
			} else {
				assert_eq!(r[i], expected_roots[i]);
			}
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(Polynomial::Constant([2.0]), vec![]);
		assert_roots(Polynomial::Constant([0.0]), vec![]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Polynomial::Linear([20.0, -4.0]), vec![Some(5.0)]);
		assert_roots(Polynomial::Linear([20.0, 4.0]), vec![Some(-5.0)]);
		assert_roots(Polynomial::Linear([-0.0, 4.0/3.0]), vec![Some(0.0)]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(
			Polynomial::Quadratic([20.0, 4.0, 7.0]),
			vec![None, None]
		);
		assert_roots(
			Polynomial::Quadratic([20.0, 4.0, -7.0]),
			vec![Some(-10.0/7.0), Some(2.0)]
		);
		assert_roots(
			Polynomial::Quadratic([ 40.0/3.0,  2.0/3.0, -17.0/100.0]),
			Polynomial::Quadratic([-40.0/3.0, -2.0/3.0,  17.0/100.0]).roots()
		);
		assert_roots(
			Polynomial::Quadratic([0.0, 4.0/6.0, -17.0/100.0]),
			vec![Some(0.0), Some(200.0/51.0)]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Polynomial::Cubic([6.0, -2077.5, -17000.0/77.0, 6712.0/70.0]),
			vec![Some(-3.64550618348), Some(0.00288720188), Some(5.94514363216)]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Polynomial {
			Polynomial::Cubic([
				a,
				b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s),
				c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + d/(2.0*3.0*s*s*s),
				d/(2.0*3.0*s*s*s)
			])
		}
		assert_roots(
			sum_poly(1000000000.0, 1.0, 1.0, 1.0, -1.0),
			vec![None, None, Some(4591405718.8520711602007)]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0, 1.0, 1.0, 1.0, -1.0),
			vec![None, None, Some(275484343229.41354766068234)]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1.0, 1.0, 1.0, -1.0),
			vec![None, None, Some(16529060593863.10213769318)]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1000.0, 300.0, 10.0, -4.0),
			vec![
				Some(-55432951711728.79099027553),
				Some(-13201315814560.382043758431976512),
				Some(95634267526286.173034033963943),
			]
		);
	}
	
	#[test]
	fn quartic() {
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Polynomial {
			Polynomial::Quartic([
				a,
				b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s),
				c/(2.0*s*s) + d/(2.0*3.0*s*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s),
				d/(2.0*3.0*s*s*s) + e/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s),
				e/(2.0*3.0*4.0*s*s*s*s),
			])
		}
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 7.5, -6.4, -10.0, 7.8, 2.7),
			vec![
				Some(-5830778969957.388671875),
				Some(-51671989204288.265625),
				Some(2846556446843.169921875),
				Some(13056211727396.474609375),
			]
		);
		let t = 1000000000.0*60.0*60.0;
		assert_roots(
			Polynomial::Quartic([7500.0, 0.0, -1000.0/(t*t), 0.0, 27.0/(t*t*t*t)]),
			vec![
				Some(-2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0))),
				Some(-2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0))),
				Some(2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0))),
				Some(2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0))),
			]
		);
	}
}

// /// `Sign(integer, numerator)`, `denominator = u64::MAX`
// #[derive(Debug, Copy, Clone)]
// enum Frac {
// 	Pos(u64, u64),
// 	Neg(u64, u64),
// }
// 
// impl_op!{ a <= b { (Frac, Frac)['total] => match(a, b) {
// 	(Frac::Pos(0, 0), Frac::Neg(0, 0)) |
// 	(Frac::Neg(0, 0), Frac::Pos(0, 0)) => Ordering::Equal,
// 	(Frac::Pos(_, _), Frac::Neg(_, _)) => Ordering::Greater,
// 	(Frac::Neg(_, _), Frac::Pos(_, _)) => Ordering::Less,
// 	(Frac::Pos(i, n), Frac::Pos(j, m)) |
// 	(Frac::Neg(j, m), Frac::Neg(i, n)) => i.cmp(j).then_with(|| n.cmp(m)),
// }}}
// 
// impl_op!{ a + b { (Frac, Frac) => match (a, b) {
// 	(Frac::Pos(i, n), Frac::Pos(j, m)) => {
// 		let (n, can_add) = n.overflowing_add(m);
// 		let i = i + j + if can_add { 1 } else { 0 };
// 		Frac::Pos(i, n)
// 	},
// 	(Frac::Neg(i, n), Frac::Neg(j, m)) => {
// 		let (n, can_add) = n.overflowing_add(m);
// 		let i = i + j + if can_add { 1 } else { 0 };
// 		Frac::Neg(i, n)
// 	},
// 	(Frac::Pos(i, n), Frac::Neg(j, m)) => Frac::Pos(i, n) - Frac::Pos(j, m),
// 	(Frac::Neg(i, n), Frac::Pos(j, m)) => Frac::Neg(i, n) - Frac::Neg(j, m),
// }}}
// 
// impl_op!{ a - b { (Frac, Frac) => match(a, b) {
// 	(Frac::Pos(i, n), Frac::Pos(j, m)) => match a.cmp(&b) {
// 		Ordering::Greater => {
// 			let (n, can_sub) = n.overflowing_sub(m);
// 			Frac::Pos(i - j - if can_sub { 1 } else { 0 }, n)
// 		},
// 		Ordering::Less => {
// 			let (m, can_sub) = m.overflowing_sub(n);
// 			Frac::Neg(j - i - if can_sub { 1 } else { 0 }, m)
// 		},
// 		Ordering::Equal => Frac::Pos(0, 0),
// 	},
// 	(Frac::Neg(i, n), Frac::Neg(j, m)) => match a.cmp(&b) {
// 		Ordering::Greater => {
// 			let (n, can_sub) = n.overflowing_sub(m);
// 			Frac::Neg(i - j - if can_sub { 1 } else { 0 }, n)
// 		},
// 		Ordering::Less => {
// 			let (m, can_sub) = m.overflowing_sub(n);
// 			Frac::Pos(j - i - if can_sub { 1 } else { 0 }, m)
// 		},
// 		Ordering::Equal => Frac::Neg(0, 0),
// 	},
// 	(Frac::Pos(i, n), Frac::Neg(j, m)) => Frac::Pos(i, n) + Frac::Pos(j, m),
// 	(Frac::Neg(i, n), Frac::Pos(j, m)) => Frac::Neg(i, n) + Frac::Neg(j, m),
// }}}
// 
// impl_op!{ a * b { (Frac, Frac) => match(a, b) {
// 	// i n/d * j m/d = i*j + i*n/d + j*m/d + n*m/d/d 
// 	(Frac::Pos(i, n), Frac::Pos(j, m)) => {
// 		let (p, p_can_add) = m.overflowing_mul(i);
// 		let (q, q_can_add) = n.overflowing_mul(j);
// 		let (r, r_can_add) = p.overflowing_add(q);
// 		let (s, s_can_add) = r.overflowing_add(n / ((u64::MAX / m) + 1));
// 		Frac::Pos(
// 			(i * j)
// 			+ if p_can_add { (i - 1) / (u64::MAX / m) } else { 0 } // these are a little wrong, due to the loss of the fractional part
// 			+ if q_can_add { (j - 1) / (u64::MAX / n) } else { 0 } // these are a little wrong, due to the loss of the fractional part
// 			+ if r_can_add { 1 } else { 0 }
// 			+ if s_can_add { 1 } else { 0 },
// 			s
// 		)
// 	}
// }}}