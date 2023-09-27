//! ...

use std::cmp::Ordering;
use std::ops::{Add, Mul, Sub};

use crate::change::{LinearValue, Scalar};

use time::{Time, TimeUnit};

use std::slice::{Iter, IterMut};
use std::vec::IntoIter;

pub type Degree = usize;

/// `a + b x + c x^2 + d x^3 + ..`
#[derive(Debug, Clone)]
pub enum Polynomial<T: LinearValue> {
	Constant ([T; 1]),
	Linear   ([T; 2]),
	Quadratic([T; 3]),
	Cubic    ([T; 4]),
	Quartic  ([T; 5]),
	// Generic(Vec<T>),
}

impl<T: LinearValue> Polynomial<T> {
	pub fn degree(&self) -> Degree {
		match self {
			Self::Constant(..)  => 0,
			Self::Linear(..)    => 1,
			Self::Quadratic(..) => 2,
			Self::Cubic(..)     => 3,
			Self::Quartic(..)   => 4,
		}
	}
	
	pub fn term(&self, degree: Degree) -> Option<&T> {
		match self {
			Self::Constant (list) => list.get(degree),
			Self::Linear   (list) => list.get(degree),
			Self::Quadratic(list) => list.get(degree),
			Self::Cubic    (list) => list.get(degree),
			Self::Quartic  (list) => list.get(degree),
		}
	}
	
	pub fn iter(&self) -> Iter<T> {
		match self {
			Self::Constant (list) => list.iter(),
			Self::Linear   (list) => list.iter(),
			Self::Quadratic(list) => list.iter(),
			Self::Cubic    (list) => list.iter(),
			Self::Quartic  (list) => list.iter(),
		}
	}
	
	pub fn iter_mut(&mut self) -> IterMut<T> {
		match self {
			Self::Constant (list) => list.iter_mut(),
			Self::Linear   (list) => list.iter_mut(),
			Self::Quadratic(list) => list.iter_mut(),
			Self::Cubic    (list) => list.iter_mut(),
			Self::Quartic  (list) => list.iter_mut(),
		}
	}
	
	pub fn into_iter(self) -> IntoIter<T> {
		match self {
			Self::Constant (list) => list.to_vec().into_iter(),
			Self::Linear   (list) => list.to_vec().into_iter(),
			Self::Quadratic(list) => list.to_vec().into_iter(),
			Self::Cubic    (list) => list.to_vec().into_iter(),
			Self::Quartic  (list) => list.to_vec().into_iter(),
		}
	}
	
	pub fn real_roots(&self) -> Result<Vec<f64>, Vec<f64>> {
		//! Returns all real-valued roots of this polynomial in ascending order.
		//! If not all roots are known, `Err` is returned.
		
		let cleanup = |mut roots: Vec<f64>| {
			roots.retain(|r| !r.is_nan());
			roots.sort_unstable_by(|a,b| a.total_cmp(b));
			roots
		};
		
		T::roots(self)
			.map(cleanup)
			.map_err(cleanup)
	}
	
	// !!! Derivative is:
	// `a + (b + (c + (d + ..) (t+2)/3) (t+1)/2) t`
	// => (b + c/2 + d/3 + ..) + ((c + d/2 + ..) + (d + ..) (t+1)/2) t
}

impl<T: LinearValue> Add for Polynomial<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		let (mut augend, addend) = if rhs.degree() > self.degree() {
			(rhs, self)
		} else {
			(self, rhs)
		};
		for (degree, term) in augend.iter_mut().enumerate() {
			if let Some(&addend_term) = addend.term(degree) {
				*term = term.clone() + addend_term;
			} else {
				break
			}
		}
		augend
	}
}

impl<T: LinearValue> Sub for Polynomial<T> {
	type Output = Self;
	fn sub(self, rhs: Self) -> Self {
		self + (rhs * Scalar(-1.0))
	}
}

impl<T: LinearValue> Mul<Scalar> for Polynomial<T> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		for term in self.iter_mut() {
			*term = term.clone() * rhs
		}
		self
	}
}

impl LinearValue for f64 {
	fn roots(poly: &Polynomial<f64>) -> Result<Vec<f64>, Vec<f64>> {
		use Polynomial as P;
		// !!! Are infinity coefficients handled properly (cubic, quartic)?
		Ok(match *poly {
			P::Constant(_) => vec![],
			
			 // Downgrade:
			P::Linear   ([a,b      ]) if b == 0. => f64::roots(&P::Constant ([a      ]))?,
			P::Quadratic([a,b,c    ]) if c == 0. => f64::roots(&P::Linear   ([a,b    ]))?,
			P::Cubic    ([a,b,c,d  ]) if d == 0. => f64::roots(&P::Quadratic([a,b,c  ]))?,
			P::Quartic  ([a,b,c,d,e]) if e == 0. => f64::roots(&P::Cubic    ([a,b,c,d]))?,
			
			 // Shift Over:
			P::Linear   ([a,b      ]) if a == 0. => [vec![0.], f64::roots(&P::Constant ([b      ]))?].concat(),
			P::Quadratic([a,b,c    ]) if a == 0. => [vec![0.], f64::roots(&P::Linear   ([b,c    ]))?].concat(),
			P::Cubic    ([a,b,c,d  ]) if a == 0. => [vec![0.], f64::roots(&P::Quadratic([b,c,d  ]))?].concat(),
			P::Quartic  ([a,b,c,d,e]) if a == 0. => [vec![0.], f64::roots(&P::Cubic    ([b,c,d,e]))?].concat(),
			
			 // Linear:
			P::Linear([a, b]) => vec![-a / b],
			
			 // Pseudo-linear:
			P::Quadratic([a,b,c]) if b == 0. => {
				let r = (-a / c).sqrt();
				vec![r, -r]
			},
			P::Cubic([a,b,c,d]) if (b,c) == (0.,0.) => {
				let r = (-a / d).cbrt();
				vec![r, r, r]
			},
			P::Quartic([a,b,c,d,e]) if (b,c,d) == (0.,0.,0.) => {
				let r = (-a / e).sqrt().sqrt();
				vec![r, r, -r, -r]
			},
			
			 // Quadratic:
			P::Quadratic([a,b,c]) => {
				let n = -b / (2.0 * c);
				let n1 = (n*n - a/c).sqrt();
				vec![n + n1, n - n1]
			},
			
			 // Depressed Cubic:
			P::Cubic([a,b,c,d]) if c == 0. => {
				let p = -a / (2.0 * d);
				let q = -b / (3.0 * d);
				let discriminant = p*p - q*q*q;
				match discriminant.partial_cmp(&0.0) {
					 // 3 Real Roots:
					Some(Ordering::Less) => {
						let sqrt_q = q.sqrt();
						debug_assert!(!sqrt_q.is_nan());
						let angle = f64::acos(p / (q * sqrt_q)) / 3.0;
						use std::f64::consts::TAU;
						vec![
							2.0 * sqrt_q * f64::cos(angle + TAU*0.0/3.0),
							2.0 * sqrt_q * f64::cos(angle + TAU*1.0/3.0),
							2.0 * sqrt_q * f64::cos(angle + TAU*2.0/3.0),
						]
					},
					
					 // 1 Real Root:
					Some(Ordering::Greater) => {
						let n = discriminant.sqrt();
						debug_assert!(!n.is_nan());
						vec![(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN]
					},
					
					_ => unreachable!()
				}
			},
			
			 // Cubic:
			P::Cubic([a,b,c,d]) => {
				let n = -c / (3.0 * d);
				let depressed_cubic = P::Cubic([
					a + (n * (b + (n * c * 2.0)/3.0)),
					b + (n * c),
					0.0,
					d,
				]);
				f64::roots(&depressed_cubic)?
					.into_iter()
					.map(|x| x + n)
					.collect()
			},
			
			 // Biquadratic:
			P::Quartic([a,b,c,d,e]) if (b,d) == (0.,0.) => {
				f64::roots(&P::Quadratic([a,c,e]))?
					.into_iter()
					.map(|r| r.sqrt())
					.flat_map(|r| [r, -r])
					.collect()
			},
			
			 // Depressed Quartic:
			P::Quartic([a,b,c,d,e]) if d == 0. => {
				let p = a / e;
				let q = b / (2.0 * e);
				let r = c / (2.0 * e);
				let resolvent_cubic = P::Cubic([
					-(q * q) / 2.0,
					(r * r) - p,
					2.0 * r,
					1.0
				]);
				let m = *resolvent_cubic.real_roots()?.last().unwrap();
				debug_assert!(m >= 0.0);
				let sqrt_2m = (2.0 * m).sqrt();
				let quad_a = P::Quadratic([ (q / sqrt_2m) + r + m, -sqrt_2m, 1.0]);
				let quad_b = P::Quadratic([-(q / sqrt_2m) + r + m,  sqrt_2m, 1.0]);
				[f64::roots(&quad_a)?, f64::roots(&quad_b)?].concat()
			},
			
			 // Quartic:
			P::Quartic([a,b,c,d,e]) => {
				let n = -d / (4.0 * e);
				let depressed_quartic = P::Quartic([
					a + (n * (b + (n * (c + (n * d * 3.0)/4.0)))),
					b + (n * (c + (n * d)) * 2.0),
					c + (n * d * 3.0)/2.0,
					0.0,
					e
				]);
				f64::roots(&depressed_quartic)?
					.into_iter()
					.map(|x| x + n)
					.collect()
			},
		})
	}
}

/*
impl LinearValue for Vec2<f64> {
	fn roots(polynomial: Polynomial<Self>) -> Vec<Option<f64>> {
		let a = (Polynomial<f64> from Vec2 primary terms).roots();
		let b = (Polynomial<f64> from Vec2 secondary terms).roots();
		sort root lists, return pairs of a and b that match (or are close?)
		// https://www.desmos.com/calculator/usraewodwz
	}
}
*/

#[cfg(test)]
mod root_tests {
	use super::*;
	
	fn assert_roots(p: Polynomial<f64>, mut expected_roots: Vec<f64>) {
		let mut r = p.real_roots().unwrap();
		for i in 0..r.len() {
			// println!("{:?} <> {:?}", r[i], expected_roots[i]);
			assert!((r[i] - expected_roots[i]).abs() < 0.1);
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(Polynomial::Constant([2.0]), vec![]);
		assert_roots(Polynomial::Constant([0.0]), vec![]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Polynomial::Linear([20.0, -4.0]), vec![5.0]);
		assert_roots(Polynomial::Linear([20.0, 4.0]), vec![-5.0]);
		assert_roots(Polynomial::Linear([-0.0, 4.0/3.0]), vec![0.0]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(
			Polynomial::Quadratic([20.0, 4.0, 7.0]),
			vec![]
		);
		assert_roots(
			Polynomial::Quadratic([20.0, 4.0, -7.0]),
			vec![-10.0/7.0, 2.0]
		);
		assert_roots(
			Polynomial::Quadratic([ 40.0/3.0,  2.0/3.0, -17.0/100.0]),
			Polynomial::Quadratic([-40.0/3.0, -2.0/3.0,  17.0/100.0]).real_roots().unwrap()
		);
		assert_roots(
			Polynomial::Quadratic([0.0, 4.0/6.0, -17.0/100.0]),
			vec![0.0, 200.0/51.0]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Polynomial::Cubic([6.0, -2077.5, -17000.0/77.0, 6712.0/70.0]),
			vec![-3.64550618348, 0.00288720188, 5.94514363216]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Polynomial<f64> {
			Polynomial::Cubic([
				a,
				b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s),
				c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + d/(2.0*3.0*s*s*s),
				d/(2.0*3.0*s*s*s)
			])
		}
		assert_roots(
			sum_poly(1000000000.0, 1.0, 1.0, 1.0, -1.0),
			vec![4591405718.8520711602007]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0, 1.0, 1.0, 1.0, -1.0),
			vec![275484343229.41354766068234]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1.0, 1.0, 1.0, -1.0),
			vec![16529060593863.10213769318]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1000.0, 300.0, 10.0, -4.0),
			vec![
				-55432951711728.79099027553,
				-13201315814560.382043758431976512,
				95634267526286.173034033963943,
			]
		);
	}
	
	#[test]
	fn quartic() {
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Polynomial<f64> {
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
				-51671989204288.265625,
				-5830778969957.388671875,
				2846556446843.169921875,
				13056211727396.474609375,
			]
		);
		let t = 1000000000.0*60.0*60.0;
		assert_roots(
			Polynomial::Quartic([7500.0, 0.0, -1000.0/(t*t), 0.0, 27.0/(t*t*t*t)]),
			vec![
				-2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
				-2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
			]
		);
	}
}