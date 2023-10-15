//! ...

use std::cmp::Ordering;
use std::ops::{Add, Mul, Sub};
use std::slice::{Iter, IterMut};

use crate::change::{LinearValue, Scalar};
use crate::degree::{Deg, IsDeg, MaxDeg};

/// A polynomial in standard form; e.g. `a + b x + c x^2 + d x^3`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<T: LinearValue, D: IsDeg>(pub T, pub D::Array<T>);

impl<T: LinearValue, D: IsDeg> Poly<T, D> {
	pub fn constant(&self) -> T {
		self.0
	}
	
	pub fn coeff(&self, coeff_index: usize) -> Option<T> {
		self.1.as_ref().get(coeff_index).copied()
	}
	
	pub fn coeff_iter(&self) -> Iter<T> {
		self.1.as_ref().iter()
	}
	
	pub fn coeff_iter_mut(&mut self) -> IterMut<T> {
		self.1.as_mut().iter_mut()
	}
	
	pub fn term(&self, degree: usize) -> Option<T> {
		if degree == 0 {
			Some(self.0)
		} else {
			self.coeff(degree - 1)
		}
	}
	
	pub fn real_roots(self) -> Result<RootList, RootList>
	where
		T: Roots<D>
	{
		//! Returns all real-valued roots of this polynomial in ascending order.
		//! If not all roots are known, `Err` is returned.
		
		let cleanup = |roots: RootList| {
			let mut roots = roots.into_vec();
			roots.retain(|r| !r.is_nan());
			roots.sort_unstable_by(|a,b| a.total_cmp(b));
			roots.into_boxed_slice()
		};
		
		T::roots(self)
			.map(cleanup)
			.map_err(cleanup)
	}
}

impl<T: LinearValue, D: IsDeg> Default for Poly<T, D> {
	fn default() -> Self {
		Self(T::zero(), D::array_from(T::zero()))
	}
}

impl<T: LinearValue, A, B> Add<Poly<T, B>> for Poly<T, A>
where
	A: IsDeg + MaxDeg<B>,
	B: IsDeg,
{
	type Output = Poly<T, <A as MaxDeg<B>>::Max>;
	fn add(self, rhs: Poly<T, B>) -> Self::Output {
		let constant = self.constant() + rhs.constant();
		let mut a = self.coeff_iter();
		let mut b = rhs.coeff_iter();
		let mut coeff_list = <<A as MaxDeg<B>>::Max as IsDeg>::array_from(T::zero());
		for i in 0..<<A as MaxDeg<B>>::Max as IsDeg>::USIZE {
			coeff_list[i] = *a.next().unwrap_or(&T::zero())
				+ *b.next().unwrap_or(&T::zero())
		}
		Poly::<T, <A as MaxDeg<B>>::Max>(constant, coeff_list)
	}
}

impl<T: LinearValue, A, B> Sub<Poly<T, B>> for Poly<T, A>
where
	A: IsDeg + MaxDeg<B>,
	B: IsDeg,
{
	type Output = Poly<T, <A as MaxDeg<B>>::Max>;
	fn sub(self, rhs: Poly<T, B>) -> Self::Output {
		self + (rhs * Scalar(-1.0))
	}
}

impl<T: LinearValue, D: IsDeg> Mul<Scalar> for Poly<T, D> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.constant() * rhs;
		for term in self.coeff_iter_mut() {
			*term = *term * rhs
		}
		self
	}
}

pub type RootList = Box<[f64]>;

/// Roots of a [`Poly`]nomial.
pub trait Roots<D: IsDeg>: LinearValue {
	fn roots(poly: Poly<Self, D>) -> Result<RootList, RootList>;
}

impl Roots<Deg<0>> for f64 {
	fn roots(_: Poly<Self, Deg<0>>) -> Result<RootList, RootList> {
		Ok([].into())
	}
}

impl Roots<Deg<1>> for f64 {
	fn roots(poly: Poly<Self, Deg<1>>) -> Result<RootList, RootList> {
		let Poly(a,[b]) = poly;
		if b == 0. {
			return f64::roots(Poly::<f64, Deg<0>>(a,[]))
		}
		Ok([-a / b].into())
	}
}

impl Roots<Deg<2>> for f64 {
	fn roots(poly: Poly<Self, Deg<2>>) -> Result<RootList, RootList> {
		let Poly(a,[b,c]) = poly;
		if c == 0. {
			return f64::roots(Poly::<f64, Deg<1>>(a,[b]))
		}
		if a == 0. {
			let roots = [[0.].into(), f64::roots(Poly::<f64, Deg<1>>(b,[c]))?];
			return Ok(roots.concat().into_boxed_slice())
		}
		
		 // Pseudo-linear:
		if b == 0. {
			let r = (-a / c).sqrt();
			return Ok([r, -r].into())
		}
		
		 // General Quadratic:
		let n = -b / (2.0 * c);
		let n1 = (n*n - a/c).sqrt();
		Ok([n + n1, n - n1].into())
	}
}

impl Roots<Deg<3>> for f64 {
	fn roots(poly: Poly<Self, Deg<3>>) -> Result<RootList, RootList> {
		let Poly(a,[b,c,d]) = poly;
		if d == 0. {
			return f64::roots(Poly::<f64, Deg<2>>(a,[b,c]))
		}
		if a == 0. {
			let roots = [[0.].into(), f64::roots(Poly::<f64, Deg<2>>(b,[c,d]))?];
			return Ok(roots.concat().into_boxed_slice())
		}
		
		 // Pseudo-linear:
		if (b,c) == (0.,0.) {
			let r = (-a / d).cbrt();
			return Ok([r, r, r].into())
		}
		
		 // Depressed Cubic:
		if c == 0. {
			let p = -a / (2.0 * d);
			let q = -b / (3.0 * d);
			let discriminant = p*p - q*q*q;
			return match discriminant.partial_cmp(&0.0) {
				 // 3 Real Roots:
				Some(Ordering::Less) => {
					let sqrt_q = q.sqrt();
					debug_assert!(!sqrt_q.is_nan());
					let angle = f64::acos(p / (q * sqrt_q)) / 3.0;
					use std::f64::consts::TAU;
					Ok([
						2.0 * sqrt_q * f64::cos(angle + TAU*0.0/3.0),
						2.0 * sqrt_q * f64::cos(angle + TAU*1.0/3.0),
						2.0 * sqrt_q * f64::cos(angle + TAU*2.0/3.0),
					].into())
				},
				
				 // 1 Real Root:
				Some(Ordering::Greater) => {
					let n = discriminant.sqrt();
					debug_assert!(!n.is_nan());
					Ok([(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN].into())
				},
				
				_ => unreachable!()
			}
		}
		
		 // General Cubic:
		let n = -c / (3.0 * d);
		let depressed_cubic = Poly::<f64, Deg<3>>(
			a + (n * (b + (n * c * (2.0/3.0)))),
			[
				b + (n * c),
				0.0,
				d,
			]
		);
		Ok(f64::roots(depressed_cubic)?
			.into_iter()
			.map(|x| x + n)
			.collect())
	}
}

impl Roots<Deg<4>> for f64 {
	fn roots(poly: Poly<Self, Deg<4>>) -> Result<RootList, RootList> {
		let Poly(a,[b,c,d,e]) = poly;
		if e == 0. {
			return f64::roots(Poly::<f64, Deg<3>>(a,[b,c,d]))
		}
		if a == 0. {
			let roots = [[0.].into(), f64::roots(Poly::<f64, Deg<3>>(b,[c,d,e]))?];
			return Ok(roots.concat().into_boxed_slice())
		}
		
		 // Pseudo-linear:
		if (b,c,d) == (0.,0.,0.) {
			let r = (-a / e).sqrt().sqrt();
			return Ok([r, r, -r, -r].into())
		}
		
		 // Biquadratic:
		if (b,d) == (0.,0.) {
			return Ok(f64::roots(Poly::<f64, Deg<2>>(a,[c,e]))?
				.into_iter()
				.map(|r| r.sqrt())
				.flat_map(|r| [r, -r])
				.collect())
		}
		
		 // Depressed Quartic:
		if d == 0. {
			let p = a / e;
			let q = b / (2.0 * e);
			let r = c / (2.0 * e);
			let resolvent_cubic = Poly::<f64, Deg<3>>(
				-(q * q) / 2.0,
				[
					(r * r) - p,
					2.0 * r,
					1.0
				]
			);
			let m = *resolvent_cubic.real_roots()?.last().unwrap();
			debug_assert!(m >= 0.0);
			let sqrt_2m = (2.0 * m).sqrt();
			let quad_a = Poly::<f64, Deg<2>>( (q / sqrt_2m) + r + m, [-sqrt_2m, 1.0]);
			let quad_b = Poly::<f64, Deg<2>>(-(q / sqrt_2m) + r + m, [ sqrt_2m, 1.0]);
			return Ok([f64::roots(quad_a)?, f64::roots(quad_b)?].concat()
				.into_boxed_slice())
		}
		
		 // General Quartic:
		let n = -d / (4.0 * e);
		let depressed_quartic = Poly::<f64, Deg<4>>(
			a + (n * (b + (n * (c + (n * d * (3.0/4.0)))))),
			[
				b + (n * (c + (n * d)) * 2.0),
				c + (n * d * (3.0/2.0)),
				0.0,
				e
			]
		);
		Ok(f64::roots(depressed_quartic)?
			.into_iter()
			.map(|x| x + n)
			.collect())
	}
}

/*
impl Roots for Poly<Vec2<f64>, D> {
	fn roots(&self) -> Result<RootList, RootList> {
		let a = (Polynomial<f64, D> from Vec2 primary terms).roots();
		let b = (Polynomial<f64, D> from Vec2 secondary terms).roots();
		sort root lists, return pairs of a and b that match (or are close?)
		// https://www.desmos.com/calculator/usraewodwz
	}
}
*/

#[cfg(test)]
mod tests {
	use crate::degree::Deg;
	use super::*;
	
	#[test]
	fn add() {
		let a = Poly::<f64, Deg<4>>(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Poly::<f64, Deg<2>>(7.1, [5.9, 3.1]);
		assert_eq!(a + b, Poly(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
		assert_eq!(b + a, Poly(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
	}
	
	#[test]
	fn sub() {
		let a = Poly::<f64, Deg<4>>(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Poly::<f64, Deg<2>>(7.1, [5.9, 3.1]);
		assert_eq!(a - b, Poly(-5.6, [-3.4000000000000004, 0.10000000000000009, 4.5, 5.7]));
		assert_eq!(b - a, Poly( 5.6, [3.4000000000000004, -0.10000000000000009, -4.5, -5.7]));
	}
	
	#[test]
	fn scalar_mul() {
		let a = Poly::<f64, Deg<4>>(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Poly::<f64, Deg<2>>(7.1, [5.9, 3.1]);
		assert_eq!(a * Scalar(1.5), Poly(2.25, [3.75, 4.800000000000001, 6.75, 8.55]));
		assert_eq!(b * Scalar(1.5), Poly(10.649999999999999, [8.850000000000001, 4.65]));
	}
	
	fn assert_roots<D: IsDeg>(p: Poly<f64, D>, expected_roots: &[f64])
	where
		f64: Roots<D>
	{
		let r = p.real_roots().unwrap();
		for i in 0..r.len() {
			// println!("{:?} <> {:?}", r[i], expected_roots[i]);
			assert!((r[i] - expected_roots[i]).abs() < 0.1);
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(Poly::<f64, Deg<0>>(2.0, []), &[]);
		assert_roots(Poly::<f64, Deg<0>>(0.0, []), &[]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Poly::<f64, Deg<1>>(20.0, [-4.0]), &[5.0]);
		assert_roots(Poly::<f64, Deg<1>>(20.0, [4.0]), &[-5.0]);
		assert_roots(Poly::<f64, Deg<1>>(-0.0, [4.0/3.0]), &[0.0]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(
			Poly::<f64, Deg<2>>(20.0, [4.0, 7.0]),
			&[]
		);
		assert_roots(
			Poly::<f64, Deg<2>>(20.0, [4.0, -7.0]),
			&[-10.0/7.0, 2.0]
		);
		assert_roots(
			Poly::<f64, Deg<2>>( 40.0/3.0, [ 2.0/3.0, -17.0/100.0]),
			Poly::<f64, Deg<2>>(-40.0/3.0, [-2.0/3.0,  17.0/100.0]).real_roots().unwrap().as_ref()
		);
		assert_roots(
			Poly::<f64, Deg<2>>(0.0, [4.0/6.0, -17.0/100.0]),
			&[0.0, 200.0/51.0]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Poly::<f64, Deg<3>>(6.0, [-2077.5, -17000.0/77.0, 6712.0/70.0]),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Poly<f64, Deg<3>> {
			Poly(a, [
				b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s),
				c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + d/(2.0*3.0*s*s*s),
				d/(2.0*3.0*s*s*s)
			])
		}
		assert_roots(
			sum_poly(1000000000.0, 1.0, 1.0, 1.0, -1.0),
			&[4591405718.8520711602007]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0, 1.0, 1.0, 1.0, -1.0),
			&[275484343229.41354766068234]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1.0, 1.0, 1.0, -1.0),
			&[16529060593863.10213769318]
		);
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 1000.0, 300.0, 10.0, -4.0),
			&[
				-55432951711728.79099027553,
				-13201315814560.382043758431976512,
				95634267526286.173034033963943,
			]
		);
	}
	
	#[test]
	fn quartic() {
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Poly<f64, Deg<4>> {
			Poly(a, [
				b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s),
				c/(2.0*s*s) + d/(2.0*3.0*s*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s),
				d/(2.0*3.0*s*s*s) + e/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s),
				e/(2.0*3.0*4.0*s*s*s*s),
			])
		}
		assert_roots(
			sum_poly(1000000000.0*60.0*60.0, 7.5, -6.4, -10.0, 7.8, 2.7),
			&[
				-51671989204288.265625,
				-5830778969957.388671875,
				2846556446843.169921875,
				13056211727396.474609375,
			]
		);
		let t = 1000000000.0*60.0*60.0;
		assert_roots(
			Poly::<f64, Deg<4>>(7500.0, [0.0, -1000.0/(t*t), 0.0, 27.0/(t*t*t*t)]),
			&[
				-2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
				-2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
			]
		);
	}
}