//! Polynomials.

use std::cmp::Ordering;
use std::ops::{Add, IndexMut, Mul, Shr, Sub};
use std::slice::{Iter, IterMut};

use crate::flux::{Deg, DegShift, FluxKind};
use crate::linear::*;

/// A polynomial in standard form; e.g. `a + b x + c x^2 + d x^3`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<K: FluxKind>(pub K::Linear, pub K::Coeffs);

impl<K: FluxKind> Poly<K> {
	pub fn constant(&self) -> K::Linear {
		self.0
	}
	
	pub fn coeff(&self, coeff_index: usize) -> Option<K> {
		self.1.as_ref().get(coeff_index).copied()
	}
	
	pub fn coeff_iter(&self) -> Iter<K> {
		self.1.as_ref().iter()
	}
	
	pub fn coeff_iter_mut(&mut self) -> IterMut<K> {
		self.1.as_mut().iter_mut()
	}
	
	pub fn real_roots(self) -> Result<RootList, RootList>
	where
		K: Roots
	{
		//! Returns all real-valued roots of this polynomial in ascending order.
		//! If not all roots are known, `Err` is returned.
		
		let cleanup = |roots: RootList| {
			let mut roots = roots.into_vec();
			roots.retain(|r| !r.is_nan());
			roots.sort_unstable_by(|a,b| a.total_cmp(b));
			roots.into_boxed_slice()
		};
		
		K::roots(self)
			.map(cleanup)
			.map_err(cleanup)
	}
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Self(K::Linear::zero(), K::zero_coeffs())
	}
}

impl<A, B> Add<Poly<B>> for Poly<A>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Linear=A::Linear>,
	<A as Add<B>>::Output: FluxKind<Linear=A::Linear>,
	A::Linear: Add<B::Linear, Output=A::Linear>,
{
	type Output = Poly<<A as Add<B>>::Output>;
	fn add(self, rhs: Poly<B>) -> Self::Output {
		let mut poly = Poly::default();
		poly.0 = self.constant() + rhs.constant();
		let mut a = self.coeff_iter();
		let mut b = rhs.coeff_iter();
		for i in 0..<<A as Add<B>>::Output as FluxKind>::DEGREE {
			let coeff = IndexMut::index_mut(&mut poly.1, i);
			*coeff = *a.next().unwrap_or(&A::zero())
				+ *b.next().unwrap_or(&B::zero())
		}
		poly
	}
}

impl<A, B> Sub<Poly<B>> for Poly<A>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Linear=A::Linear>,
	<A as Add<B>>::Output: FluxKind<Linear=A::Linear>,
	A::Linear: Add<B::Linear, Output=A::Linear>,
{
	type Output = Poly<<A as Add<B>>::Output>;
	fn sub(self, rhs: Poly<B>) -> Self::Output {
		self + (rhs * Scalar(-1.0))
	}
}

impl<K: FluxKind> Mul<Scalar> for Poly<K> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.constant() * rhs;
		for term in self.coeff_iter_mut() {
			*term = *term * rhs
		}
		self
	}
}

impl<K: FluxKind> Shr<DegShift> for Poly<K>
where
	K: Shr<DegShift>,
	<K as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
	K: From<K::Linear>,
{
	type Output = Poly<<K as Shr<DegShift>>::Output>;
	fn shr(self, rhs: DegShift) -> Self::Output {
		let mut poly = Poly::default();
		let coeff = IndexMut::index_mut(&mut poly.1, 0);
		*coeff = K::from(self.0) >> rhs;
		for d in 0..K::DEGREE {
			let coeff = IndexMut::index_mut(&mut poly.1, d + 1);
			*coeff = self.1[d] >> rhs;
		}
		poly
	}
}

pub type RootList = Box<[f64]>;

/// Roots of a [`Poly`]nomial.
pub trait Roots: FluxKind {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList>;
}

impl Roots for Deg<f64, 0> {
	fn roots(_: Poly<Self>) -> Result<RootList, RootList> {
		Ok([].into())
	}
}

impl Roots for Deg<f64, 1> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a,[Deg(b)]) = poly;
		if b == 0. {
			return Deg::<f64, 0>::roots(Poly(a,[]))
		}
		Ok([-a / b].into())
	}
}

impl Roots for Deg<f64, 2> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a,[Deg(b),Deg(c)]) = poly;
		if c == 0. {
			return Deg::<f64, 1>::roots(Poly(a,[Deg(b)]))
		}
		if a == 0. {
			let roots = [[0.].into(), Deg::<f64, 1>::roots(Poly(b,[Deg(c)]))?];
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

impl Roots for Deg<f64, 3> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a,[Deg(b),Deg(c),Deg(d)]) = poly;
		if d == 0. {
			return Deg::<f64, 2>::roots(Poly(a,[Deg(b),Deg(c)]))
		}
		if a == 0. {
			let roots = [[0.].into(), Deg::<f64, 2>::roots(Poly(b,[Deg(c),Deg(d)]))?];
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
		let depressed_cubic = Poly::<Deg<f64, 3>>(
			a + (n * (b + (n * c * (2.0/3.0)))),
			[
				Deg(b + (n * c)),
				Deg(0.0),
				Deg(d),
			]
		);
		Ok(Roots::roots(depressed_cubic)?
			.into_iter()
			.map(|x| x + n)
			.collect())
	}
}

impl Roots for Deg<f64, 4> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a,[Deg(b),Deg(c),Deg(d),Deg(e)]) = poly;
		if e == 0. {
			return Deg::<f64, 3>::roots(Poly(a,[Deg(b),Deg(c),Deg(d)]))
		}
		if a == 0. {
			let roots = [[0.].into(), Deg::<f64, 3>::roots(Poly(b,[Deg(c),Deg(d),Deg(e)]))?];
			return Ok(roots.concat().into_boxed_slice())
		}
		
		 // Pseudo-linear:
		if (b,c,d) == (0.,0.,0.) {
			let r = (-a / e).sqrt().sqrt();
			return Ok([r, r, -r, -r].into())
		}
		
		 // Biquadratic:
		if (b,d) == (0.,0.) {
			return Ok(Deg::<f64, 2>::roots(Poly(a,[Deg(c),Deg(e)]))?
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
			let resolvent_cubic = Poly::<Deg<f64, 3>>(
				-(q * q) / 2.0,
				[
					Deg((r * r) - p),
					Deg(2.0 * r),
					Deg(1.0),
				]
			);
			let m = *resolvent_cubic.real_roots()?.last().unwrap();
			debug_assert!(m >= 0.0);
			let sqrt_2m = (2.0 * m).sqrt();
			let quad_a = Poly::<Deg<f64, 2>>( (q / sqrt_2m) + r + m, [Deg(-sqrt_2m), Deg(1.0)]);
			let quad_b = Poly::<Deg<f64, 2>>(-(q / sqrt_2m) + r + m, [Deg( sqrt_2m), Deg(1.0)]);
			return Ok([Roots::roots(quad_a)?, Roots::roots(quad_b)?].concat()
				.into_boxed_slice())
		}
		
		 // General Quartic:
		let n = -d / (4.0 * e);
		let depressed_quartic = Poly::<Deg<f64, 4>>(
			a + (n * (b + (n * (c + (n * d * (3.0/4.0)))))),
			[
				Deg(b + (n * (c + (n * d)) * 2.0)),
				Deg(c + (n * d * (3.0/2.0))),
				Deg(0.0),
				Deg(e),
			]
		);
		Ok(Roots::roots(depressed_quartic)?
			.into_iter()
			.map(|x| x + n)
			.collect())
	}
}

/*
impl Roots for Deg<Vec2<f64>, K> {
	fn roots(&self) -> Result<RootList, RootList> {
		let a = (Poly from Vec2 primary terms).roots();
		let b = (Poly from Vec2 secondary terms).roots();
		sort root lists, return pairs of a and b that match (or are close?)
		// https://www.desmos.com/calculator/usraewodwz
	}
}
*/

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn add() {
		let a = Poly::<Deg<f64, 4>>(1.5, [Deg(2.5), Deg(3.2), Deg(4.5), Deg(5.7)]);
		let b = Poly::<Deg<f64, 2>>(7.1, [Deg(5.9), Deg(3.1)]);
		assert_eq!(a + b, Poly(8.6, [Deg(8.4), Deg(6.300000000000001), Deg(4.5), Deg(5.7)]));
		assert_eq!(b + a, Poly(8.6, [Deg(8.4), Deg(6.300000000000001), Deg(4.5), Deg(5.7)]));
	}
	
	#[test]
	fn sub() {
		let a = Poly::<Deg<f64, 4>>(1.5, [Deg(2.5), Deg(3.2), Deg(4.5), Deg(5.7)]);
		let b = Poly::<Deg<f64, 2>>(7.1, [Deg(5.9), Deg(3.1)]);
		assert_eq!(a - b, Poly(-5.6, [Deg(-3.4000000000000004), Deg(0.10000000000000009), Deg(4.5), Deg(5.7)]));
		assert_eq!(b - a, Poly( 5.6, [Deg(3.4000000000000004), Deg(-0.10000000000000009), Deg(-4.5), Deg(-5.7)]));
	}
	
	#[test]
	fn scalar_mul() {
		let a = Poly::<Deg<f64, 4>>(1.5, [Deg(2.5), Deg(3.2), Deg(4.5), Deg(5.7)]);
		let b = Poly::<Deg<f64, 2>>(7.1, [Deg(5.9), Deg(3.1)]);
		assert_eq!(a * Scalar(1.5), Poly(2.25, [Deg(3.75), Deg(4.800000000000001), Deg(6.75), Deg(8.55)]));
		assert_eq!(b * Scalar(1.5), Poly(10.649999999999999, [Deg(8.850000000000001), Deg(4.65)]));
	}
	
	fn assert_roots<K: FluxKind>(p: Poly<K>, expected_roots: &[f64])
	where
		K: Roots
	{
		let r = p.real_roots().unwrap();
		for i in 0..r.len() {
			// println!("{:?} <> {:?}", r[i], expected_roots[i]);
			assert!((r[i] - expected_roots[i]).abs() < 0.1);
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(Poly::<Deg<f64, 0>>(2.0, []), &[]);
		assert_roots(Poly::<Deg<f64, 0>>(0.0, []), &[]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Poly::<Deg<f64, 1>>(20.0, [Deg(-4.0)]), &[5.0]);
		assert_roots(Poly::<Deg<f64, 1>>(20.0, [Deg(4.0)]), &[-5.0]);
		assert_roots(Poly::<Deg<f64, 1>>(-0.0, [Deg(4.0/3.0)]), &[0.0]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(
			Poly::<Deg<f64, 2>>(20.0, [Deg(4.0), Deg(7.0)]),
			&[]
		);
		assert_roots(
			Poly::<Deg<f64, 2>>(20.0, [Deg(4.0), Deg(-7.0)]),
			&[-10.0/7.0, 2.0]
		);
		assert_roots(
			Poly::<Deg<f64, 2>>( 40.0/3.0, [Deg( 2.0/3.0), Deg(-17.0/100.0)]),
			Poly::<Deg<f64, 2>>(-40.0/3.0, [Deg(-2.0/3.0), Deg( 17.0/100.0)]).real_roots().unwrap().as_ref()
		);
		assert_roots(
			Poly::<Deg<f64, 2>>(0.0, [Deg(4.0/6.0), Deg(-17.0/100.0)]),
			&[0.0, 200.0/51.0]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Poly::<Deg<f64, 3>>(6.0, [Deg(-2077.5), Deg(-17000.0/77.0), Deg(6712.0/70.0)]),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Poly<Deg<f64, 3>> {
			Poly(a, [
				Deg(b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s)),
				Deg(c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + d/(2.0*3.0*s*s*s)),
				Deg(d/(2.0*3.0*s*s*s))
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
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Poly<Deg<f64, 4>> {
			Poly(a, [
				Deg(b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Deg(c/(2.0*s*s) + d/(2.0*3.0*s*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Deg(d/(2.0*3.0*s*s*s) + e/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Deg(e/(2.0*3.0*4.0*s*s*s*s)),
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
			Poly::<Deg<f64, 4>>(7500.0, [Deg(0.0), Deg(-1000.0/(t*t)), Deg(0.0), Deg(27.0/(t*t*t*t))]),
			&[
				-2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
				-2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
			]
		);
	}
}