//! Summation over time.

use std::cmp::Ordering;
use std::ops::{Add, Index, IndexMut, Mul, Sub};
use crate::{*, kind::*, linear::*, exp::*};

/// Summation over time.
/// 
/// - `Sum<0>`: `a`
/// - `Sum<1>`: `a + bx`
/// - `Sum<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Sum<T, const DEGREE: usize>(T, [T; DEGREE]);

impl<T: Linear, const D: usize> Sum<T, D> {
	pub fn new(value: T, coeffs: [T; D]) -> Self {
		Self(value, coeffs)
	}
	
	pub fn iter(&self) -> impl Iterator<Item = &T> {
		std::iter::once(&self.0).chain(&self.1)
	}
	
	pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
		std::iter::once(&mut self.0).chain(&mut self.1)
	}
	
	pub fn into_iter(self) -> impl Iterator<Item = T> {
		// !!! impl IntoIterator
		std::iter::once(self.0).chain(self.1)
	}
}

impl<T: Linear, const D: usize> From<T> for Sum<T, D> {
	fn from(value: T) -> Self {
		Self {
			0: value,
			..Sum::zero()
		}
	}
}

impl<T, const DEG: usize> Index<usize> for Sum<T, DEG> {
	type Output = T;
	fn index(&self, index: usize) -> &Self::Output {
		if index == 0 {
			&self.0
		} else {
			&self.1[index - 1]
		}
	}
}

impl<T, const D: usize> IndexMut<usize> for Sum<T, D> {
	fn index_mut(&mut self, index: usize) -> &mut Self::Output {
		if index == 0 {
			&mut self.0
		} else {
			&mut self.1[index - 1]
		}
	}
}

impl<T: Linear, const D: usize> Mul<Scalar> for Sum<T, D> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.0 * rhs;
		for i in 0..D {
			self.1[i] = self.1[i] * rhs;
		}
		self
	}
}

impl<T: Linear, const D: usize> FluxKind for Sum<T, D> {
	type Value = T;
	
	type Accum<'a> = SumAccum<'a, Self>;
	
	fn at(&self, time: Scalar) -> Self::Value {
		if time == Scalar(0.) {
			return self.0
		}
		let mut value = T::zero();
		for degree in 1..=D {
			value = value*time + self.1[D - degree]
		}
		value*time + self.0
	}
	
	fn initial_order(&self) -> Option<Ordering>
	where
		Self::Value: PartialOrd
	{
		let order = self.iter()
			.enumerate()
			.find_map(|(degree, x)| {
				let order = x.partial_cmp(&T::zero());
				if order != Some(Ordering::Equal) {
					if degree % 2 == 0 {
						order
					} else {
						order.map(|o| o.reverse())
					}
				} else {
					None
				}
			});
		
		if order.is_none() && self.is_zero() {
			Some(Ordering::Equal)
		} else {
			order
		}
	}
	
	fn zero() -> Self {
		Self(T::zero(), [T::zero(); D])
	}
}

impl<T: Linear> Into<Constant<T>> for Sum<T, 0> {
	fn into(self) -> Constant<T> {
		Constant::from(self[0])
	}
}

impl<T: Linear> From<Constant<T>> for Sum<T, 0> {
	fn from(value: Constant<T>) -> Self {
		Self::from(*value)
	}
}

impl<K: FluxKind, const D: usize> Add<K> for Sum<K::Value, D>
where
	K: Into<Constant<K::Value>>
{
	type Output = Self;
	fn add(mut self, rhs: K) -> Self {
		self.0 = self.0 + rhs.value();
		self
	}
}

/// Used for upgrading the degree of a [`Sum`].
trait SumShiftUp: FluxKind {
	type Up: FluxKind<Value=Self::Value>;
	fn shift_up(self) -> <Self as SumShiftUp>::Up;
}

impl<K: FluxKind> SumShiftUp for K
where
	K: Into<Constant<K::Value>>
{
	type Up = Sum<Self::Value, 1>;
	fn shift_up(self) -> <Self as SumShiftUp>::Up {
		Sum(Self::Value::zero(), [self.value()])
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
		impl<T: Linear> SumShiftUp for Sum<T, { $($num +)+ 0 }> {
			type Up = Sum<T, { $($num +)+ 0 + 1 }>;
			fn shift_up(self) -> <Self as SumShiftUp>::Up {
				let mut sum = Sum::zero();
				sum.1[0] = self.0;
				sum.1[1..].copy_from_slice(self.1.as_slice());
				sum
			}
		}
		impl<T: Linear> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, 0> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				rhs
			}
		}
		impl<T: Linear> Add for Sum<T, { $($num +)+ 0 }> {
			type Output = Self;
			fn add(mut self, rhs: Self) -> Self {
				self.0 = self.0 + rhs.0;
				for i in 0..($($num +)+ 0) {
					self.1[i] = self.1[i] + rhs.1[i];
				}
				self
			}
		}
		impl<T: Linear> Mul for Sum<T, { $($num +)+ 0 }> // Squaring
		where
			T: Mul<Output = T>
		{
			type Output = Sum<T, { 2 * ($($num +)+ 0) }>;
			fn mul(self, rhs: Self) -> Self::Output {
				const SIZE: usize = $($num +)+ 0;
				let Sum(a_value, a_terms) = self;
				let Sum(b_value, b_terms) = rhs;
				let mut terms = [T::zero(); { 2 * SIZE }];
				for i in 0..SIZE {
					terms[i] = terms[i] + a_terms[i]*b_value + a_value*b_terms[i];
					for j in 0..SIZE {
						terms[i+j+1] = terms[i+j+1] + a_terms[i]*b_terms[j];
					}
				}
				Sum(a_value*b_value, terms)
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
		impl<T: Linear> Add<Sum<T, $a>> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(mut self, rhs: Sum<T, $a>) -> Self::Output {
				self.0 = self.0 + rhs.0;
				for i in 0..($a) {
					self.1[i] = self.1[i] + rhs.1[i];
				}
				self
			}
		}
		impl<T: Linear> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, $a> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, mut rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				rhs.0 = rhs.0 + self.0;
				for i in 0..($a) {
					rhs.1[i] = rhs.1[i] + self.1[i];
				}
				rhs
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);

impl Roots for Sum<f64, 0> {
	type Output = [f64; 0];
	type IntoIter = <[f64; 0] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		[]
	}
}

impl Roots for Sum<f64, 1> {
	type Output = [f64; 1];
	type IntoIter = <[f64; 1] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let Sum(a, [b]) = self;
		[-a / b]
	}
}

impl Roots for Sum<f64, 2> {
	type Output = [f64; 2];
	type IntoIter = <[f64; 2] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let Sum(a, [b, c]) = self;
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

impl Roots for Sum<f64, 3> {
	type Output = [f64; 3];
	type IntoIter = <[f64; 3] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let Sum(a, [b, c, d]) = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x   * c) / b).abs() < 1e-5 && // ??? Adjust as needed
			((x*x * d) / b).abs() < 1e-16
		) {
			let [y, z] = Sum(b, [c, d]).roots();
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
			let [x, y] = Sum(a, [b, c]).roots();
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
			let discriminant = q.mul_add(-q*q, p*p);
			// if discriminant.abs() < f64::EPSILON {
			// 	discriminant = 0.;
			// }
			return match discriminant.partial_cmp(&0.) {
				 // 3 Real Roots:
				Some(Ordering::Less) => {
					let sqrt_q = q.sqrt();
					debug_assert!(!sqrt_q.is_nan());
					let angle = f64::acos(p / (q * sqrt_q));
					use std::f64::consts::TAU;
					[
						2. * sqrt_q * f64::cos((angle         ) / 3.),
						2. * sqrt_q * f64::cos((angle +    TAU) / 3.),
						2. * sqrt_q * f64::cos((angle + 2.*TAU) / 3.),
					] 
				},
				
				 // 1 Real Root:
				Some(Ordering::Greater) => {
					let n = discriminant.sqrt();
					debug_assert!(!n.is_nan());
					[(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN]
				},
				
				 // 2 Real Roots:
				_ => {
					let r = q / p;
					[-r, -r, 2.*r]
				}
			}
		}
		
		 // General Cubic:
		n /= 3.;
		let depressed_cubic = Sum(
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

impl Roots for Sum<f64, 4> {
	type Output = [f64; 4];
	type IntoIter = <[f64; 4] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let Sum(a, [b, c, d, e]) = self;
		
		 // Weak Constant:
		let x = -a / b;
		if x.is_nan() || (
			((x     * c) / b).abs() < 1e-4  && // ??? Adjust as needed
			((x*x   * d) / b).abs() < 1e-12 &&
			((x*x*x * e) / b).abs() < 1e-20
		) {
			let [y, z, w] = Sum(b, [c, d, e]).roots();
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
			let [x, y, z] = Sum(a, [b, c, d]).roots();
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
				let [x, y] = Sum(a, [c, e]).roots();
				let (x, y) = (x.sqrt(), y.sqrt());
				return [-x, x, -y, y];
			}
			
			 // Depressed Quartic:
			let p = a / e;
			let q = b / (2. * e);
			let r = c / (2. * e);
			let resolvent_cubic = Poly::from(Sum(
				-(q * q) / 2.,
				[
					r.mul_add(r, -p),
					2. * r,
					1.,
				]
			));
			let mut m = resolvent_cubic.real_roots()
				.find(|&r| r > 0.)
				.unwrap_or_else(|| panic!(
					"expected a positive root from: {:?}",
					resolvent_cubic
				));
			let y = (*resolvent_cubic).at(Scalar(m))
				/ m.mul_add(m.mul_add(3., 4.*r), resolvent_cubic.1[0]);
			if m > y && y.is_finite() {
				m -= y; // Newton-Raphson step
			}
			let sqrt_2m = (2. * m).sqrt();
			let [x, y] = Sum( (q / sqrt_2m) + r + m, [-sqrt_2m, 1.]).roots();
			let [z, w] = Sum(-(q / sqrt_2m) + r + m, [ sqrt_2m, 1.]).roots();
			return [x, y, z, w]
		}
		
		 // General Quartic:
		n /= 4.;
		let depressed_quartic = Sum(
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

#[cfg(feature = "glam")]
impl<const D: usize> Roots for Sum<glam::DVec2, D>
where
	Sum<f64, D>: FluxKind<Value=f64> + Roots
{
	type Output = [f64; D];
	type IntoIter = <[f64; D] as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let mut root_list = [f64::NAN; D];
		let mut root_count = 0;
		
		let mut x_poly = Sum::zero();
		let mut y_poly = Sum::zero();
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

impl<T: Linear, const D: usize> Roots for Sum<Exp<T>, D>
where
	Sum<T, D>: FluxKind<Value=T> + Roots
{
	type Output = <Sum<T, D> as Roots>::Output;
	type IntoIter = <<Self as Roots>::Output as IntoIterator>::IntoIter;
	fn roots(self) -> <Self as Roots>::Output {
		let mut b_poly = Sum::zero();
		let mut iter = self.into_iter();
		for item in b_poly.iter_mut() {
			*item = iter.next().unwrap().0;
		}
		b_poly.roots()
	}
}

/// Nested summation change accumulator.
pub struct SumAccum<'a, K: FluxKind> {
	poly: &'a mut K,
	depth: usize,
	time: time::Time,
}

impl<'a, K: FluxKind> FluxAccum<'a, K> for SumAccum<'a, K> {
	fn from_kind(poly: &'a mut K, depth: usize, time: time::Time) -> Self {
		Self { poly, depth, time }
	}
}

impl<K: FluxKind, V: Flux> Add<Change<'_, V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>
{
	type Output = Self;
	fn add(mut self, rhs: Change<'_, V>) -> Self {
		<(K, V::Kind)>::eval(&mut self, 1., rhs.rate, rhs.unit);
		self
	}
}

impl<K: FluxKind, V: Flux> Sub<Change<'_, V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>
{
	type Output = Self;
	fn sub(mut self, rhs: Change<'_, V>) -> Self {
		<(K, V::Kind)>::eval(&mut self, -1., rhs.rate, rhs.unit);
		self
	}
}

/// Used to remove redundant trait bounds.
#[doc(hidden)]
pub trait SumAccumHelper<A: FluxKind, B: FluxKind> {
	fn eval<V: Flux<Kind=B>>(
		kind: &mut SumAccum<'_, A>,
		scalar: f64,
		value: &V,
		unit: time::Time,
	);
}

impl<A, B> SumAccumHelper<A, B> for (A, B)
where
	A: FluxKind,
	B: FluxKind<Value=A::Value> + SumShiftUp,
	A: Add<B, Output=A> + Add<<B as SumShiftUp>::Up, Output=A>,
{
	fn eval<V: Flux<Kind=B>>(
		kind: &mut SumAccum<'_, A>,
		scalar: f64,
		flux: &V,
		unit: time::Time,
	) {
		let mut sub_poly = B::from(flux.value(kind.time));
		flux.change(B::Accum::from_kind(
			&mut sub_poly,
			kind.depth + 1,
			kind.time,
		));
		let time_scale = unit.as_secs_f64().recip();
		*kind.poly = *kind.poly + (sub_poly.shift_up()
			* Scalar((time_scale / ((kind.depth + 1) as f64)) * scalar));
		// https://www.desmos.com/calculator/mhlpjakz32
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn add() {
		let a = Sum(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Sum(7.1, [5.9, 3.1]);
		assert_eq!(a + b, Sum(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
		assert_eq!(b + a, Sum(8.6, [8.4, 6.300000000000001, 4.5, 5.7]));
	}
	
	#[test]
	fn sub() {
		let a = Sum(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Sum(7.1, [5.9, 3.1]);
		assert_eq!(ops::Sub::sub(a, b), Sum(
			-5.6,
			[-3.4000000000000004, 0.10000000000000009, 4.5, 5.7]
		));
		assert_eq!(ops::Sub::sub(b, a), Sum(
			5.6,
			[3.4000000000000004, -0.10000000000000009, -4.5, -5.7]
		));
	}
	
	#[test]
	fn scalar_mul() {
		let a = Sum(1.5, [2.5, 3.2, 4.5, 5.7]);
		let b = Sum(7.1, [5.9, 3.1]);
		assert_eq!(a * Scalar(1.5), Sum(2.25, [3.75, 4.800000000000001, 6.75, 8.55]));
		assert_eq!(b * Scalar(1.5), Sum(10.649999999999999, [8.850000000000001, 4.65]));
	}
	
	fn assert_roots<K: FluxKind>(p: K, expected_roots: &[f64])
	where
		K: Roots
	{
		let mut r = Poly::from(p).real_roots().into_iter()
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
				(x - y).abs() < 0.1,
				"{:?} vs {:?}",
				x, y
			);
		}
	}
	
	#[test]
	fn constant() {
		assert_roots(Sum(2., []), &[]);
		assert_roots(Sum(0., []), &[]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Sum(20., [-4.]), &[5.]);
		assert_roots(Sum(20., [4.]), &[-5.]);
		assert_roots(Sum(-0., [4./3.]), &[0.]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(Sum(20., [4., 7.]), &[]);
		assert_roots(
			Sum(20., [4., -7.]),
			&[-10./7., 2.]
		);
		let r = Poly::from(Sum(-40./3., [-2./3., 17./100.])).real_roots()
			.into_iter().collect::<Vec<_>>();
		assert_roots(
			Sum(40./3., [2./3., -17./100.]),
			r.as_slice()
		);
		assert_roots(
			Sum(0., [4./6., -17./100.]),
			&[0., 200./51.]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Sum(6., [-2077.5, -17000./77., 6712./70.]),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		assert_roots(
			Sum(-3.6350710512584225e-33, [
				0.10240000000000019,
				-0.5455003017809628,
				1.
			]),
			&[3.54987407349455e-32]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Sum<f64, 3> {
			Sum(a, [
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
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Sum<f64, 4> {
			Sum(a, [
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
			Sum(
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
			Sum(
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
			Sum(-35.99999999999999, [8.046918796812349e-6, 2999.99988981303, -54772.25541522836, 250000.0]),
			&[-0.0677022, 0.177247]
		);
	}
	
	#[cfg(feature = "glam")]
	#[test]
	fn vec() {
		#[derive(PartialEq)]
		#[flux(Sum<glam::DVec2, 1> = {value} + spd.per(time::SEC), crate = "crate")]
		struct Pos {
			value: glam::IVec2,
			spd: Spd,
		}
		
		#[derive(PartialEq)]
		#[flux(Sum<glam::DVec2, 0> = {value}, crate = "crate")]
		struct Spd {
			value: glam::IVec2
		}
		
		let a_pos = Pos {
			value: glam::IVec2::new(10, 10),
			spd: Spd {
				value: glam::IVec2::new(6, 4)
			}
		}.to_flux(time::Time::default());
		
		let b_pos = Pos {
			value: glam::IVec2::new(14, 18),
			spd: Spd {
				value: glam::IVec2::new(4, 0)
			}
		}.to_flux(time::Time::default());
		
		// println!("{:?}", a_pos.poly() - b_pos.poly());
		
		assert_eq!(
			a_pos.when_eq(&b_pos).collect::<Vec<time::Time>>(),
			[2*time::SEC]
		);
	}
}