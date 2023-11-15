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
		for item in self.iter_mut() {
			*item = *item * rhs;
		}
		self
	}
}

impl<T: Linear, const D: usize> FluxKind for Sum<T, D> {
	type Value = T;
	
	type Accum<'a> = SumAccum<'a, Self>;
	
	fn value(&self) -> Self::Value {
		self[0]
	}
	
	fn initial_order(&self) -> Option<Ordering>
	where
		Self::Value: PartialOrd
	{
		let mut degree = D;
		loop {
			let order = self[degree].partial_cmp(&T::zero());
			if degree == 0 || order != Some(Ordering::Equal) {
				break if degree % 2 == 0 {
					order
				} else {
					order.map(|o| o.reverse())
				}
			}
			degree -= 1;
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
		Self::from(value.value())
	}
}

impl<K: FluxKind, const D: usize> Add<K> for Sum<K::Value, D>
where
	K: Into<Constant<K::Value>>
{
	type Output = Self;
	fn add(mut self, rhs: K) -> Self {
		self[0] = self[0] + rhs.value();
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
				let mut iter = self.into_iter();
				for item in sum.iter_mut().skip(1) {
					*item = iter.next().unwrap();
				}
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
				let mut iter = rhs.into_iter();
				for item in self.iter_mut() {
					*item = *item + iter.next().unwrap_or_else(T::zero);
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
				let mut sum = Sum::zero();
				for i in 0..=SIZE {
				for j in 0..=SIZE {
					sum[i+j] = sum[i+j] + self[i]*rhs[j];
				}}
				sum
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
			fn add(self, rhs: Sum<T, $a>) -> Self::Output {
				let mut sum = Sum::zero();
				let mut iter1 = self.into_iter();
				let mut iter2 = rhs.into_iter();
				for item in sum.iter_mut() {
					*item = iter1.next().unwrap_or_else(T::zero)
						+ iter2.next().unwrap_or_else(T::zero);
				}
				sum
			}
		}
		impl<T: Linear> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, $a> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				let mut sum = Sum::zero();
				let mut iter1 = self.into_iter();
				let mut iter2 = rhs.into_iter();
				for item in sum.iter_mut() {
					*item = iter1.next().unwrap_or_else(T::zero)
						+ iter2.next().unwrap_or_else(T::zero);
				}
				sum
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);

const LIMIT: f64 = 0.01;

impl Roots for Sum<f64, 0> {
	fn roots(self) -> Result<RootList, RootList> {
		Ok([].into())
	}
}

impl Roots for Sum<f64, 1> {
	fn roots(self) -> Result<RootList, RootList> {
		let Sum(a, [b]) = self;
		if b == 0. {
			return Sum(a, []).roots()
		}
		Ok([-a / b].into())
	}
}

impl Roots for Sum<f64, 2> {
	fn roots(self) -> Result<RootList, RootList> {
		let Sum(a, [b, c]) = self;
		if c == 0. {
			return Sum(a, [b]).roots()
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum(b, [c]).roots()?
			];
			return Ok(roots.concat().into())
		}
		
		 // Pseudo-linear:
		if b == 0. {
			let r = (-a / c).sqrt();
			return Ok([r, -r].into())
		}
		
		 // Precision Breakdown:
		if (4.0 * a * c/(b*b)).abs() < LIMIT {
			// https://www.desmos.com/calculator/coxloe79ea
			let roots = [
				[-b / c].into(),
				Sum(a, [b]).roots()?
			];
			return Ok(roots.concat().into())
		}
		
		 // General Quadratic:
		let n = -b / (2.0 * c);
		let n1 = (n*n - a/c).sqrt();
		Ok([n + n1, n - n1].into())
	}
}

impl Roots for Sum<f64, 3> {
	fn roots(self) -> Result<RootList, RootList> {
		let Sum(a, [b, c, d]) = self;
		if d == 0. {
			return Sum(a, [b, c]).roots()
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum(b, [c, d]).roots()?
			];
			return Ok(roots.concat().into())
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
				
				 // This shouldn't really ever run, but just in case:
				_ => {
					let roots = [
						[0.].into(),
						Sum(a, [b, c]).roots()?
					];
					return Ok(roots.concat().into())
				}
			}
		}
		
		 // Precision Breakdown:
		if
			(4.0  * b *   d/(c*c))   .abs() < LIMIT &&
			(27.0 * a * d*d/(c*c*c)) .abs() < LIMIT
			// https://www.desmos.com/calculator/4kgmefjtqz
		{
			let roots = [
				[-c / d].into(),
				Sum(a, [b, c]).roots()?
			];
			return Ok(roots.concat().into())
		}
		
		 // General Cubic:
		let n = -c / (3.0 * d);
		let depressed_cubic = Sum(
			a + (n * (b + (n * c * (2.0/3.0)))),
			[
				b + (n * c),
				0.0,
				d,
			]
		);
		Ok(depressed_cubic.roots()?
			.iter()
			.map(|x| x + n)
			.collect())
	}
}

impl Roots for Sum<f64, 4> {
	fn roots(self) -> Result<RootList, RootList> {
		let Sum(a, [b, c, d, e]) = self;
		if e == 0. {
			return Sum(a, [b, c, d]).roots()
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum(b, [c, d, e]).roots()?
			];
			return Ok(roots.concat().into())
		}
		
		 // Pseudo-linear:
		if (b,c,d) == (0.,0.,0.) {
			let r = (-a / e).sqrt().sqrt();
			return Ok([r, r, -r, -r].into())
		}
		
		 // Biquadratic:
		if (b,d) == (0.,0.) {
			return Ok(Sum(a, [c, e]).roots()?
				.iter()
				.map(|r| r.sqrt())
				.flat_map(|r| [r, -r])
				.collect())
		}
		
		 // Depressed Quartic:
		if d == 0. {
			let p = a / e;
			let q = b / (2.0 * e);
			let r = c / (2.0 * e);
			let resolvent_cubic = Sum(
				-(q * q) / 2.0,
				[
					(r * r) - p,
					2.0 * r,
					1.0,
				]
			);
			let m = *Poly::from(resolvent_cubic).real_roots()
				.unwrap_or_default() 
				.iter().find(|&r| *r >= 0.0)
				.expect("this shouldn't happen, probably a precision issue");
			let sqrt_2m = (2.0 * m).sqrt();
			let quad_a = Sum( (q / sqrt_2m) + r + m, [-sqrt_2m, 1.0]);
			let quad_b = Sum(-(q / sqrt_2m) + r + m, [ sqrt_2m, 1.0]);
			return Ok([Roots::roots(quad_a)?, Roots::roots(quad_b)?].concat().into())
		}
		
		 // Precision Breakdown:
		if
			(4.0   * c *     e/(d*d))     .abs() < LIMIT &&
			(27.0  * b *   e*e/(d*d*d))   .abs() < LIMIT &&
			(256.0 * a * e*e*e/(d*d*d*d)) .abs() < LIMIT
			// https://www.desmos.com/calculator/09upa4bqod
		{
			let roots = [
				[-d / e].into(),
				Sum(a, [b, c, d]).roots()?
			];
			return Ok(roots.concat().into())
		}
		
		 // General Quartic:
		let n = -d / (4.0 * e);
		let depressed_quartic = Sum(
			a + (n * (b + (n * (c + (n * d * (3.0/4.0)))))),
			[
				b + (n * (c + (n * d)) * 2.0),
				c + (n * d * (3.0/2.0)),
				0.0,
				e,
			]
		);
		Ok(depressed_quartic.roots()?
			.iter()
			.map(|x| x + n)
			.collect())
	}
}

#[cfg(feature = "glam")]
impl<const D: usize> Roots for Sum<glam::DVec2, D>
where
	Sum<f64, D>: FluxKind<Value=f64> + Roots
{
	fn roots(self) -> Result<RootList, RootList> {
		let mut x_poly = Sum::zero();
		let mut y_poly = Sum::zero();
		for (index, coeff) in self.into_iter().enumerate() {
			x_poly[index] = coeff.x;
			y_poly[index] = coeff.y;
		}
		
		type Output = fn(RootList) -> Result<RootList, RootList>;
		let (output, mut x_roots, mut y_roots): (Output, _, _) = match (
			x_poly.roots(),
			y_poly.roots()
		) {
			(Ok(x),        Ok(y)       ) => (Ok,  x.to_vec(), y.to_vec()),
			(Ok(x)|Err(x), Ok(y)|Err(y)) => (Err, x.to_vec(), y.to_vec()),
		};
		x_roots.sort_unstable_by(f64::total_cmp);
		y_roots.sort_unstable_by(f64::total_cmp);
		
		let mut x_iter = x_roots.into_iter();
		let mut y_iter = y_roots.into_iter();
		
		let mut x = x_iter.next();
		let mut y = y_iter.next();
		
		let mut root_list = Vec::new();
		
		while let (Some(a), Some(b)) = (x.as_ref(), y.as_ref()) {
			match a.total_cmp(b) {
				Ordering::Less    => x = x_iter.next(),
				Ordering::Greater => y = y_iter.next(),
				Ordering::Equal   => {
					root_list.push(*a);
					x = x_iter.next();
					y = y_iter.next();
					// !!! What if the next value(s) are equal to the current?
					// Should only one of these be advanced?
				}
			}
		}
		
		output(root_list.into())
	}
}

impl<T: Linear, const D: usize> Roots for Sum<Exp<T>, D>
where
	Sum<T, D>: FluxKind<Value=T> + Roots
{
	fn roots(self) -> Result<RootList, RootList> {
		let mut b_poly = Sum::zero();
		let mut iter = self.into_iter();
		for item in b_poly.iter_mut() {
			*item = iter.next().unwrap().0;
		}
		b_poly.roots()
	}
}

/// Nested summation change accumulator.
pub struct SumAccum<'a, K: FluxKind>(FluxAccumKind<'a, K>);

impl<'a, K: FluxKind> FluxAccum<'a, K> for SumAccum<'a, K> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self {
		Self(kind)
	}
}

impl<K: FluxKind> SumAccum<'_, K> {
	fn accum<V: Flux>(mut self, scalar: Scalar, change: Change<'_, V>) -> Self
	where
		(K, V::Kind): SumAccumHelper<K, V::Kind>,
	{
		<(K, V::Kind)>::eval(&mut self.0, scalar, change.rate, change.unit);
		self
	}
}

impl<K: FluxKind, V: Flux> Add<Change<'_, V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>
{
	type Output = Self;
	fn add(self, rhs: Change<'_, V>) -> Self {
		self.accum(Scalar(1.), rhs)
	}
}

impl<K: FluxKind, V: Flux> Sub<Change<'_, V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>
{
	type Output = Self;
	fn sub(self, rhs: Change<'_, V>) -> Self {
		self.accum(Scalar(-1.), rhs)
	}
}

/// Used to remove redundant trait bounds.
#[doc(hidden)]
pub trait SumAccumHelper<A: FluxKind, B: FluxKind> {
	fn eval<V: Flux<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		value: &V,
		unit: TimeUnit,
	);
}

impl<A, B> SumAccumHelper<A, B> for (A, B)
where
	A: FluxKind,
	B: FluxKind<Value=A::Value> + SumShiftUp,
	A: Add<B, Output=A> + Add<<B as SumShiftUp>::Up, Output=A>,
{
	fn eval<V: Flux<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		flux: &V,
		unit: TimeUnit,
	) {
		match kind {
			FluxAccumKind::Value { value, depth, time, base_time } => {
				let mut sub_value = flux.value(*base_time);
				flux.change(B::Accum::from_kind(FluxAccumKind::Value {
					value: &mut sub_value,
					depth: *depth + 1,
					time: *time,
					base_time: *base_time,
				}));
				let depth = *depth as f64;
				let time_scale = (time.as_secs_f64() - base_time.as_secs_f64())
					/ (1*unit).as_secs_f64();
				**value = **value
					+ (sub_value * Scalar((time_scale+depth) / (depth+1.)) * scalar);
			},
			FluxAccumKind::Poly { poly, depth, time, base_time } => {
				let mut sub_poly = Poly::with_value(flux.value(*time));
				flux.change(B::Accum::from_kind(FluxAccumKind::Poly {
					poly: &mut sub_poly,
					depth: *depth + 1,
					time: *time,
					base_time: *base_time,
				}));
				let sup_poly = Poly::from(sub_poly.shift_up());
				let depth = *depth as f64;
				let time_scale = (1*unit).as_secs_f64().recip();
				**poly = **poly
					+ (sub_poly * (Scalar(depth / (depth + 1.)) * scalar))
					+ (sup_poly * (Scalar(time_scale / (depth + 1.)) * scalar));
				// https://www.desmos.com/calculator/mhlpjakz32
			},
		}
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
		let r = Poly::from(p).real_roots().unwrap_or_default();
		assert_eq!(
			r.len(), expected_roots.len(),
			"{:?} vs {:?}",
			r, expected_roots
		);
		for i in 0..r.len() {
			assert!(
				(r[i] - expected_roots[i]).abs() < 0.1,
				"{:?} vs {:?}",
				r[i], expected_roots[i]
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
		assert_roots(
			Sum(40./3., [2./3., -17./100.]),
			Poly::from(Sum(-40./3., [-2./3., 17./100.])).real_roots().unwrap().as_ref()
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
	}
	
	#[cfg(feature = "glam")]
	#[test]
	fn vec() {
		#[derive(PartialEq)]
		#[flux(Sum<glam::DVec2, 1> = {value} + spd.per(TimeUnit::Secs), crate = "crate")]
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
		}.to_flux(Time::default());
		
		let b_pos = Pos {
			value: glam::IVec2::new(14, 18),
			spd: Spd {
				value: glam::IVec2::new(4, 0)
			}
		}.to_flux(Time::default());
		
		// println!("{:?}", a_pos.poly() - b_pos.poly());
		
		assert_eq!(
			a_pos.when_eq(&b_pos).collect::<Vec<Time>>(),
			[2*TimeUnit::Secs]
		);
	}
}