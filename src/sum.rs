//! Summation over time.

use std::cmp::Ordering;
use std::ops::{Add, Index, IndexMut, Mul, Sub};
use crate::{*, kind::*, linear::*, exp::*};
use crate::Flux;

/// Summation over time.
/// 
/// - `Sum<0>`: `a`
/// - `Sum<1>`: `a + bx`
/// - `Sum<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug)]
#[repr(C)]
pub struct Sum<T, const DEGREE: usize>(T, [T; DEGREE]);

impl<T: LinearPlus, const D: usize> Sum<T, D> {
	pub fn new(value: T, coeffs: [T; D]) -> Self {
		Self(value, coeffs)
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

impl<T: LinearPlus, const D: usize> From<T> for Sum<T, D> {
	fn from(value: T) -> Self {
		Self::from_value(value.into_inner())
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

impl<T, const D: usize> PartialEq for Sum<T, D>
where
	T: LinearPlus,
	T::Inner: PartialEq,
{
	fn eq(&self, other: &Self) -> bool {
		if !self.0.inner_eq(&other.0) {
			return false
		}
		for i in 0..D {
			if !self.1[i].inner_eq(&other.1[i]) {
				return false
			}
		}
		true
	}
}

impl<T: LinearPlus, const D: usize> Mul<Scalar> for Sum<T, D> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = T::from_inner(self.0.into_inner().mul_scalar(rhs));
		self.1 = self.1.map(|x| T::from_inner(x.into_inner().mul_scalar(rhs)));
		self
	}
}

impl<T: LinearPlus, const D: usize> FluxKind for Sum<T, D> {
	type Value = T;
	
	type Accum<'a> = SumAccum<'a, Self>;
	
	type OutAccum<'a> = SumAccum<'a, Self>;
	
	const DEGREE: usize = D;
	
	fn from_value(value: <Self::Value as LinearPlus>::Inner) -> Self {
		Self(T::from_inner(value), std::array::from_fn(|_| T::zero()))
	}

	fn deriv(mut self) -> Self {
		std::mem::swap(&mut self.0, &mut self.1[0]);
		self.1.rotate_left(1);
		self.1[D-1] = T::zero();
		let mut d = 1.;
		self.1 = self.1.map(|x| {
			d += 1.;
			T::from_inner(x.into_inner().mul_scalar(Scalar::from(d)))
		});
		self
	}
	
	fn as_accum(&mut self, depth: usize, base_time: time::Time, time: time::Time) -> Self::Accum<'_> {
		SumAccum {
			poly: self,
			depth,
			base_time,
			time,
		}
	}
	
	fn eval(&self, time: Scalar) -> <Self::Value as LinearPlus>::Inner {
		if time == Scalar::from(0.) {
			return self.0.clone().into_inner()
		}
		let mut value = <T::Inner as Linear>::zero();
		for degree in 1..=D {
			value = value.mul_scalar(time).add(self.1[D - degree].clone().into_inner());
		}
		value.mul_scalar(time).add(self.0.clone().into_inner())
	}
	
	fn to_time(mut self, time: Scalar) -> Self {
		if time == Scalar::from(0.) {
			return self
		}
		let mut deriv = self.clone();
		self.0 = T::from_inner(deriv.eval(time));
		for degree in 1..=D {
			deriv = deriv.deriv().mul(Scalar::from(1. / (degree as f64)));
			self.1[degree-1] = T::from_inner(deriv.eval(time));
		}
		self
	}
}

impl<T, const DEGREE: usize, const SIZE: usize> Vector<SIZE> for Sum<T, DEGREE>
where
	T: Vector<SIZE>
{
	type Output = Sum<T::Output, DEGREE>;
	fn index(&self, index: usize) -> Self::Output {
		Sum(
			self.0.index(index),
			std::array::from_fn(|i| self.1[i].index(index))
		)
	}
}

#[allow(clippy::from_over_into)]
impl<T: LinearPlus> Into<Constant<T>> for Sum<T, 0> {
	fn into(self) -> Constant<T> {
		Constant::from(self.0)
	}
}

impl<T: LinearPlus> From<Constant<T>> for Sum<T, 0> {
	fn from(Constant(value): Constant<T>) -> Self {
		Self::from(value)
	}
}

impl<A, B, const D: usize> Add<B> for Sum<A, D>
where
	A: LinearPlus,
	B: FluxKind + Into<Constant<B::Value>>,
	B::Value: LinearPlus<Inner = A::Inner>,
{
	type Output = Self;
	fn add(mut self, rhs: B) -> Self {
		self.0 = A::from_inner(self.0.into_inner().add(rhs.value().into_inner()));
		self
	}
}

impl<T: LinearPlus> FluxIntegral for Constant<T> {
	type Integ = Sum<T, 1>;
	fn integ(self) -> Self::Integ {
		let Constant(value) = self;
		Sum::new(T::zero(), [value])
	}
}

impl<T: LinearPlus> FluxIntegral for Sum<T, 0> {
	type Integ = Sum<T, 1>;
	fn integ(self) -> Self::Integ {
		let Sum(value, _) = self;
		Sum(T::zero(), [value])
	}
}

/// Used for upgrading the degree of a [`Sum`].
trait SumShiftUp: FluxKind {
	type Up: FluxKind<Value: LinearPlus<Inner = KindLinear<Self>>>;
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
		impl<T: LinearPlus> SumShiftUp for Sum<T, { $($num +)+ 0 }> {
			type Up = Sum<T, { $($num +)+ 0 + 1 }>;
			fn shift_up(self) -> <Self as SumShiftUp>::Up {
				debug_assert_eq!(
					std::mem::size_of::<Self>(),
					std::mem::size_of::<[T; { $($num +)+ 0 + 1 }]>(),
				);
				Sum(
					T::zero(),
					unsafe {
						// SAFETY: I don't know if the memory layout of arrays
						// is guaranteed, but it's simple & fast. Sorry boss.
						std::mem::transmute_copy(&self)
					}
				)
			}
		}
		impl<T: LinearPlus> FluxIntegral for Sum<T, { $($num +)+ 0 }> {
			type Integ = Sum<T, { $($num +)+ 0 + 1 }>;
			fn integ(self) -> Self::Integ {
				debug_assert_eq!(
					std::mem::size_of::<Self>(),
					std::mem::size_of::<[T; { $($num +)+ 0 + 1 }]>(),
				);
				let mut terms = unsafe {
					// SAFETY: I don't know if the memory layout of arrays
					// is guaranteed, but it's simple & fast. Sorry boss.
					std::mem::transmute_copy::<Self, [T; { $($num +)+ 0 + 1 }]>(&self)
				};
				let mut i = 0.;
				terms = terms.map(|term| {
					i += 1.;
					T::from_inner(term.into_inner().mul_scalar(Scalar::from(1. / i)))
				});
				Sum(T::zero(), terms)
			}
		}
		impl<A, B> Add<Sum<B, { $($num +)+ 0 }>> for Sum<A, 0>
		where
			A: LinearPlus,
			B: LinearPlus<Inner = A::Inner>,
		{
			type Output = Sum<A, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<B, { $($num +)+ 0 }>) -> Self::Output {
				Sum::new(
					A::from_inner(rhs.0.into_inner()),
					rhs.1.map(|x| A::from_inner(LinearPlus::into_inner(x))),
				)
			}
		}
		impl<A, B> Add<Sum<B, { $($num +)+ 0 }>> for Sum<A, { $($num +)+ 0 }>
		where
			A: LinearPlus,
			B: LinearPlus<Inner = A::Inner>,
		{
			type Output = Sum<A, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<B, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				Sum(
					A::from_inner(self.0.into_inner().add(rhs.0.into_inner())),
					std::array::from_fn(|_| unsafe {
						// SAFETY: Sizes of all input & output arrays are equal.
						A::from_inner(a.next().unwrap_unchecked().into_inner()
							.add(b.next().unwrap_unchecked().into_inner()))
					}),
				)
			}
		}
		impl<T> Mul for Sum<T, { $($num +)+ 0 }> // Squaring
		where
			T: LinearPlus,
			T::Inner: Mul<Output = T::Inner>,
		{
			type Output = Sum<T, { 2 * ($($num +)+ 0) }>;
			#[allow(clippy::suspicious_arithmetic_impl)]
			fn mul(self, rhs: Self) -> Self::Output {
				const SIZE: usize = $($num +)+ 0;
				let Sum(a_value, a_terms) = self;
				let Sum(b_value, b_terms) = rhs;
				let a_value = a_value.into_inner();
				let b_value = b_value.into_inner();
				let a_terms = a_terms.map(T::into_inner);
				let b_terms = b_terms.map(T::into_inner);
				let terms = std::array::from_fn(|i| T::from_inner(
					if i < SIZE {
						let mut x = (a_terms[i].clone() * b_value.clone())
							.add(a_value.clone() * b_terms[i].clone());
						for j in 0..i {
							x = x.add(a_terms[j].clone() * b_terms[i-j-1].clone());
						}
						x
					} else {
						let mut x = T::zero().into_inner();
						for j in (i - SIZE)..SIZE {
							x = x.add(a_terms[j].clone() * b_terms[i-j-1].clone());
						}
						x
					}
				));
				Sum(T::from_inner(a_value.into_inner()*b_value.into_inner()), terms)
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
		impl<A, B> Add<Sum<B, $a>> for Sum<A, { $($num +)+ 0 }>
		where
			A: LinearPlus,
			B: LinearPlus<Inner = A::Inner>,
		{
			type Output = Sum<A, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<B, $a>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				Sum(
					A::from_inner(self.0.into_inner().add(rhs.0.into_inner())),
					std::array::from_fn(|i| unsafe {
						// SAFETY: `a` is the same size as the output array, and
						// `b` is the size of `$a` (bad naming).
						if i < $a {
							A::from_inner(a.next().unwrap_unchecked().into_inner()
								.add(b.next().unwrap_unchecked().into_inner()))
						} else {
							a.next().unwrap_unchecked()
						}
					}),
				)
			}
		}
		impl<A, B> Add<Sum<B, { $($num +)+ 0 }>> for Sum<A, $a>
		where
			A: LinearPlus,
			B: LinearPlus<Inner = A::Inner>,
		{
			type Output = Sum<A, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<B, { $($num +)+ 0 }>) -> Self::Output {
				let mut a = self.1.into_iter();
				let mut b = rhs.1.into_iter();
				Sum::new(
					A::from_inner(self.0.into_inner().add(rhs.0.into_inner())),
					std::array::from_fn(|i| unsafe {
						// SAFETY: `b` is the same size as the output array, and
						// `a` is the size of `$a`.
						if i < $a {
							A::from_inner(a.next().unwrap_unchecked().into_inner()
								.add(b.next().unwrap_unchecked().into_inner()))
						} else {
							A::from_inner(b.next().unwrap_unchecked().into_inner())
						}
					}),
				)
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);

impl Roots for Sum<f64, 0> {
	type Output = [f64; 0];
	fn roots(self) -> <Self as Roots>::Output {
		[]
	}
}

impl Roots for Sum<f64, 1> {
	type Output = [f64; 1];
	fn roots(self) -> <Self as Roots>::Output {
		let Sum(a, [b]) = self;
		[-a / b]
	}
}

impl Roots for Sum<f64, 2> {
	type Output = [f64; 2];
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
			let resolvent_cubic = Sum(
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
			let y = resolvent_cubic.eval(Scalar::from(m))
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

impl<A, B, const D: usize> Roots for Sum<Iso<A, B>, D>
where
	A: Linear,
	B: LinearIso<A>,
	Sum<A, D>: Roots,
{
	type Output = <Sum<A, D> as Roots>::Output;
	fn roots(self) -> <Self as Roots>::Output {
		let sum = Sum::<A, D>::new(
			self.0.into_inner(),
			self.1.map(|x| x.into_inner())
		);
		<Sum<A, D> as Roots>::roots(sum)
	}
}

#[cfg(feature = "glam")]
impl<const D: usize> Roots for Sum<glam::DVec2, D>
where
	Sum<f64, D>: FluxKind<Value=f64> + Roots,
	<Sum<f64, D> as Roots>::Output: IntoIterator<Item=f64>,
{
	type Output = [f64; D];
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
	base_time: time::Time,
	time: time::Time,
}

impl<K: FluxKind, V: Flux> Add<Change<&V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>,
	V::Moment: Moment<Flux=V>,
{
	type Output = Self;
	fn add(mut self, rhs: Change<&V>) -> Self {
		<(K, V::Kind)>::eval(&mut self, 1., rhs.rate, rhs.unit);
		self
	}
}

impl<K: FluxKind, V: Flux> Sub<Change<&V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>,
	V::Moment: Moment<Flux=V>,
{
	type Output = Self;
	fn sub(mut self, rhs: Change<&V>) -> Self {
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
		flux: &V,
		unit: time::Time,
	)
	where
		V::Moment: Moment<Flux=V>;
}

impl<A, B> SumAccumHelper<A, B> for (A, B)
where
	A: FluxKind,
	B: FluxKind<Value: LinearPlus<Inner = KindLinear<A>>> + SumShiftUp,
	<B as SumShiftUp>::Up: Mul<Scalar, Output = <B as SumShiftUp>::Up>,
	A: Add<B, Output=A> + Add<<B as SumShiftUp>::Up, Output=A>,
{
	fn eval<V: Flux<Kind=B>>(
		kind: &mut SumAccum<'_, A>,
		scalar: f64,
		flux: &V,
		unit: time::Time,
	)
	where
		V::Moment: Moment<Flux=V>
	{
		let flux_ref = FluxValue::new(flux, kind.base_time);
		let mut sub_poly = B::from_value(flux_ref.eval(kind.time));
		flux_ref.change(sub_poly.as_accum(kind.depth + 1, kind.base_time, kind.time));
		let time_scale = unit.as_secs_f64().recip();
		*kind.poly = kind.poly.clone() + (sub_poly.shift_up()
			* Scalar::from((time_scale / ((kind.depth + 1) as f64)) * scalar));
		// https://www.desmos.com/calculator/mhlpjakz32
	}
}

#[cfg(test)]
mod tests {
	use crate::time::*;
	use crate::pred::Prediction;
	use super::*;
	
	fn real_roots<K>(poly: Poly<K>) -> impl Iterator<Item=f64>
	where
		K: Roots,
		<K as Roots>::Output: IntoIterator<Item=f64>,
	{
		poly.into_inner().roots().into_iter()
			.filter(|r| r.is_finite())
	}
	
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
		assert_eq!(a * Scalar::from(1.5), Sum(2.25, [3.75, 4.800000000000001, 6.75, 8.55]));
		assert_eq!(b * Scalar::from(1.5), Sum(10.649999999999999, [8.850000000000001, 4.65]));
	}
	
	fn assert_roots<K: FluxKind>(p: K, expected_roots: &[f64])
	where
		K: Roots,
		<K as Roots>::Output: IntoIterator<Item=f64>,
	{
		let mut r = real_roots(Poly::from(p)).into_iter()
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
		let r = real_roots(Poly::from(Sum(-40./3., [-2./3., 17./100.])))
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
		assert_eq!(
			real_roots(Poly::from(Sum(
				-236263115684.8131,
				[-9476965815.566229, -95034754.784949, 1.0]
			))).count(),
			3
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
			&[-0.067702232, 0.177246743]
		);
		assert_eq!(Poly::new(Sum(-181.99999999999994, [-7.289202428347473e-6, -500.0]), 1209618449*NANOSEC).at(1209618450*NANOSEC), -181.99999999999994);
		assert_roots(
			Sum(9.094947017729282e-13, [-1.7967326421342023e-5, -8983.663173028655, 997.5710159206409, 250000.0]),
			&[-0.191570016, -1.11113e-8, 9.11132e-9, 0.187579734]
		);
	}
	
	#[cfg(feature = "glam")]
	#[test]
	fn vec() {
		#[flux(
			kind = Sum<Iso<glam::DVec2, glam::IVec2>, 1>,
			value = value,
			change = |c| c + spd.per(SEC),
			crate = crate,
		)]
		#[derive(PartialEq)]
		struct Pos {
			value: Iso<glam::DVec2, glam::IVec2>,
			spd: Spd,
		}
		
		#[flux(
			kind = Sum<Iso<glam::DVec2, glam::IVec2>, 0>,
			value = value,
			crate = crate,
		)]
		#[derive(PartialEq)]
		struct Spd {
			value: Iso<glam::DVec2, glam::IVec2>
		}
		
		let a_pos = Pos {
			value: glam::IVec2::new(10, 10).into(),
			spd: Spd {
				value: glam::IVec2::new(6, 4).into()
			}
		}.to_flux_value(Time::default());
		
		let b_pos = Pos {
			value: glam::IVec2::new(14, 18).into(),
			spd: Spd {
				value: glam::IVec2::new(4, 0).into()
			}
		}.to_flux_value(Time::default());
		
		// println!("{:?}", a_pos.poly() - b_pos.poly());
		
		assert_eq!(
			Vec::from_iter(a_pos.when_eq(&b_pos)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[(2*SEC - 83333333*NANOSEC, 2*SEC + 83333333*NANOSEC)]
		);
	}
	
	#[test]
	fn precise() {
		// https://www.desmos.com/calculator/1z97cqlopx
		
		 // Basic Check:
		let a = Poly::new(Sum::new(-193.99999999999997, [4.481238217799146e-6, -500.]), SEC);
		let b = Poly::new(Constant::from(-194.), SEC);
		assert_eq!(
			Vec::from_iter(crate::pred::WhenEq::when_eq(a, b)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(SEC-5*NANOSEC, SEC-3*NANOSEC),
				(SEC+12*NANOSEC, SEC+14*NANOSEC)
			]
		);
		
		 // Distance Check:
		let a = Poly::new([
			Sum::new(0.0036784761334161292, [1.1687626970174242e-7, 0.]),
			Sum::new(-182.00000057575835, [-7.537214753878195e-7, -500.])
		], SEC);
		let b = Poly::new([
			Constant::from(-3.8808053943969956e-5),
			Constant::from(-193.99999999999997)
		], SEC);
		let dis = Poly::new(Constant::from(12.), SEC);
		assert_eq!(
			Vec::from_iter(crate::pred::WhenDisEq::when_dis_eq(a, b, dis)
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