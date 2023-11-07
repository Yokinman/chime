//! Summation over time.

use std::ops::{Add, Mul, Sub, Shl, Shr};
use crate::{*, kind::*, linear::*, exp::*};

/// Summation over time.
/// 
/// - `Sum<0>`: `a`
/// - `Sum<1>`: `a + bx`
/// - `Sum<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Sum<T: Linear, const DEG: usize>(pub T);

impl<T: Linear, const DEG: usize> From<T> for Sum<T, DEG> {
	fn from(value: T) -> Self {
		Self(value)
	}
}

impl<T: Linear, const DEG: usize> Mul<Scalar> for Sum<T, DEG> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self::Output {
		Self(self.0 * rhs)
	}
}

impl<T: Linear, const DEG: usize> FluxKind for Sum<T, DEG> {
	type Value = T;
	type Coeffs = [Self; DEG];
	type Accum<'a> = SumAccum<'a, Self>;
	
	const DEGREE: usize = DEG;
	
	fn zero() -> Self {
		Self(T::zero())
	}
	
	fn zero_coeffs() -> Self::Coeffs {
		[Self::zero(); DEG]
	}
}

/// Degree sequential ordering.
macro_rules! impl_deg_order {
	(1  1  $($num:tt)*) => { impl_deg_order!(2  $($num)*); };
	(2  2  $($num:tt)*) => { impl_deg_order!(4  $($num)*); };
	(4  4  $($num:tt)*) => { impl_deg_order!(8  $($num)*); };
	(8  8  $($num:tt)*) => { impl_deg_order!(16 $($num)*); };
	// (16 16 $($num:tt)*) => { impl_deg_order!(32 $($num)*); };
	// (32 32 $($num:tt)*) => { impl_deg_order!(64 $($num)*); };
	(16) => {/* break */};
	($($num:tt)+) => {
		impl<T: Linear> Shr<DegShift> for Sum<T, { $($num +)+ 0 - 1 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn shr(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: Linear> Shl<DegShift> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 - 1 }>;
			fn shl(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: Linear> Add for Sum<T, { $($num +)+ 0 }> {
			type Output = Self;
			fn add(self, rhs: Self) -> Self {
				Self(self.0 + rhs.0)
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
	($a:tt, 8  8  $($num:tt)*) => { impl_deg_add!($a, 16 $($num)*); };
	// ($a:tt, 16 16 $($num:tt)*) => { impl_deg_add!($a, 32 $($num)*); };
	// ($a:tt, 32 32 $($num:tt)*) => { impl_deg_add!($a, 64 $($num)*); };
	($a:tt, 16) => {/* break */};
	($a:tt, $($num:tt)+) => {
		impl<T: Linear> Add<Sum<T, $a>> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, $a>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl<T: Linear> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, $a> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);
impl<T: Linear> Add for Sum<T, 0> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}
impl_deg_add!(0, 1);

const LIMIT: f64 = 0.01;

impl Roots for Sum<f64, 0> {
	fn roots(_: Poly<Self>) -> Result<RootList, RootList> {
		Ok([].into())
	}
}

impl Roots for Sum<f64, 1> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a, [Sum(b)]) = poly;
		if b == 0. {
			return Sum::<f64, 0>::roots(Poly(a, []))
		}
		Ok([-a / b].into())
	}
}

impl Roots for Sum<f64, 2> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a, [Sum(b), Sum(c)]) = poly;
		if c == 0. {
			return Sum::<f64, 1>::roots(Poly(a, [Sum(b)]))
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum::<f64, 1>::roots(Poly(b, [Sum(c)]))?
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
				Sum::<f64, 1>::roots(Poly(a, [Sum(b)]))?
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
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a, [Sum(b), Sum(c), Sum(d)]) = poly;
		if d == 0. {
			return Sum::<f64, 2>::roots(Poly(a, [Sum(b), Sum(c)]))
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum::<f64, 2>::roots(Poly(b, [Sum(c), Sum(d)]))?
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
				Some(std::cmp::Ordering::Less) => {
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
				Some(std::cmp::Ordering::Greater) => {
					let n = discriminant.sqrt();
					debug_assert!(!n.is_nan());
					Ok([(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN].into())
				},
				
				 // This shouldn't really ever run, but just in case:
				_ => {
					let roots = [
						[0.].into(),
						Sum::<f64, 2>::roots(Poly(a, [Sum(b), Sum(c)]))?
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
				Sum::<f64, 2>::roots(Poly(a, [Sum(b), Sum(c)]))?
			];
			return Ok(roots.concat().into())
		}
		
		 // General Cubic:
		let n = -c / (3.0 * d);
		let depressed_cubic = Poly::<Sum<f64, 3>>(
			a + (n * (b + (n * c * (2.0/3.0)))),
			[
				Sum(b + (n * c)),
				Sum(0.0),
				Sum(d),
			]
		);
		Ok(Roots::roots(depressed_cubic)?
			.iter()
			.map(|x| x + n)
			.collect())
	}
}

impl Roots for Sum<f64, 4> {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let Poly(a, [Sum(b), Sum(c), Sum(d), Sum(e)]) = poly;
		if e == 0. {
			return Sum::<f64, 3>::roots(Poly(a, [Sum(b), Sum(c), Sum(d)]))
		}
		if a == 0. {
			let roots = [
				[0.].into(),
				Sum::<f64, 3>::roots(Poly(b, [Sum(c), Sum(d), Sum(e)]))?
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
			return Ok(Sum::<f64, 2>::roots(Poly(a, [Sum(c), Sum(e)]))?
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
			let resolvent_cubic = Poly::<Sum<f64, 3>>(
				-(q * q) / 2.0,
				[
					Sum((r * r) - p),
					Sum(2.0 * r),
					Sum(1.0),
				]
			);
			let m = *resolvent_cubic.real_roots().unwrap_or_default() 
				.iter().find(|&r| *r >= 0.0)
				.expect("this shouldn't happen, probably a precision issue");
			let sqrt_2m = (2.0 * m).sqrt();
			let quad_a = Poly::<Sum<f64, 2>>( (q / sqrt_2m) + r + m, [Sum(-sqrt_2m), Sum(1.0)]);
			let quad_b = Poly::<Sum<f64, 2>>(-(q / sqrt_2m) + r + m, [Sum( sqrt_2m), Sum(1.0)]);
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
				Sum::<f64, 3>::roots(Poly(a, [Sum(b), Sum(c), Sum(d)]))?
			];
			return Ok(roots.concat().into())
		}
		
		 // General Quartic:
		let n = -d / (4.0 * e);
		let depressed_quartic = Poly::<Sum<f64, 4>>(
			a + (n * (b + (n * (c + (n * d * (3.0/4.0)))))),
			[
				Sum(b + (n * (c + (n * d)) * 2.0)),
				Sum(c + (n * d * (3.0/2.0))),
				Sum(0.0),
				Sum(e),
			]
		);
		Ok(Roots::roots(depressed_quartic)?
			.iter()
			.map(|x| x + n)
			.collect())
	}
}

#[cfg(feature = "glam")]
impl<const DEG: usize> Roots for Sum<glam::DVec2, DEG>
where
	Sum<f64, DEG>: FluxKind<Value=f64> + Roots
{
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let mut x_poly = Poly::<Sum<f64, DEG>> {
			0: poly.0.x,
			..Default::default()
		};
		let mut y_poly = Poly::<Sum<f64, DEG>> {
			0: poly.0.y,
			..Default::default()
		};
		for (index, coeff) in poly.coeff_iter().enumerate() {
			x_poly.1[index].0 = coeff.0.x;
			y_poly.1[index].0 = coeff.0.y;
		}
		
		type Output = fn(RootList) -> Result<RootList, RootList>;
		let (output, mut x_roots, mut y_roots): (Output, _, _) = match (
			Roots::roots(x_poly),
			Roots::roots(y_poly)
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
				std::cmp::Ordering::Less    => x = x_iter.next(),
				std::cmp::Ordering::Greater => y = y_iter.next(),
				std::cmp::Ordering::Equal   => {
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

impl<T: Linear, const DEG: usize> Roots for Sum<Exp<T>, DEG>
where
	Sum<T, DEG>: FluxKind<Value=T> + Roots
{
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let mut b_poly = Poly::<Sum<T, DEG>> {
			0: poly.constant().0,
			..Default::default()
		};
		let mut b_coeff_iter = b_poly.coeff_iter_mut();
		for &Sum(Exp(coeff)) in poly.coeff_iter() {
			*b_coeff_iter.next().unwrap() = Sum(coeff);
		}
		Roots::roots(b_poly)
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
		self.accum(Scalar(1.0), rhs)
	}
}

impl<K: FluxKind, V: Flux> Sub<Change<'_, V>> for SumAccum<'_, K>
where
	(K, V::Kind): SumAccumHelper<K, V::Kind>
{
	type Output = Self;
	fn sub(self, rhs: Change<'_, V>) -> Self {
		self.accum(Scalar(-1.0), rhs)
	}
}

/// Used to remove redundant trait bounds.
#[doc(hidden)]
pub trait SumAccumHelper<A: FluxKind, B: FluxKind> {
	fn eval<C: Flux<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		value: &C,
		unit: TimeUnit,
	);
}

impl<A, B> SumAccumHelper<A, B> for (A, B)
where
	A: FluxKind,
	B: FluxKind<Value=A::Value>,
	A: Add<B, Output=A> + Add<<B as Shr<DegShift>>::Output, Output=A>,
	B: Add<A, Output=A> + Shr<DegShift> + From<B::Value>,
	<B as Shr<DegShift>>::Output: FluxKind<Value=A::Value>,
{
	fn eval<V: Flux<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		flux: &V,
		unit: TimeUnit,
	) {
		match kind {
			FluxAccumKind::Value { value, depth, time, offset } => {
				let mut sub_value = flux.value();
				flux.change(B::Accum::from_kind(FluxAccumKind::Value {
					value: &mut sub_value,
					depth: *depth + 1,
					time: *time,
					offset: flux.time(),
				}));
				let depth = *depth as f64;
				let time_scale = (time.as_secs_f64() - offset.as_secs_f64()) / (1*unit).as_secs_f64();
				**value = **value + (sub_value * Scalar((time_scale + depth) / (depth + 1.0)) * scalar);
			}
			FluxAccumKind::Poly { poly, depth } => {
				let mut sub_poly = Poly {
					0: flux.value(),
					..Default::default()
				};
				flux.change(B::Accum::from_kind(FluxAccumKind::Poly {
					poly: &mut sub_poly,
					depth: *depth + 1,
				}));
				let depth = *depth as f64;
				let time_scale = (1*unit).as_secs_f64().recip();
				**poly = **poly
					+ (sub_poly * (Scalar(depth / (depth + 1.0)) * scalar))
					+ ((sub_poly >> DegShift) * (Scalar(time_scale / (depth + 1.0)) * scalar));
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
		let a = Poly::<Sum<f64, 4>>(1.5, [Sum(2.5), Sum(3.2), Sum(4.5), Sum(5.7)]);
		let b = Poly::<Sum<f64, 2>>(7.1, [Sum(5.9), Sum(3.1)]);
		assert_eq!(a + b, Poly(8.6, [Sum(8.4), Sum(6.300000000000001), Sum(4.5), Sum(5.7)]));
		assert_eq!(b + a, Poly(8.6, [Sum(8.4), Sum(6.300000000000001), Sum(4.5), Sum(5.7)]));
	}
	
	#[test]
	fn sub() {
		let a = Poly::<Sum<f64, 4>>(1.5, [Sum(2.5), Sum(3.2), Sum(4.5), Sum(5.7)]);
		let b = Poly::<Sum<f64, 2>>(7.1, [Sum(5.9), Sum(3.1)]);
		assert_eq!(a - b, Poly(-5.6, [Sum(-3.4000000000000004), Sum(0.10000000000000009), Sum(4.5), Sum(5.7)]));
		assert_eq!(b - a, Poly( 5.6, [Sum(3.4000000000000004), Sum(-0.10000000000000009), Sum(-4.5), Sum(-5.7)]));
	}
	
	#[test]
	fn scalar_mul() {
		let a = Poly::<Sum<f64, 4>>(1.5, [Sum(2.5), Sum(3.2), Sum(4.5), Sum(5.7)]);
		let b = Poly::<Sum<f64, 2>>(7.1, [Sum(5.9), Sum(3.1)]);
		assert_eq!(a * Scalar(1.5), Poly(2.25, [Sum(3.75), Sum(4.800000000000001), Sum(6.75), Sum(8.55)]));
		assert_eq!(b * Scalar(1.5), Poly(10.649999999999999, [Sum(8.850000000000001), Sum(4.65)]));
	}
	
	fn assert_roots<K: FluxKind>(p: Poly<K>, expected_roots: &[f64])
	where
		K: Roots
	{
		let r = p.real_roots().unwrap_or_default();
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
		assert_roots(Poly::<Sum<f64, 0>>(2.0, []), &[]);
		assert_roots(Poly::<Sum<f64, 0>>(0.0, []), &[]);
	}
	
	#[test]
	fn linear() {
		assert_roots(Poly::<Sum<f64, 1>>(20.0, [Sum(-4.0)]), &[5.0]);
		assert_roots(Poly::<Sum<f64, 1>>(20.0, [Sum(4.0)]), &[-5.0]);
		assert_roots(Poly::<Sum<f64, 1>>(-0.0, [Sum(4.0/3.0)]), &[0.0]);
	}
	
	#[test]
	fn quadratic() {
		assert_roots(
			Poly::<Sum<f64, 2>>(20.0, [Sum(4.0), Sum(7.0)]),
			&[]
		);
		assert_roots(
			Poly::<Sum<f64, 2>>(20.0, [Sum(4.0), Sum(-7.0)]),
			&[-10.0/7.0, 2.0]
		);
		assert_roots(
			Poly::<Sum<f64, 2>>( 40.0/3.0, [Sum( 2.0/3.0), Sum(-17.0/100.0)]),
			Poly::<Sum<f64, 2>>(-40.0/3.0, [Sum(-2.0/3.0), Sum( 17.0/100.0)]).real_roots().unwrap().as_ref()
		);
		assert_roots(
			Poly::<Sum<f64, 2>>(0.0, [Sum(4.0/6.0), Sum(-17.0/100.0)]),
			&[0.0, 200.0/51.0]
		);
	}
	
	#[test]
	fn cubic() {
		assert_roots(
			Poly::<Sum<f64, 3>>(6.0, [Sum(-2077.5), Sum(-17000.0/77.0), Sum(6712.0/70.0)]),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> Poly<Sum<f64, 3>> {
			Poly(a, [
				Sum(b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s)),
				Sum(c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + d/(2.0*3.0*s*s*s)),
				Sum(d/(2.0*3.0*s*s*s))
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
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> Poly<Sum<f64, 4>> {
			Poly(a, [
				Sum(b/s + c/(2.0*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Sum(c/(2.0*s*s) + d/(2.0*3.0*s*s*s) + (2.0*d)/(2.0*3.0*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*3.0*e)/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Sum(d/(2.0*3.0*s*s*s) + e/(2.0*3.0*4.0*s*s*s*s) + (2.0*e)/(2.0*3.0*4.0*s*s*s*s) + (3.0*e)/(2.0*3.0*4.0*s*s*s*s)),
				Sum(e/(2.0*3.0*4.0*s*s*s*s)),
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
			Poly::<Sum<f64, 4>>(7500.0, [Sum(0.0), Sum(-1000.0/(t*t)), Sum(0.0), Sum(27.0/(t*t*t*t))]),
			&[
				-2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
				-2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 - 6.0*f64::sqrt(19.0)),
				 2000000000000.0 * f64::sqrt(60.0 + 6.0*f64::sqrt(19.0)),
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
		// https://www.desmos.com/calculator/usraewodwz
	}
}