//! Defining a *kind* of change over time.

use super::*;

use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::ops::{Add, Index, IndexMut, Mul, Shl, Shr};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind: Copy + Clone + Default + Debug + Mul<Scalar, Output=Self> {
	const DEGREE: usize;
	
	type Linear: LinearValue;
	
	type Accum<'a>: FluxAccum<'a, Self> where Self::Linear: 'a;
	
	type Coeffs:
		Copy + Clone + Debug
		+ IntoIterator<Item=Self>
		+ AsRef<[Self]> + AsMut<[Self]>
		+ Borrow<[Self]> + BorrowMut<[Self]>
		+ Index<usize, Output=Self> + IndexMut<usize>;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
	
	fn zero_coeffs() -> Self::Coeffs;
}

/// Summation over time.
/// 
/// - `Sum<0>`: `a`
/// - `Sum<1>`: `a + bx`
/// - `Sum<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Sum<T: LinearValue, const DEG: usize>(pub T);

impl<T: LinearValue, const DEG: usize> From<T> for Sum<T, DEG> {
	fn from(value: T) -> Self {
		Self(value)
	}
}

impl<T: LinearValue, const DEG: usize> Default for Sum<T, DEG> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: LinearValue, const DEG: usize> Mul<Scalar> for Sum<T, DEG> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self::Output {
		Self(self.0 * rhs)
	}
}

impl<T: LinearValue, const DEG: usize> FluxKind for Sum<T, DEG> {
	const DEGREE: usize = DEG;
	
	type Linear = T;
	type Accum<'a> = SumAccum<'a, Self> where T: 'a;
	
	type Coeffs = [Self; DEG];
	
	fn zero() -> Self {
		Self(T::zero())
	}
	
	fn zero_coeffs() -> Self::Coeffs {
		[Self::zero(); DEG]
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
	pub fn add<V: FluxValue>(self, value: &V, unit: TimeUnit) -> Self
	where
		(K, V::Kind): SumAccumHelper<K, V::Kind>,
	{
		self.accum(Scalar(1.0), value, unit)
	}
	
	pub fn sub<V: FluxValue>(self, value: &V, unit: TimeUnit) -> Self
	where
		(K, V::Kind): SumAccumHelper<K, V::Kind>,
	{
		self.accum(Scalar(-1.0), value, unit)
	}
	
	fn accum<V: FluxValue>(mut self, scalar: Scalar, value: &V, unit: TimeUnit) -> Self
	where
		(K, V::Kind): SumAccumHelper<K, V::Kind>,
	{
		<(K, V::Kind)>::eval(&mut self.0, scalar, value, unit);
		self
	}
}

/// Used to remove redundant trait bounds.
#[doc(hidden)]
pub trait SumAccumHelper<A: FluxKind, B: FluxKind> {
	fn eval<C: FluxValue<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		value: &C,
		unit: TimeUnit,
	);
}

impl<A, B> SumAccumHelper<A, B> for (A, B)
where
	A: FluxKind,
	B: FluxKind<Linear=A::Linear>,
	A: Add<B, Output=A> + Add<<B as Shr<DegShift>>::Output, Output=A>,
	B: Add<A, Output=A> + Shr<DegShift> + From<B::Linear>,
	<B as Shr<DegShift>>::Output: FluxKind<Linear=A::Linear>,
{
	fn eval<V: FluxValue<Kind=B>>(
		kind: &mut FluxAccumKind<'_, A>,
		scalar: Scalar,
		value: &V,
		unit: TimeUnit,
	) {
		match kind {
			FluxAccumKind::Sum { sum, depth, time, offset } => {
				let mut sub_sum = value.value();
				value.change(B::Accum::from_kind(FluxAccumKind::Sum {
					sum: &mut sub_sum,
					depth: *depth + 1,
					time: *time,
					offset: value.time(),
				}));
				let depth = *depth as f64;
				let time_scale = (time.as_secs_f64() - offset.as_secs_f64()) / (1*unit).as_secs_f64();
				**sum = **sum + (sub_sum * Scalar((time_scale + depth) / (depth + 1.0)) * scalar);
			},
			FluxAccumKind::Poly { poly, depth } => {
				let mut sub_poly = Poly::default();
				sub_poly.0 = value.value();
				value.change(B::Accum::from_kind(FluxAccumKind::Poly {
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

/// Used with `Shr` and `Shl` for upgrading/downgrading a [`FluxKind`].
#[derive(Copy, Clone, Debug)]
pub struct DegShift;
// !!! assert { K::DEGREE + 1 } == <K as Shr<DegShift>>::Output::DEGREE ?

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
		impl<T: LinearValue> Shr<DegShift> for Sum<T, { $($num +)+ 0 - 1 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn shr(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: LinearValue> Shl<DegShift> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 - 1 }>;
			fn shl(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: LinearValue> Add for Sum<T, { $($num +)+ 0 }> {
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
		impl<T: LinearValue> Add<Sum<T, $a>> for Sum<T, { $($num +)+ 0 }> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, $a>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl<T: LinearValue> Add<Sum<T, { $($num +)+ 0 }>> for Sum<T, $a> {
			type Output = Sum<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Sum<T, { $($num +)+ 0 }>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}
impl_deg_order!(1);
impl<T: LinearValue> Add for Sum<T, 0> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}
impl_deg_add!(0, 1);