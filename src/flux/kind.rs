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

impl<K: FluxKind> SumAccum<'_, K>
where
	K: Add<Output=K>,
{
	pub fn add<C: FluxValue>(self, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		self.accum(Scalar(1.0), value, unit)
	}
	
	pub fn sub<C: FluxValue>(self, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		self.accum(Scalar(-1.0), value, unit)
	}
	
	fn accum<C: FluxValue>(mut self, scalar: Scalar, value: &C, unit: TimeUnit) -> Self
	where
		C::Kind: FluxKind<Linear=K::Linear> + Add<K, Output=K> + Shr<DegShift>,
		<C::Kind as Shr<DegShift>>::Output: FluxKind<Linear=K::Linear>,
		C::Kind: From<<C::Kind as FluxKind>::Linear>,
		K: Add<C::Kind, Output=K> + Add<<C::Kind as Shr<DegShift>>::Output, Output=K>,
	{
		match &mut self.0 {
			FluxAccumKind::Sum { sum, depth, time } => {
				let mut sub_sum = value.value();
				value.change(<C::Kind as FluxKind>::Accum::from_kind(FluxAccumKind::Sum {
					sum: &mut sub_sum,
					depth: *depth + 1,
					time: *time,
				}));
				let depth = *depth as f64;
				let time_scale = (time.as_nanos() as f64) / ((unit >> TimeUnit::Nanosecs) as f64);
				**sum = **sum + (sub_sum * Scalar((time_scale + depth) / (depth + 1.0)) * scalar);
			},
			FluxAccumKind::Poly { poly, depth } => {
				let mut sub_poly = Poly::default();
				sub_poly.0 = value.value();
				value.change(<C::Kind as FluxKind>::Accum::from_kind(FluxAccumKind::Poly {
					poly: &mut sub_poly,
					depth: *depth + 1,
				}));
				let depth = *depth as f64;
				let unit_scale = ((unit >> TimeUnit::Nanosecs) as f64).recip();
				**poly = **poly
					+ (sub_poly * (Scalar(depth / (depth + 1.0)) * scalar))
					+ ((sub_poly >> DegShift) * (Scalar(unit_scale / (depth + 1.0)) * scalar));
				// https://www.desmos.com/calculator/mhlpjakz32
			}
		}
		self
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