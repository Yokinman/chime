//! Compile-time assured degree of change over time.

use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::ops::{Add, Index, IndexMut, Mul, Shl, Shr};
use crate::change::{LinearValue, Scalar};
use crate::flux::{FluxAccum, SumAccum};

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

/// Linear change type.
/// 
/// - `Deg<0>`: `a`
/// - `Deg<1>`: `a + bx`
/// - `Deg<2>`: `a + bx + cx^2`
/// - etc.
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Deg<T: LinearValue, const DEG: usize>(pub T);

impl<T: LinearValue, const DEG: usize> From<T> for Deg<T, DEG> {
	fn from(value: T) -> Self {
		Self(value)
	}
}

impl<T: LinearValue, const DEG: usize> Default for Deg<T, DEG> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: LinearValue, const DEG: usize> Mul<Scalar> for Deg<T, DEG> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self::Output {
		Self(self.0 * rhs)
	}
}

impl<T: LinearValue, const DEG: usize> FluxKind for Deg<T, DEG> {
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
		impl<T: LinearValue> Shr<DegShift> for Deg<T, { $($num +)+ 0 - 1 }> {
			type Output = Deg<T, { $($num +)+ 0 }>;
			fn shr(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: LinearValue> Shl<DegShift> for Deg<T, { $($num +)+ 0 }> {
			type Output = Deg<T, { $($num +)+ 0 - 1 }>;
			fn shl(self, _: DegShift) -> Self::Output {
				self.0.into()
			}
		}
		impl<T: LinearValue> Add for Deg<T, { $($num +)+ 0 }> {
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
		impl<T: LinearValue> Add<Deg<T, $a>> for Deg<T, { $($num +)+ 0 }> {
			type Output = Deg<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Deg<T, $a>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl<T: LinearValue> Add<Deg<T, { $($num +)+ 0 }>> for Deg<T, $a> {
			type Output = Deg<T, { $($num +)+ 0 }>;
			fn add(self, rhs: Deg<T, { $($num +)+ 0 }>) -> Self::Output {
				Self::Output::from(self.0 + rhs.0)
			}
		}
		impl_deg_add!($a, 1 $($num)+);
	};
}

impl_deg_order!(1);
impl<T: LinearValue> Add for Deg<T, 0> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}
impl_deg_add!(0, 1);