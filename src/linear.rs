//! Utilities for working with vector-like values.

use std::fmt::Debug;
use std::ops::{Add, Mul};

use impl_op::impl_op;

/// A scalar value, used for multiplication with any [`Linear`] value.
#[derive(Copy, Clone, Debug)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) => Scalar(a.0 * b.0),
	(f64, Scalar)['commut] => a * b.0,
	(f32, Scalar)['commut] => a * (b.0 as f32),
	// (u64, Scalar)['commut] |
	// (u32, Scalar)['commut] |
	// (u16, Scalar)['commut] |
	// (u8,  Scalar)['commut] |
	// (i64, Scalar)['commut] |
	// (i32, Scalar)['commut] |
	// (i16, Scalar)['commut] |
	// (i8,  Scalar)['commut] => ((a as f64) * b.0).round() as Self,
}}

/// Any vector type that has addition and [`Scalar`] multiplication.
pub trait Linear:
	'static + Copy + Clone + Debug
	+ Add<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

impl<T> Linear for T
where T:
	'static + Copy + Clone + Debug
	+ Add<Output=T>
	+ Mul<Scalar, Output=T>
	+ Default
{
	fn zero() -> Self {
		Self::default()
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<T: InvLinearIso<Self> + ?Sized>: Linear {
	fn map(self) -> T;
}

pub trait InvLinearIso<T: LinearIso<Self>> {
	fn inv_map(self) -> T;
}

impl<T: Linear> LinearIso<T> for T {
	fn map(self) -> T {
		self
	}
}

impl<T: Linear> InvLinearIso<T> for T {
	fn inv_map(self) -> T {
		self
	}
}

macro_rules! impl_iso_for_int {
	($a:ty: $($b:ty),+) => {$(
		impl LinearIso<$b> for $a {
			fn map(self) -> $b {
				self.round() as $b
			}
		}
		impl InvLinearIso<$a> for $b {
			fn inv_map(self) -> $a {
				self as $a
			}
		}
	)+}
}
impl_iso_for_int!(f32: u8, u16, u32, u64, i8, i16, i32, i64);
impl_iso_for_int!(f64: u8, u16, u32, u64, i8, i16, i32, i64);