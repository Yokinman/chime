//! Utilities for working with vector-like values.

use std::fmt::Debug;
use std::ops::{Add, Mul};

use impl_op::impl_op;

use crate::flux::{Sum, FluxKind};
use crate::poly::*;

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
	Sized + Copy + Clone + Debug
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
	Sized + Copy + Clone + Debug
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

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: Linear>(T);

impl<T: Linear> Add for Exp<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}

impl<T: Linear> Mul<Scalar> for Exp<T> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self {
		Self(self.0 * rhs)
	}
}

impl<T: Linear> Default for Exp<T> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: Linear> From<T> for Exp<T> {
	fn from(value: T) -> Exp<T> {
		Exp(value)
	}
}

impl<T: Linear, const DEG: usize> Roots for Sum<Exp<T>, DEG>
where
	Sum<T, DEG>: FluxKind<Linear=T> + Roots
{
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let mut b_poly = Poly::<Sum<T, DEG>>::default();
		let mut b_coeff_iter = b_poly.coeff_iter_mut();
		for &Sum(Exp(coeff)) in poly.coeff_iter() {
			*b_coeff_iter.next().unwrap() = Sum(coeff);
		}
		b_poly.0 = poly.constant().0;
		Roots::roots(b_poly)
	}
}

impl LinearIso<f64> for Exp<f64> {
	fn map(self) -> f64 {
		self.0.exp()
	}
}

impl InvLinearIso<Exp<f64>> for f64 {
	fn inv_map(self) -> Exp<f64> {
		Exp(self.ln())
	}
}

impl LinearIso<u64> for Exp<f64> {
	fn map(self) -> u64 {
		self.0.exp().round() as u64
	}
}

impl InvLinearIso<Exp<f64>> for u64 {
	fn inv_map(self) -> Exp<f64> {
		Exp((self as f64).ln())
	}
}