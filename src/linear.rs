//! Utilities for working with vector-like values.

use std::fmt::Debug;
use std::ops::{Add, Mul, Sub};

/// A scalar value, used for multiplication with any [`Linear`] value.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Scalar(pub f64);

impl Mul<Scalar> for Scalar {
	type Output = Scalar;
	fn mul(self, rhs: Scalar) -> Self::Output {
		Scalar(self.0 * rhs.0)
	}
}

impl Mul<Scalar> for f64 {
	type Output = f64;
	fn mul(self, rhs: Scalar) -> Self::Output {
		self * rhs.0
	}
}

impl Mul<Scalar> for f32 {
	type Output = f32;
	fn mul(self, rhs: Scalar) -> Self::Output {
		((self as f64) * rhs.0) as f32
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::Vec2 {
	type Output = glam::Vec2;
	fn mul(self, rhs: Scalar) -> Self::Output {
		(self.as_dvec2() * rhs).as_vec2()
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::Vec3 {
	type Output = glam::Vec3;
	fn mul(self, rhs: Scalar) -> Self::Output {
		(self.as_dvec3() * rhs).as_vec3()
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::Vec4 {
	type Output = glam::Vec4;
	fn mul(self, rhs: Scalar) -> Self::Output {
		(self.as_dvec4() * rhs).as_vec4()
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::DVec2 {
	type Output = glam::DVec2;
	fn mul(self, rhs: Scalar) -> glam::DVec2 {
		self * rhs.0
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::DVec3 {
	type Output = glam::DVec3;
	fn mul(self, rhs: Scalar) -> glam::DVec3 {
		self * rhs.0
	}
}

#[cfg(feature = "glam")]
impl Mul<Scalar> for glam::DVec4 {
	type Output = glam::DVec4;
	fn mul(self, rhs: Scalar) -> glam::DVec4 {
		self * rhs.0
	}
}

/// Any vector type that has addition and [`Scalar`] multiplication.
pub trait Linear:
	'static + Copy + Clone + Debug + Send + Sync
	+ Add<Output=Self>
	+ Sub<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	fn sqrt(self) -> Self;
	
	fn sign(self) -> Self;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

impl Linear for f64 {
	fn sqrt(self) -> Self {
		self.sqrt()
	}
	
	fn sign(self) -> Self {
		// What if sign is 0
		self.signum()
	}
	
	fn zero() -> Self {
		0.
	}
}

impl Linear for f32 {
	fn sqrt(self) -> Self {
		self.sqrt()
	}
	
	fn sign(self) -> Self {
		self.signum()
	}
	
	fn zero() -> Self {
		0.
	}
}

/// Multidimensional vector type.
pub trait LinearVec<const SIZE: usize>: Linear {
	type Value: Linear;
	fn index(&self, index: usize) -> Self::Value;
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<T: LinearIsoInv<Self>>: Linear {
	fn map(value: Self) -> T;
}

/// The inverse map of a linear isomorphism.
pub trait LinearIsoInv<T: LinearIso<Self>>: Sized {
	fn inv_map(value: Self) -> T;
}

impl<T: Linear> LinearIso<T> for T {
	fn map(value: Self) -> T {
		value
	}
}
impl<T: Linear> LinearIsoInv<T> for T {
	fn inv_map(value: Self) -> T {
		value
	}
}

impl LinearIso<f64> for f32 {
	fn map(value: Self) -> f64 {
		value as f64
	}
}
impl LinearIsoInv<f32> for f64 {
	fn inv_map(value: Self) -> f32 {
		value as f32
	}
}

impl LinearIso<f32> for f64 {
	fn map(value: Self) -> f32 {
		value as f32
	}
}
impl LinearIsoInv<f64> for f32 {
	fn inv_map(value: Self) -> f64 {
		value as f64
	}
}

macro_rules! impl_iso_for_int {
	($a:ty: $($b:ty),+) => {$(
		impl LinearIso<$b> for $a {
			fn map(value: Self) -> $b {
				value.round() as $b
			}
		}
		impl LinearIsoInv<$a> for $b {
			fn inv_map(value: Self) -> $a {
				value as $a
			}
		}
	)+}
}
impl_iso_for_int!(f32: u8, u16, u32, u64, i8, i16, i32, i64);
impl_iso_for_int!(f64: u8, u16, u32, u64, i8, i16, i32, i64);

#[cfg(feature = "glam")]
mod glam_stuff {
	use glam::*;
	use crate::linear::*;
	
	macro_rules! impl_linear_for_vec {
		($a_vec:tt $(, $b_vec:tt)+) => {
			impl_linear_for_vec!($a_vec);
			$(impl_linear_for_vec!($b_vec);)+
		};
		(($vec:ty, $size:literal, $value:ty)) => {
			impl Linear for $vec {
				fn sqrt(self) -> Self {
					self.powf(0.5)
				}
				fn sign(self) -> Self {
					self.signum()
				}
				fn zero() -> Self {
					Self::ZERO
				}
			}
			impl LinearVec<$size> for $vec {
				type Value = $value;
				fn index(&self, index: usize) -> Self::Value {
					if index >= $size {
						panic!("index out of bounds")
					}
					self[index]
				}
			}
		};
	}
	impl_linear_for_vec!(
		(Vec2, 2, f32),
		(Vec3, 3, f32),
		(Vec4, 4, f32),
		(DVec2, 2, f64),
		(DVec3, 3, f64),
		(DVec4, 4, f64)
	);
	
	macro_rules! impl_iso_for_vec {
		($a:ty, $as_b:ident: $b:ty, $as_a:ident) => {
			impl LinearIso<$b> for $a {
				fn map(value: Self) -> $b {
					value.round().$as_b()
				}
			}
			impl LinearIsoInv<$a> for $b {
				fn inv_map(value: Self) -> $a {
					value.$as_a()
				}
			}
		}
	}
	impl_iso_for_vec!(Vec2, as_uvec2: UVec2, as_vec2);
	impl_iso_for_vec!(Vec3, as_uvec3: UVec3, as_vec3);
	impl_iso_for_vec!(Vec4, as_uvec4: UVec4, as_vec4);
	impl_iso_for_vec!(Vec2, as_ivec2: IVec2, as_vec2);
	impl_iso_for_vec!(Vec3, as_ivec3: IVec3, as_vec3);
	impl_iso_for_vec!(Vec4, as_ivec4: IVec4, as_vec4);
	impl_iso_for_vec!(Vec2, as_dvec2: DVec2, as_vec2);
	impl_iso_for_vec!(Vec3, as_dvec3: DVec3, as_vec3);
	impl_iso_for_vec!(Vec4, as_dvec4: DVec4, as_vec4);
	impl_iso_for_vec!(DVec2, as_uvec2: UVec2, as_dvec2);
	impl_iso_for_vec!(DVec3, as_uvec3: UVec3, as_dvec3);
	impl_iso_for_vec!(DVec4, as_uvec4: UVec4, as_dvec4);
	impl_iso_for_vec!(DVec2, as_ivec2: IVec2, as_dvec2);
	impl_iso_for_vec!(DVec3, as_ivec3: IVec3, as_dvec3);
	impl_iso_for_vec!(DVec4, as_ivec4: IVec4, as_dvec4);
	impl_iso_for_vec!(DVec2, as_vec2: Vec2, as_dvec2);
	impl_iso_for_vec!(DVec3, as_vec3: Vec3, as_dvec3);
	impl_iso_for_vec!(DVec4, as_vec4: Vec4, as_dvec4);
}

#[cfg(feature = "glam")]
pub use glam_stuff::*;