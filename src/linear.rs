//! Utilities for working with vector-like values.

use std::fmt::Debug;
use std::ops::{Add, Mul};

use impl_op::impl_op;

/// A scalar value, used for multiplication with any [`Linear`] value.
#[derive(Copy, Clone, Debug)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) => Scalar(a.0 * b.0),
	(f64, Scalar) => a * b.0,
	(f32, Scalar) => ((a as f64) * b.0) as f32,
	// (u64, Scalar) |
	// (u32, Scalar) |
	// (u16, Scalar) |
	// (u8,  Scalar) |
	// (i64, Scalar) |
	// (i32, Scalar) |
	// (i16, Scalar) |
	// (i8,  Scalar) => ((a as f64) * b.0).round() as Self,
}}

#[cfg(feature = "glam")]
impl_op!{ a * b {
	(glam::Vec2,  Scalar) => (a.as_dvec2() * b).as_vec2(),
	(glam::Vec3,  Scalar) => (a.as_dvec3() * b).as_vec3(),
	(glam::Vec4,  Scalar) => (a.as_dvec4() * b).as_vec4(),
	(glam::DVec2, Scalar) |
	(glam::DVec3, Scalar) |
	(glam::DVec4, Scalar) => a * b.0,
	// (glam::IVec2,   Scalar) |
	// (glam::IVec3,   Scalar) |
	// (glam::IVec4,   Scalar) |
	// (glam::UVec2,   Scalar) |
	// (glam::UVec3,   Scalar) |
	// (glam::UVec4,   Scalar) |
	// (glam::I64Vec2, Scalar) |
	// (glam::I64Vec3, Scalar) |
	// (glam::I64Vec4, Scalar) |
	// (glam::U64Vec2, Scalar) |
	// (glam::U64Vec3, Scalar) |
	// (glam::U64Vec4, Scalar) => a * b.0,
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
				// !!! ^ Unsure if this should round or floor.
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

#[cfg(feature = "glam")]
mod glam_stuff {
	use super::*;
	use glam::*;
	
	macro_rules! impl_iso_for_vec {
		($a:ty, $as_b:ident: $b:ty, $as_a:ident) => {
			impl LinearIso<$b> for $a {
				fn map(self) -> $b {
					self.round().$as_b()
					// !!! ^ Unsure if this should round or floor.
				}
			}
			impl InvLinearIso<$a> for $b {
				fn inv_map(self) -> $a {
					self.$as_a()
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
	impl_iso_for_vec!(DVec2, as_uvec2: UVec2, as_dvec2);
	impl_iso_for_vec!(DVec3, as_uvec3: UVec3, as_dvec3);
	impl_iso_for_vec!(DVec4, as_uvec4: UVec4, as_dvec4);
	impl_iso_for_vec!(DVec2, as_ivec2: IVec2, as_dvec2);
	impl_iso_for_vec!(DVec3, as_ivec3: IVec3, as_dvec3);
	impl_iso_for_vec!(DVec4, as_ivec4: IVec4, as_dvec4);
}

#[cfg(feature = "glam")]
pub use glam_stuff::*;