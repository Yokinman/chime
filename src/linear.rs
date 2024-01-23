//! Utilities for working with vector-like values.

use std::fmt::Debug;
use std::ops::{Add, Mul, Sub};

use impl_op::impl_op;

/// A scalar value, used for multiplication with any [`Linear`] value.
#[derive(Copy, Clone, Debug, PartialEq)]
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
	'static + Copy + Clone + Debug + Send + Sync
	+ Add<Output=Self>
	+ Sub<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	fn sqrt(self) -> Self;
	
	fn sign(self) -> Self;
	
	fn next_up(self) -> Self;
	
	fn next_down(self) -> Self;
	
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
	
	fn next_up(self) -> Self {
		//! https://doc.rust-lang.org/stable/std/primitive.f64.html#method.next_up
		const TINY_BITS: u64 = 0x1; // Smallest positive f64.
		const CLEAR_SIGN_MASK: u64 = 0x7fff_ffff_ffff_ffff;
		
		let bits = self.to_bits();
		if self.is_nan() || bits == Self::INFINITY.to_bits() {
		    return self;
		}
		
		let abs = bits & CLEAR_SIGN_MASK;
		let next_bits = if abs == 0 {
		    TINY_BITS
		} else if bits == abs {
		    bits + 1
		} else {
		    bits - 1
		};
		Self::from_bits(next_bits)
	}
	
	fn next_down(self) -> Self {
		//! https://doc.rust-lang.org/stable/std/primitive.f64.html#method.next_down
		const NEG_TINY_BITS: u64 = 0x8000_0000_0000_0001; // Smallest (in magnitude) negative f64.
		const CLEAR_SIGN_MASK: u64 = 0x7fff_ffff_ffff_ffff;
		
		let bits = self.to_bits();
		if self.is_nan() || bits == Self::NEG_INFINITY.to_bits() {
		    return self;
		}
		
		let abs = bits & CLEAR_SIGN_MASK;
		let next_bits = if abs == 0 {
		    NEG_TINY_BITS
		} else if bits == abs {
		    bits - 1
		} else {
		    bits + 1
		};
		Self::from_bits(next_bits)
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
	
	fn next_up(self) -> Self {
		// https://doc.rust-lang.org/stable/std/primitive.f32.html#method.next_up
		const TINY_BITS: u32 = 0x1; // Smallest positive f32.
		const CLEAR_SIGN_MASK: u32 = 0x7fff_ffff;
		
		let bits = self.to_bits();
		if self.is_nan() || bits == Self::INFINITY.to_bits() {
		    return self;
		}
		
		let abs = bits & CLEAR_SIGN_MASK;
		let next_bits = if abs == 0 {
		    TINY_BITS
		} else if bits == abs {
		    bits + 1
		} else {
		    bits - 1
		};
		Self::from_bits(next_bits)
	}
	
	fn next_down(self) -> Self {
		// https://doc.rust-lang.org/stable/std/primitive.f32.html#method.next_down
		const NEG_TINY_BITS: u32 = 0x8000_0001; // Smallest (in magnitude) negative f32.
		const CLEAR_SIGN_MASK: u32 = 0x7fff_ffff;
		
		let bits = self.to_bits();
		if self.is_nan() || bits == Self::NEG_INFINITY.to_bits() {
		    return self;
		}
		
		let abs = bits & CLEAR_SIGN_MASK;
		let next_bits = if abs == 0 {
		    NEG_TINY_BITS
		} else if bits == abs {
		    bits - 1
		} else {
		    bits + 1
		};
		Self::from_bits(next_bits)
	}
	
	fn zero() -> Self {
		0.
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<T: LinearIsoInv<Self> + ?Sized>: Linear {
	fn map(self) -> T;
}

/// The inverse map of a linear isomorphism.
pub trait LinearIsoInv<T: LinearIso<Self>> {
	fn inv_map(self) -> T;
}

impl<T: Linear> LinearIso<T> for T {
	fn map(self) -> T {
		self
	}
}

impl<T: Linear> LinearIsoInv<T> for T {
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
		impl LinearIsoInv<$a> for $b {
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
	use glam::*;
	use crate::linear::{LinearIso, LinearIsoInv};
	
	macro_rules! impl_iso_for_vec {
		($a:ty, $as_b:ident: $b:ty, $as_a:ident) => {
			impl LinearIso<$b> for $a {
				fn map(self) -> $b {
					self.round().$as_b()
					// !!! ^ Unsure if this should round or floor.
				}
			}
			impl LinearIsoInv<$a> for $b {
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