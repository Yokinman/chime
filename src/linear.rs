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
	'static + Copy + Clone + Debug
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

/// Underlying [`Linear`] type paired with an interfacing [`LinearIso`] type.
pub trait LinearPlus: Copy + Debug + 'static {
	type Inner: Linear;
	type Outer: LinearIso<Self::Inner>;
	fn from_inner(inner: Self::Inner) -> Self;
	fn into_inner(self) -> Self::Inner;
	fn zero() -> Self {
		Self::from_inner(Linear::zero())
	}
}

impl<T> LinearPlus for T
where
	T: Linear,
{
	type Inner = T;
	type Outer = T;
	fn from_inner(inner: Self::Inner) -> Self {
		inner
	}
	fn into_inner(self) -> Self::Inner {
		self
	}
}

impl<A, B> LinearPlus for Iso<A, B>
where
	A: Linear,
	B: LinearIso<A>,
{
	type Inner = A;
	type Outer = B;
	fn from_inner(inner: Self::Inner) -> Self {
		Iso(Some(inner), LinearIso::<A>::from_linear(inner))
	}
	fn into_inner(self) -> Self::Inner {
		let Iso(inner, outer) = self;
		inner.unwrap_or_else(|| LinearIso::<A>::into_linear(outer))
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<T: Linear>: Sized + Copy + Debug + 'static {
	fn into_linear(value: Self) -> T;
	fn from_linear(value: T) -> Self;
	// fn identity(value: Self) -> Self {
	// 	Self::inv_map(Self::map(value))
	// }
	fn linear_id(value: T) -> T {
		Self::into_linear(Self::from_linear(value))
	}
}

mod _linear_iso_impls {
	use super::{Linear, LinearIso};
	
	impl<T: Linear> LinearIso<T> for T {
		fn into_linear(value: T) -> T {
			value
		}
		fn from_linear(value: T) -> T {
			value
		}
		fn linear_id(value: T) -> T {
			value
		}
	}
	
	impl LinearIso<f64> for f32 {
		fn into_linear(value: f32) -> f64 {
			value as f64
		}
		fn from_linear(value: f64) -> f32 {
			value as f32
		}
	}
	
	impl LinearIso<f32> for f64 {
		fn into_linear(value: f64) -> f32 {
			value as f32
		}
		fn from_linear(value: f32) -> f64 {
			value as f64
		}
	}
	
	macro_rules! impl_iso_for_int {
		($b:ty: $($a:ty),+) => {$(
			impl LinearIso<$b> for $a {
				fn into_linear(value: $a) -> $b {
					value as $b
				}
				fn from_linear(value: $b) -> $a {
					value.round() as $a
				}
			}
		)+}
	}
	impl_iso_for_int!(f32: u8, u16, u32, u64, usize, i8, i16, i32, i64, isize);
	impl_iso_for_int!(f64: u8, u16, u32, u64, usize, i8, i16, i32, i64, isize);
}

/// ...
pub trait Vector<const SIZE: usize> {
	type Output;
	fn index(&self, index: usize) -> Self::Output;
}

impl<A, B, const SIZE: usize> Vector<SIZE> for Iso<A, B>
where
	A: Vector<SIZE, Output: Linear>,
	B: Vector<SIZE, Output: LinearIso<A::Output>>,
{
	type Output = Iso<A::Output, B::Output>;
	fn index(&self, index: usize) -> Self::Output {
		let Iso(inner, outer) = self;
		let inner = inner.as_ref()
			.map(|x| x.index(index))
			.unwrap_or_else(|| LinearIso::into_linear(outer.index(index)));
		Iso::from_inner(inner)
	}
}

#[cfg(feature = "glam")]
mod glam_stuff {
	use glam::*;
	use crate::linear::*;
	
	macro_rules! impl_vector {
		($vec:ty, $size:literal, $value:ty) => {
			impl Vector<$size> for $vec {
				type Output = $value;
				fn index(&self, index: usize) -> Self::Output {
					if index >= $size {
						panic!("index out of bounds")
					}
					self[index]
				}
			}
		};
	}
	impl_vector!(Vec2, 2, f32);
	impl_vector!(Vec3, 3, f32);
	impl_vector!(Vec4, 4, f32);
	impl_vector!(DVec2, 2, f64);
	impl_vector!(DVec3, 3, f64);
	impl_vector!(DVec4, 4, f64);
	impl_vector!(IVec2, 2, i32);
	impl_vector!(IVec3, 3, i32);
	impl_vector!(IVec4, 4, i32);
	impl_vector!(UVec2, 2, u32);
	impl_vector!(UVec3, 3, u32);
	impl_vector!(UVec4, 4, u32);
	
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
		($b:ty, $b_as_a:ident : $a:ty, $a_as_b:ident : $size:literal, $a_value:ty) => {
			impl LinearIso<$b> for $a {
				fn into_linear(value: $a) -> $b {
					value.$a_as_b()
				}
				fn from_linear(value: $b) -> $a {
					value.round().$b_as_a()
				}
			}
		};
	}
	impl_iso_for_vec!(Vec2, as_uvec2  : UVec2, as_vec2  : 2, u32);
	impl_iso_for_vec!(Vec3, as_uvec3  : UVec3, as_vec3  : 3, u32);
	impl_iso_for_vec!(Vec4, as_uvec4  : UVec4, as_vec4  : 4, u32);
	impl_iso_for_vec!(Vec2, as_ivec2  : IVec2, as_vec2  : 2, i32);
	impl_iso_for_vec!(Vec3, as_ivec3  : IVec3, as_vec3  : 3, i32);
	impl_iso_for_vec!(Vec4, as_ivec4  : IVec4, as_vec4  : 4, i32);
	impl_iso_for_vec!(Vec2, as_dvec2  : DVec2, as_vec2  : 2, f64);
	impl_iso_for_vec!(Vec3, as_dvec3  : DVec3, as_vec3  : 3, f64);
	impl_iso_for_vec!(Vec4, as_dvec4  : DVec4, as_vec4  : 4, f64);
	impl_iso_for_vec!(DVec2, as_uvec2 : UVec2, as_dvec2 : 2, u32);
	impl_iso_for_vec!(DVec3, as_uvec3 : UVec3, as_dvec3 : 3, u32);
	impl_iso_for_vec!(DVec4, as_uvec4 : UVec4, as_dvec4 : 4, u32);
	impl_iso_for_vec!(DVec2, as_ivec2 : IVec2, as_dvec2 : 2, i32);
	impl_iso_for_vec!(DVec3, as_ivec3 : IVec3, as_dvec3 : 3, i32);
	impl_iso_for_vec!(DVec4, as_ivec4 : IVec4, as_dvec4 : 4, i32);
	impl_iso_for_vec!(DVec2, as_vec2  : Vec2, as_dvec2  : 2, f32);
	impl_iso_for_vec!(DVec3, as_vec3  : Vec3, as_dvec3  : 3, f32);
	impl_iso_for_vec!(DVec4, as_vec4  : Vec4, as_dvec4  : 4, f32);
}

/// ...
#[derive(Copy, Clone, Debug, Default, PartialOrd, PartialEq)]
pub struct Iso<A, B>(Option<A>, B);

mod _iso_impls {
	use std::ops::{Deref, DerefMut};
	use super::Iso;
	
	impl<A, B> From<B> for Iso<A, B> {
		fn from(value: B) -> Self {
			Iso(None, value)
		}
	}
	
	impl<A, B> Iso<A, B> {
		pub fn new(outer: B) -> Self {
			Iso(None, outer)
		}
	}
	
	impl<A, B> Deref for Iso<A, B> {
		type Target = B;
		fn deref(&self) -> &Self::Target {
			let Iso(_, outer) = self;
			outer
		}
	}
	
	impl<A, B> DerefMut for Iso<A, B> {
		fn deref_mut(&mut self) -> &mut Self::Target {
			let Iso(inner, outer) = self;
			*inner = None;
			outer
		}
	}
}