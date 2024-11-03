//! Utilities for working with vector-like values.

use std::fmt::Debug;

/// A scalar value, used for multiplication with any [`Linear`] value.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Scalar(f64);

mod _scalar_impls {
	use std::ops::Mul;
	use super::Scalar;
	
	impl From<f64> for Scalar {
		fn from(value: f64) -> Self {
			Scalar(value)
		}
	}
	
	impl From<Scalar> for f64 {
		fn from(Scalar(value): Scalar) -> Self {
			value
		}
	}
	
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
}

/// Any vector type that has addition and [`Scalar`] multiplication.
/// 
/// This basically just represents floating-point numbers and vector types.
/// For vectors, operations are applied per component in parallel.
pub trait Linear: Basis<Inner=Self> + Clone + Debug + LinearIso<Self> + 'static {
	fn add(self, other: Self) -> Self;
	
	fn sub(self, other: Self) -> Self;
	
	fn mul_scalar(self, scalar: Scalar) -> Self;
	
	fn sqr(self) -> Self;
	
	fn sqrt(self) -> Self;
	
	fn sign(&self) -> Self;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&<Self as Linear>::zero())
	}
}

mod _linear_impls {
	use super::{Linear, Scalar};
	
	impl Linear for f64 {
		fn add(self, other: Self) -> Self {
			self + other
		}
		fn sub(self, other: Self) -> Self {
			self - other
		}
		fn mul_scalar(self, scalar: Scalar) -> Self {
			self * scalar
		}
		fn sqr(self) -> Self {
			self * self
		}
		fn sqrt(self) -> Self {
			self.sqrt()
		}
		fn sign(&self) -> Self {
			// What if sign is 0
			self.signum()
		}
		fn zero() -> Self {
			0.
		}
	}
	
	impl Linear for f32 {
		fn add(self, other: Self) -> Self {
			self + other
		}
		fn sub(self, other: Self) -> Self {
			self - other
		}
		fn mul_scalar(self, scalar: Scalar) -> Self {
			self * scalar
		}
		fn sqr(self) -> Self {
			self * self
		}
		fn sqrt(self) -> Self {
			self.sqrt()
		}
		fn sign(&self) -> Self {
			self.signum()
		}
		fn zero() -> Self {
			0.
		}
	}
	
	impl<T: Linear, const SIZE: usize> Linear for [T; SIZE] {
		fn add(self, other: Self) -> Self {
			let mut iter = other.into_iter();
			self.map(|x| unsafe {
				// SAFETY: `self` is the same length as `other`.
				x.add(iter.next().unwrap_unchecked())
			})
		}
		fn sub(self, other: Self) -> Self {
			let mut iter = other.into_iter();
			self.map(|x| unsafe {
				// SAFETY: `self` is the same length as `other`.
				x.sub(iter.next().unwrap_unchecked())
			})
		}
		fn mul_scalar(self, scalar: Scalar) -> Self {
			self.map(|x| x.mul_scalar(scalar))
		}
		fn sqr(self) -> Self {
			self.map(T::sqr)
		}
		fn sqrt(self) -> Self {
			self.map(T::sqrt)
		}
		fn sign(&self) -> Self {
			self.each_ref().map(T::sign)
		}
		fn zero() -> Self {
			std::array::from_fn(|_| Linear::zero())
		}
	}
}

/// A [`Linear`] type packaged with extra information (e.g. [`Iso`]).
pub trait Basis: Clone + Debug + 'static {
	type Inner: Linear;
	
	fn from_inner(inner: Self::Inner) -> Self;
	fn into_inner(self) -> Self::Inner;
	
	fn with<R>(&self, f: impl FnOnce(&Self::Inner) -> R) -> R;
	
	fn inner_id(inner: Self::Inner) -> Self::Inner;
	
	fn inner_eq(&self, other: &Self) -> bool
	where
		Self::Inner: PartialEq	
	{
		self.with(|x| other.with(|y| x == y))
	}
	
	fn zero() -> Self {
		Self::from_inner(Linear::zero())
	}
}

mod _linear_plus_impls {
	use super::{Iso, Linear, LinearIso, Basis, Simple};
	
	impl<T> Basis for T
	where
		T: Linear + Simple,
	{
		type Inner = Self;
		fn from_inner(inner: Self::Inner) -> Self {
			inner
		}
		fn into_inner(self) -> Self::Inner {
			self
		}
		fn with<R>(&self, f: impl FnOnce(&Self::Inner) -> R) -> R {
			f(&self)
		}
		fn inner_id(inner: Self::Inner) -> Self::Inner {
			inner
		}
	}
	
	impl<T, const N: usize> Basis for [T; N]
	where
		T: Basis,
	{
		type Inner = [T::Inner; N];
		fn from_inner(inner: Self::Inner) -> Self {
			inner.map(T::from_inner)
		}
		fn into_inner(self) -> Self::Inner {
			self.map(T::into_inner)
		}
		fn with<R>(&self, f: impl FnOnce(&Self::Inner) -> R) -> R {
			f(&self.clone().into_inner())
		}
		fn inner_id(inner: Self::Inner) -> Self::Inner {
			inner.map(T::inner_id)
		}
	}
	
	impl<A, B> Basis for Iso<A, B>
	where
		A: Linear,
		B: LinearIso<A>,
	{
		type Inner = A;
		fn from_inner(inner: Self::Inner) -> Self {
			Iso(Some(inner.clone()), LinearIso::<A>::from_linear(inner))
		}
		fn into_inner(self) -> Self::Inner {
			let Iso(inner, outer) = self;
			inner.unwrap_or_else(|| LinearIso::<A>::into_linear(outer))
		}
		fn with<R>(&self, f: impl FnOnce(&Self::Inner) -> R) -> R {
			let Iso(inner, outer) = self;
			if let Some(inner) = inner.as_ref() {
				f(inner)
			} else {
				f(&LinearIso::<A>::into_linear(outer.clone()))
			}
		}
		fn inner_id(inner: Self::Inner) -> Self::Inner {
			B::linear_id(inner)
		}
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<T: Linear>: Sized + Clone + Debug + 'static {
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
	
	impl LinearIso<f64> for f64 {
		fn into_linear(value: Self) -> f64 {
			value
		}
		fn from_linear(value: f64) -> Self {
			value
		}
	}
	
	impl LinearIso<f32> for f32 {
		fn into_linear(value: Self) -> f32 {
			value
		}
		fn from_linear(value: f32) -> Self {
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
	
	impl<A, B, const SIZE: usize> LinearIso<[A; SIZE]> for [B; SIZE]
	where
		A: Linear,
		B: LinearIso<A>,
	{
		fn into_linear(value: Self) -> [A; SIZE] {
			value.map(B::into_linear)
		}
		fn from_linear(value: [A; SIZE]) -> Self {
			value.map(B::from_linear)
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
	
	// !!! The integer isomorphisms should not round by default. Instead, they
	// should do the default `as` cast behavior and then interface types like
	// `Rounded<T>` can be used for specifying that behavior. However, before
	// that can be done: I don't think flooring is handled by prediction yet.
}

/// Multidimensional values. Essentially [`std::ops::Index<usize>`] but with a
/// definable size and a by-value output.
pub trait Vector<const SIZE: usize> {
	type Output;
	fn index(&self, index: usize) -> Self::Output;
	// !!! I want to add a lifetime parameter to `Self::Output` so that I can
	// return values by reference when desired. However, this will also require
	// refactoring systems to account for bounds like: `Temporal<&T>: Flux`,
	// `Sum<&T, D>: FluxKind`, etc.
}

impl<T, const SIZE: usize> Vector<SIZE> for [T; SIZE]
where
	T: Clone,
{
	type Output = T;
	fn index(&self, index: usize) -> Self::Output {
		self[index].clone()
	}
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
				fn add(self, other: Self) -> Self {
					self + other
				}
				fn sub(self, other: Self) -> Self {
					self - other
				}
				fn mul_scalar(self, scalar: Scalar) -> Self {
					self * scalar
				}
				fn sqr(self) -> Self {
					self.powf(2.)
				}
				fn sqrt(self) -> Self {
					self.powf(0.5)
				}
				fn sign(&self) -> Self {
					self.signum()
				}
				fn zero() -> Self {
					Self::ZERO
				}
			}
			impl LinearIso<$vec> for $vec {
				fn into_linear(value: Self) -> $vec {
			        value
			    }
				fn from_linear(value: $vec) -> Self {
			        value
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
#[derive(Copy, Clone, Debug, Default)]
pub struct Iso<A, B>(Option<A>, B);

mod _iso_impls {
	use std::cmp::Ordering;
	use std::ops::{Deref, DerefMut};
	use crate::kind::constant::Constant;
	use super::{Iso, Linear, LinearIso, Scalar};
	
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
	
	impl<A, B, X, Y> PartialEq<Iso<X, Y>> for Iso<A, B>
	where
		B: PartialEq<Y>,
	{
		fn eq(&self, other: &Iso<X, Y>) -> bool {
			self.1 == other.1
		}
	}
	
	impl<A, B, X, Y> PartialOrd<Iso<X, Y>> for Iso<A, B>
	where
		B: PartialOrd<Y>,
	{
		fn partial_cmp(&self, other: &Iso<X, Y>) -> Option<Ordering> {
			self.1.partial_cmp(&other.1)
		}
	}
	
	impl<A, B> crate::Flux for Iso<A, B>
	where
		A: Linear,
		B: LinearIso<A>,
	{
		type Basis = Self;
		type Kind = Constant<Self>;
		fn basis(&self) -> Self::Basis {
			self.clone()
		}
		fn change(&self, basis: Self::Basis) -> Self::Kind {
			basis.into()
		}
	}
	
	impl<A, B> crate::ToMoment for Iso<A, B>
	where
		A: Linear,
		B: LinearIso<A>,
	{
		type Moment<'a> = Self;
		fn to_moment(&self, _time: Scalar) -> Self::Moment<'_> {
			self.clone()
		}
	}
}

/// Types that represent a specific concept.
/// 
/// This trait is generally implemented for non-generic types, or generic types
/// that are bounded by specific traits.
/// 
/// This generally excludes types with unbounded generic parameters, such as:
/// `Option<T>`, `[T; N]`, `(T,U,..)`, etc.
/// 
/// Used to support blanket impls with exceptions for generic types, allowing
/// for special case impls.
pub trait Simple {}