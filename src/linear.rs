//! Utilities for working with vector-like values.

/// Any vector type that has addition and scalar multiplication.
/// 
/// This basically just represents floating-point numbers and vector types.
/// For vectors, operations are applied per component in parallel.
pub trait Linear: Copy + Clone + PartialEq + PartialOrd {
	fn add(self, other: Self) -> Self;
	
	fn sub(self, other: Self) -> Self;
	
	fn mul(self, other: Self) -> Self;
	
	fn div(self, other: Self) -> Self;
	
	fn pow(self, other: Self) -> Self;
	
	fn exp(self) -> Self;
	
	fn ln(self) -> Self;
	
	fn sqr(self) -> Self;
	
	fn sqrt(self) -> Self;
	
	fn sign(&self) -> Self;
	
	fn from_f64(n: f64) -> Self;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&<Self as Linear>::zero())
	}
}

mod _linear_impls {
	use super::Linear;
	
	impl Linear for f64 {
		fn add(self, other: Self) -> Self {
			self + other
		}
		fn sub(self, other: Self) -> Self {
			self - other
		}
		fn mul(self, other: Self) -> Self {
			self * other
		}
		fn div(self, other: Self) -> Self {
			self / other
		}
		fn pow(self, other: Self) -> Self {
			self.powf(other)
		}
		fn exp(self) -> Self {
			f64::exp(self)
		}
		fn ln(self) -> Self {
			f64::ln(self)
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
		fn from_f64(n: f64) -> Self {
			n
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
		fn mul(self, other: Self) -> Self {
			self * other
		}
		fn div(self, other: Self) -> Self {
			self / other
		}
		fn pow(self, other: Self) -> Self {
			self.powf(other)
		}
		fn exp(self) -> Self {
			f32::exp(self)
		}
		fn ln(self) -> Self {
			f32::ln(self)
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
		fn from_f64(n: f64) -> Self {
			n as f32
		}
		fn zero() -> Self {
			0.
		}
	}
}

/// A [`Linear`] type packaged with extra information (e.g. [`Iso`]).
pub trait Basis: Clone {
	type Inner: Linear;
	
	fn from_inner(inner: Self::Inner) -> Self;
	
	fn map_inner<F>(self, f: F) -> Self
	where
		F: Fn(Self::Inner) -> Self::Inner
	{
		Self::each_map_inner([self], |[x]| f(x))
	}
	
	fn zip_map_inner<F>(self, other: Self, f: F) -> Self
	where
		F: Fn(Self::Inner, Self::Inner) -> Self::Inner
	{
		Self::each_map_inner([self, other], |[a, b]| f(a, b))
	}
	
	fn each_map_inner<F, const N: usize>(items: [Self; N], f: F) -> Self
	where
		F: Fn([Self::Inner; N]) -> Self::Inner;
	
	// !!! Can probably get rid of this method with some clever particular usage
	// of `Basis::map`. TBD
	fn inner_id(self) -> Self;
	
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
		fn each_map_inner<F, const N: usize>(items: [Self; N], f: F) -> Self
		where
			F: Fn([Self::Inner; N]) -> Self::Inner
		{
			f(items)
		}
		fn inner_id(self) -> Self {
			self
		}
	}
	
	impl<T, const N: usize> Basis for [T; N]
	where
		T: Basis,
	{
		type Inner = T::Inner;
		fn from_inner(inner: Self::Inner) -> Self {
			std::array::from_fn(|_| T::from_inner(inner.clone()))
		}
		fn each_map_inner<F, const M: usize>(items: [Self; M], f: F) -> Self
		where
			F: Fn([Self::Inner; M]) -> Self::Inner
		{
			let mut item_iters = items.map(IntoIterator::into_iter);
			std::array::from_fn(|_| T::each_map_inner(
				std::array::from_fn(|i| item_iters[i].next().unwrap()),
				&f
			))
		}
		fn inner_id(self) -> Self {
			self.map(T::inner_id)
		}
	}
	
	impl<A, B> Basis for Iso<A, B>
	where
		A: Basis,
		B: LinearIso<A>,
	{
		type Inner = A::Inner;
		fn from_inner(inner: Self::Inner) -> Self {
			let inner = A::from_inner(inner);
			Iso(Some(inner.clone()), LinearIso::<A>::from_linear(inner))
		}
		fn each_map_inner<F, const N: usize>(items: [Self; N], f: F) -> Self
		where
			F: Fn([Self::Inner; N]) -> Self::Inner
		{
			let inner = A::each_map_inner(items.map(Iso::into_inner), f);
			Iso(Some(inner.clone()), LinearIso::<A>::from_linear(inner))
		}
		fn inner_id(self) -> Self {
			Self::from_inner(B::linear_id(self.into_inner()))
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
pub trait LinearIso<T>: Sized + Clone {
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
	use super::{Linear, LinearIso, Simple};
	
	impl<T> LinearIso<T> for T
	where
		T: Linear + Simple
	{
		fn into_linear(value: Self) -> T {
			value
		}
		fn from_linear(value: T) -> Self {
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
	// `SumPoly<&T, D>: Poly`, etc.
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
	A: Vector<SIZE, Output: Basis>,
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
			impl Basis for $vec {
				type Inner = $value;
				
				fn from_inner(inner: Self::Inner) -> Self {
					<$vec>::splat(inner)
				}
				
				fn each_map_inner<F, const N: usize>(items: [Self; N], f: F) -> Self
				where
					F: Fn([Self::Inner; N]) -> Self::Inner
				{
					<$vec>::from_array(std::array::from_fn(|i| {
						let list = std::array::from_fn(|j| items[j][i]);
						Basis::each_map_inner(list, &f)
					}))
				}
				
				fn inner_id(self) -> Self {
					self
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
	use super::{Basis, Iso, Linear, LinearIso, Simple};
	
	impl<A, B> Iso<A, B>
	where
		A: Basis,
		B: LinearIso<A>,
	{
		pub fn into_inner(self) -> A {
			let Iso(inner, outer) = self;
			inner.unwrap_or_else(|| LinearIso::<A>::into_linear(outer))
		}
		
		pub fn from_inner(inner: A) -> Self {
			Iso(Some(inner.clone()), LinearIso::<A>::from_linear(inner))
		}
	}
	
	impl<A, B> Simple for Iso<A, B> {}
	
	impl<A, B> From<B> for Iso<A, B> {
		fn from(value: B) -> Self {
			Iso(None, value)
		}
	}
	
	impl<A, B> Iso<A, B> {
		pub fn new(outer: B) -> Self {
			Iso(None, outer)
		}
		
		pub fn map<T>(self, f: impl FnOnce(B) -> T) -> Iso<A, T> {
			Iso(self.0, f(self.1))
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
	
	impl<A, B> std::ops::Mul for Iso<A, B>
	where
		A: Basis,
		B: LinearIso<A>,
	{
		type Output = Self;
		fn mul(self, rhs: Self) -> Self::Output {
			self.zip_map_inner(rhs, Linear::mul)
		}
	}
	
	impl<A, B> crate::Flux for Iso<A, B>
	where
		A: Basis,
		B: LinearIso<A>,
	{
		type Basis = Self;
		type Change = crate::change::constant::Nil<Self>;
		fn basis(&self) -> Self::Basis {
			self.clone()
		}
		fn change(&self) -> Self::Change {
			self.accum().into_change()
		}
	}
	
	impl<A, B> crate::ToMoment for Iso<A, B>
	where
		A: Basis,
		B: LinearIso<A>,
	{
		type Moment<'a> = Self where Self: 'a;
		fn to_moment(&self, _time: f64) -> Self::Moment<'_> {
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