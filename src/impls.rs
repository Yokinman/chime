//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use crate::{Constant, Flux, ToMoment, ToMomentMut};
use crate::linear::{BasisArray, Scalar};

macro_rules! impl_ {
    ($trait_:ident for $type_:ty) => {
		impl $trait_ for $type_ {
			impl_!(@impl $trait_);
		}
    };
	
	 // Split Type Path Tree:
    ($trait_:tt for $($p:ident::)* {$($t:ident),* $(,)?}) => {
		impl_!(@type_split $trait_ ($($p::)*): {$($t),*});
    };
	(@type_split $trait_:tt $p:tt: {$($t:ident),*}) => {
		$(impl_!(@type_join $trait_ $p: $t);)*
	};
	(@type_join $trait_:tt ($($p:ident::)*): $t:ident) => {
		impl_!($trait_ for $($p::)* $t);
	};
	
	 // Split Trait Path Tree:
    ({$($trait_:ident),* $(,)?} for $type_:ty) => {
	    $(impl_!($trait_ for $type_);)*
    };
	
	 // Impl Cases:
    (@impl Flux) => {
		type Basis = Self;
		type Kind = Constant<Self>;
		fn basis(&self) -> Self::Basis {
			*self
		}
		fn change(&self, basis: Constant<Self::Basis>) -> Self::Kind {
			basis
		}
    };
    (@impl ToMoment) => {
		type Moment<'a> = Self;
		fn to_moment(&self, _time: Scalar) -> Self::Moment<'_> {
			*self
		}
    };
    (@impl ToMomentMut) => {
		type MomentMut<'a> = &'a mut Self;
		fn to_moment_mut(&mut self, _time: Scalar) -> Self::MomentMut<'_> {
			self
		}
    };
}

impl_!({ToMoment, ToMomentMut} for bool);
impl_!({ToMoment, ToMomentMut} for char);
impl_!({Flux, ToMoment, ToMomentMut} for {f32, f64});
impl_!({ToMoment, ToMomentMut} for {i8, i16, i32, i64, i128, isize});
impl_!({ToMoment, ToMomentMut} for {u8, u16, u32, u64, u128, usize});

impl_!({ToMoment, ToMomentMut} for std::time::Duration);

#[cfg(feature = "glam")]
impl_!({Flux, ToMoment, ToMomentMut} for glam::{
	Vec2, Vec3, Vec4,
	DVec2, DVec3, DVec4,
});

#[cfg(feature = "glam")]
impl_!({ToMoment, ToMomentMut} for glam::{
	IVec2, IVec3, IVec4,
	UVec2, UVec3, UVec4,
});

mod _reference_impls {
	use crate::Constant;
	use super::*;
	
	impl<'t, T> Flux for &'t T
	where
		T: Flux,
	{
		type Basis = T::Basis;
		type Kind = T::Kind;
		fn basis(&self) -> Self::Basis {
			T::basis(self)
		}
		fn change(&self, basis: Constant<Self::Basis>) -> Self::Kind {
			T::change(self, basis)
		}
	}
	
	impl<'t, T> ToMoment for &'t T
	where
		T: ToMoment,
	{
		type Moment<'a> = T::Moment<'a> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			T::to_moment(self, time)
		}
	}
	
	impl<'t, T> Flux for &'t mut T
	where
		T: Flux,
	{
		type Basis = T::Basis;
		type Kind = T::Kind;
		fn basis(&self) -> Self::Basis {
			T::basis(self)
		}
		fn change(&self, basis: Constant<Self::Basis>) -> Self::Kind {
			T::change(self, basis)
		}
	}
	
	impl<'t, T> ToMoment for &'t mut T
	where
		T: ToMoment,
	{
		type Moment<'a> = T::Moment<'a> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			T::to_moment(self, time)
		}
	}
	
	impl<'t, T> ToMomentMut for &'t mut T
	where
		T: ToMomentMut,
	{
		type MomentMut<'a> = T::MomentMut<'a> where Self: 'a;
		fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
			T::to_moment_mut(self, time)
		}
	}
}

impl<T> ToMoment for Option<T>
where
	T: ToMoment,
{
	type Moment<'a> = Option<T::Moment<'a>> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		self.as_ref().map(|x| x.to_moment(time))
	}
}

impl<T, E> ToMoment for Result<T, E>
where
	T: ToMoment,
	E: ToMoment,
{
	type Moment<'a> = Result<T::Moment<'a>, E::Moment<'a>> where Self: 'a;
	fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
		self.as_ref()
			.map(|x| x.to_moment(time))
			.map_err(|x| x.to_moment(time))
	}
}

mod _range_impls {
	use super::*;
	
	impl_!(ToMoment for std::ops::RangeFull);
	
	impl<T> ToMoment for std::ops::Range<T>
	where
		T: ToMoment,
	{
		type Moment<'a> = std::ops::Range<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			std::ops::Range {
				start: self.start.to_moment(time),
				end: self.end.to_moment(time),
			}
		}
	}
	
	impl<T> ToMoment for std::ops::RangeFrom<T>
	where
		T: ToMoment,
	{
		type Moment<'a> = std::ops::RangeFrom<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			std::ops::RangeFrom {
				start: self.start.to_moment(time),
			}
		}
	}
	
	impl<T> ToMoment for std::ops::RangeInclusive<T>
	where
		T: ToMoment,
	{
		type Moment<'a> = std::ops::RangeInclusive<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			std::ops::RangeInclusive::new(
				self.start().to_moment(time),
				self.end().to_moment(time),
			)
		}
	}
	
	impl<T> ToMoment for std::ops::RangeTo<T>
	where
		T: ToMoment,
	{
		type Moment<'a> = std::ops::RangeTo<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			std::ops::RangeTo {
				end: self.end.to_moment(time),
			}
		}
	}
	
	impl<T> ToMoment for std::ops::RangeToInclusive<T>
	where
		T: ToMoment,
	{
		type Moment<'a> = std::ops::RangeToInclusive<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			std::ops::RangeToInclusive {
				end: self.end.to_moment(time),
			}
		}
	}
}

mod _array_impls {
	use crate::Constant;
	use super::*;
	
	impl<T: Flux, const N: usize> Flux for [T; N] {
		type Basis = BasisArray<T::Basis, N>;
		type Kind = [T::Kind; N];
		fn basis(&self) -> Self::Basis {
			BasisArray(self.each_ref().map(T::basis))
		}
		fn change(&self, basis: Constant<Self::Basis>) -> Self::Kind {
			let mut basis_iter = basis.into_iter();
			self.each_ref()
				.map(|x| x.change(basis_iter.next().unwrap()))
		}
	}
	
	impl<T: ToMoment, const N: usize> ToMoment for [T; N] {
		type Moment<'a> = [T::Moment<'a>; N] where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			self.each_ref().map(|x| x.to_moment(time))
		}
	}
	
	impl<T: ToMomentMut, const N: usize> ToMomentMut for [T; N] {
		type MomentMut<'a> = [T::MomentMut<'a>; N] where Self: 'a;
		fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
			self.each_mut().map(|x| x.to_moment_mut(time))
		}
	}
}

mod _vec_impls {
	use std::vec::Vec;
	use super::*;
	
	impl<T: ToMoment> ToMoment for Vec<T> {
		type Moment<'a> = Vec<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			self.iter()
				.map(|x| x.to_moment(time))
				.collect()
		}
	}
	
	impl<T: ToMomentMut> ToMomentMut for Vec<T> {
		type MomentMut<'a> = Vec<T::MomentMut<'a>> where Self: 'a;
		fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
			self.iter_mut()
				.map(|x| x.to_moment_mut(time))
				.collect()
		}
	}
}

mod _hashmap_impls {
	use std::collections::HashMap;
	use super::*;
	
	impl<K, V> ToMoment for HashMap<K, V>
	where
		K: std::hash::Hash + Eq,
		V: ToMoment,
	{
		type Moment<'a> = HashMap<&'a K, V::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			self.iter()
				.map(|(k, v)| (k, v.to_moment(time)))
				.collect()
		}
	}
	
	impl<K, V> ToMomentMut for HashMap<K, V>
	where
		K: std::hash::Hash + Eq,
		V: ToMomentMut,
	{
		type MomentMut<'a> = HashMap<&'a K, V::MomentMut<'a>> where Self: 'a;
		fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
			self.iter_mut()
				.map(|(k, v)| (k, v.to_moment_mut(time)))
				.collect()
		}
	}
}

mod _tuple_impls {
	use super::*;
	
	macro_rules! impl_to_moment_for_tuples {
	    (@impl $($t:ident),*) => {
			impl<$($t,)*> ToMoment for ($($t,)*)
			where
				$($t: ToMoment,)*
			{
				type Moment<'a> = ($($t::Moment<'a>,)*) where Self: 'a;
				#[allow(unused_variables)]
				fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
					let ($($t,)*) = self;
					($($t.to_moment(time),)*)
				}
			}
	    };
		($($a:ident $(, $b:ident)*)?) => {
			$(impl_to_moment_for_tuples!($($b),*);)?
			impl_to_moment_for_tuples!(@impl $($a $(, $b)*)?);
		};
	}
	
	impl_to_moment_for_tuples!(_0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);
}