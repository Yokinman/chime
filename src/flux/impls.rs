//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::{Flux, ToMoment, ToMomentMut};
use crate::linear::{BasisArray, Scalar};

macro_rules! impl_to_moment {
	($($p:ident::)* {$($t:ident),* $(,)?}) => {
		impl_to_moment!(@pack_path ($($p::)*): {$($t),*});
	};
	(@pack_path $p:tt: {$($t:ident),*}) => {
		$(impl_to_moment!(@unpack_path $p: $t);)*
	};
	(@unpack_path ($($p:ident::)*): $t:ident) => {
		impl_to_moment!($($p::)* $t);
	};
    ($t:ty) => {
		impl ToMoment for $t {
			type Moment<'a> = Self;
			fn to_moment(&self, _time: Scalar) -> Self::Moment<'_> {
				*self
			}
		}
    };
}

impl_to_moment!(bool);
impl_to_moment!(char);
impl_to_moment!({f32, f64});
impl_to_moment!({i8, i16, i32, i64, i128, isize});
impl_to_moment!({u8, u16, u32, u64, u128, usize});

impl_to_moment!(std::time::Duration);

#[cfg(feature = "glam")]
impl_to_moment!(glam::{
	 Vec2,  Vec3,  Vec4,
	DVec2, DVec3, DVec4,
	IVec2, IVec3, IVec4,
	UVec2, UVec3, UVec4,
});

mod _reference_impls {
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
		fn change(&self, kind: Self::Kind) -> Self::Kind {
			T::change(self, kind)
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
		fn change(&self, kind: Self::Kind) -> Self::Kind {
			T::change(self, kind)
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
	
	impl_to_moment!(std::ops::RangeFull);
	
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
	use super::*;
	
	impl<T: Flux, const SIZE: usize> Flux for [T; SIZE] {
		type Basis = BasisArray<T::Basis, SIZE>;
		type Kind = [T::Kind; SIZE];
		fn basis(&self) -> Self::Basis {
			BasisArray(self.each_ref().map(T::basis))
		}
		fn change(&self, kind: Self::Kind) -> Self::Kind {
			let mut kind_iter = kind.into_iter();
			self.each_ref()
				.map(|x| x.change(kind_iter.next().unwrap()))
		}
	}
	
	impl<T: ToMoment, const SIZE: usize> ToMoment for [T; SIZE] {
		type Moment<'a> = [T::Moment<'a>; SIZE] where Self: 'a;
		fn to_moment(&self, time: Scalar) -> Self::Moment<'_> {
			self.each_ref().map(|x| x.to_moment(time))
		}
	}
	
	impl<T: ToMomentMut, const SIZE: usize> ToMomentMut for [T; SIZE] {
		type MomentMut<'a> = [T::MomentMut<'a>; SIZE] where Self: 'a;
		fn to_moment_mut(&mut self, time: Scalar) -> Self::MomentMut<'_> {
			self.each_mut().map(|x| x.to_moment_mut(time))
		}
	}
}

mod _vec_impls {
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