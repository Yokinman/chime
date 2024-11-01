//! Convenient [`Flux`] implementations (`Vec<T>`, `[T; S]`, ??? tuples, etc.).

use std::collections::HashMap;
use std::vec::Vec;

use crate::{Flux, ToMoment, ToMomentMut};
use crate::linear::{BasisArray, Scalar};

macro_rules! impl_to_moment {
	({$($t:ty),*}) => {
		$(impl_to_moment!($t);)*
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

impl_to_moment!({u8, u16, u32, u64, u128, usize});
impl_to_moment!({i8, i16, i32, i64, i128, isize});
impl_to_moment!({f32, f64});

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