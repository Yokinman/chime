//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::Mul;

use crate::linear::{Linear, Basis, Scalar};
use crate::time::Time;
use crate::{Flux, ToMomentMut};

pub mod constant;
pub mod sum;

/// An abstract description of change over time.
/// 
/// Used to define the standard representations of [`Flux`] types. In
/// other words, the layout of a polynomial.
pub trait FluxKind: Flux<Kind=Self> + FromChange<Self::Change> + ToMomentMut + Clone + Debug + 'static {
	const DEGREE: usize;
	
	fn with_basis(value: Self::Basis) -> Self;
	
	fn add_basis(self, value: Self::Basis) -> Self;
	
	fn sub_basis(self, value: Self::Basis) -> Self {
		self.add_basis(Basis::from_inner(value.into_inner()
			.mul_scalar(Scalar::from(-1.))))
	}
	
	fn eval(&self, time: Scalar) -> Self::Basis;
	
	/// The order at or immediately preceding the value at a time.
	/// 
	/// This should be the first non-zero [`FluxKind::eval`] value of this kind
	/// or its derivatives; reversed for odd derivatives.
	fn initial_order(&self, time: Scalar) -> Option<Ordering>
	where
		<Self::Basis as Basis>::Inner: PartialOrd
	{
		// !!! Alternative: Translate polynomial using `to_moment_mut` and then
		// check leading terms in order. Unknown which is more precise/faster.
		
		let order = self.eval(time)
			.into_inner()
			.partial_cmp(&Linear::zero());
			
		if order != Some(Ordering::Equal) || Self::DEGREE == 0 {
			return order
		}
		
		self.change()
			.initial_order(time)
			.map(Ordering::reverse)
	}
	
	fn zero() -> Self {
		Self::with_basis(Basis::zero())
	}
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

/// A [`FluxKind`] that can be integrated into a higher degree of change.
/// 
/// e.g. `Sum<T, 1>::Integ == Sum<T, 2>`, `Cos<T>::Integ == Sin<T>`.
/// 
/// Used for the `std::ops::{Add, Sub}` impls of [`FluxAccum`].
pub trait FluxIntegral: FluxKind + Mul<Scalar, Output=Self> {
	type Integ: FluxKind<Basis = Self::Basis>;
	fn integ(self, basis: Self::Basis) -> Self::Integ;
}

impl<T: FluxKind, const SIZE: usize> FluxKind for [T; SIZE] {
	const DEGREE: usize = T::DEGREE;
	fn with_basis(value: Self::Basis) -> Self {
		value.map(T::with_basis)
	}
	fn add_basis(self, value: Self::Basis) -> Self {
		let mut values = value.into_iter();
		self.map(|x| x.add_basis(values.next().unwrap()))
	}
	fn eval(&self, time: Scalar) -> Self::Basis {
		self.each_ref().map(|x| T::eval(x, time))
	}
}

impl<F, T, const N: usize> FromChange<[T; N]> for [F; N]
where
	F: FromChange<T>,
	T: FluxKind,
{
	fn from_change(basis: [T::Basis; N], change: [T; N]) -> Self {
		let mut change_iter = change.into_iter();
		basis.map(|x| F::from_change(x, change_iter.next().unwrap()))
	}
}

/// Shortcut for the inner [`Linear`] type of a [`FluxKind`].
pub(crate) type KindLinear<T> = <<T as Flux>::Basis as Basis>::Inner;

/// ...
pub trait FromChange<T: FluxKind> {
	fn from_change(basis: <T as Flux>::Basis, change: T) -> Self;
}

/// Combining [`FluxKind`] types.
/// 
/// Primarily this serves as a way to put two kinds of change-over-time into
/// the same space, for combination or comparison purposes.
pub mod ops {
	use std::ops;
	use super::{FluxKind, KindLinear};
	use crate::linear::Basis;
	
	/// Squaring a kind of change.
	pub trait Sqr: FluxKind {
		type Output: FluxKind<Basis: Basis<Inner = KindLinear<Self>>>;
		fn sqr(self) -> <Self as Sqr>::Output;
	}
	
	impl<K: FluxKind> Sqr for K
	where
		K: ops::Mul,
		<K as ops::Mul>::Output: FluxKind<Basis: Basis<Inner = KindLinear<K>>>,
	{
		type Output = <K as ops::Mul>::Output;
		fn sqr(self) -> <Self as Sqr>::Output {
			self.clone() * self
		}
	}
}

/// Roots of a [`FluxKind`] polynomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: FluxKind {
	type Output: IntoTimes;
	fn roots(self) -> <Self as Roots>::Output;
}

/// Conversion from some [`Roots::Output`] into an iterator of time.
pub trait IntoTimes {
	type TimeIter: Iterator<Item=LinearTime>;
	fn into_times(self) -> Self::TimeIter;
}

impl<const N: usize> IntoTimes for [f64; N] {
	type TimeIter = std::array::IntoIter<LinearTime, N>;
	fn into_times(self) -> Self::TimeIter {
		let mut times = self.map(LinearTime::from_secs_f64);
		times.sort_unstable();
		times.into_iter()
	}
}

/// Interface for creating [`Time`] values to override conversion.
#[derive(Copy, Clone, Debug, Default)]
pub struct LinearTime(f64);

impl LinearTime {
	pub fn from_secs_f64(secs: f64) -> Self {
		Self(secs)
	}
	
	/// Conversion into [`Time`], but always rounds down.
	/// 
	/// Consistently rounding down is important so that tiny ranges of time
	/// aren't ignored due to rounding. It's also better if predictions catch
	/// events before they happen rather than after.
	pub(crate) fn try_into_time(self, basis: Time) -> Result<Time, Time> {
		let LinearTime(mut t) = self;
		let sign = t.signum();
		t *= sign;
		
		const MANT_MASK: u64 = (1 << 52) - 1;
		const EXP_MASK: u64 = (1 << 11) - 1;
		
		let bits = t.to_bits();
		let mant = (bits & MANT_MASK) | (MANT_MASK + 1);
		let exp = (((bits >> 52) & EXP_MASK) as i16) - 1023;
		
		let time = if exp < -30 {
			// Too small; `1ns < 2s^-30`.
			if sign == -1. && t != 0. {
				Time::new(0, 1)
			} else {
				Time::ZERO
			}
		} else if exp < 0 {
			// No integer part.
			let nanos_tmp = ((mant as u128) << (44 + exp)) * 1_000_000_000;
			let mut nanos = (nanos_tmp >> (44 + 52)) as u32;
			if sign == -1. && (t * 1e9).fract() != 0. {
				nanos += 1;
			}
			Time::new(0, nanos)
		} else if exp < 52 {
			let secs = mant >> (52 - exp);
			let nanos_tmp = (((mant << exp) & MANT_MASK) as u128) * 1_000_000_000;
			let mut nanos = (nanos_tmp >> 52) as u32;
			if sign == -1. && (t * 1e9).fract() != 0. {
				nanos += 1;
			}
			Time::new(secs, nanos)
		} else if exp < 64 {
			// No fractional part.
			Time::from_secs(mant << (exp - 52))
		} else {
			// Too big.
			return if sign == -1. {
				Err(Time::ZERO)
			} else {
				Err(Time::MAX)
			}
		};
		
		if sign == -1. {
			basis.checked_sub(time).ok_or(Time::ZERO)
		} else {
			basis.checked_add(time).ok_or(Time::MAX)
		}
	}
}

impl Ord for LinearTime {
	fn cmp(&self, other: &Self) -> Ordering {
		self.0.total_cmp(&other.0)
	}
}

impl PartialOrd for LinearTime {
	fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		Some(self.cmp(other))
	}
}

impl Eq for LinearTime {}

impl PartialEq for LinearTime {
	fn eq(&self, other: &Self) -> bool {
		self.cmp(other) == Ordering::Equal
	}
}

/// ...
#[derive(Clone)]
pub struct RootFilterMap<T> {
	pub(crate) times: T,
	pub(crate) basis: Time,
	pub(crate) prev_time: Time,
}

impl<T> Iterator for RootFilterMap<T>
where
	T: Iterator<Item = LinearTime>,
{
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		for root in self.times.by_ref() {
			if let Ok(time) = root.try_into_time(self.basis) {
				let prev_time = self.prev_time;
				self.prev_time = time;
				return Some(time - prev_time)
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.times.size_hint().1)
	}
}