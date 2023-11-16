//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Mul, Sub};

use crate::linear::{Linear, Scalar};
use crate::{Time, Times, TimeRanges};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind:
	'static + Copy + Clone + Debug
	+ From<<Self as FluxKind>::Value>
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Accum<'a>: FluxAccum<'a, Self>;
	
	fn value(&self) -> Self::Value;
	
	// ??? Could instead make the polynomial evaluable at `-infinity`.
	fn initial_order(&self) -> Option<Ordering>
	where
		Self::Value: PartialOrd;
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

/// Combining [`FluxKind`] types.
/// 
/// Primarily this serves as a way to put two kinds of change-over-time into
/// the same space, for combination or comparison purposes.
pub mod ops {
	use std::ops;
	use super::FluxKind;
	use crate::linear::Scalar;
	
	/// Adding two kinds of change.
	pub trait Add<K: FluxKind = Self>: FluxKind<Value=K::Value> {
		type Output: FluxKind<Value=K::Value>;
		fn add(self, kind: K) -> <Self as Add<K>>::Output;
	}
	
	impl<A, B> Add<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value=A::Value>,
		<A as ops::Add<B>>::Output: FluxKind<Value=A::Value>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn add(self, kind: B) -> <A as ops::Add<B>>::Output {
			self + kind
		}
	}
	
	/// Differentiating two kinds of change.
	pub trait Sub<K: FluxKind = Self>: FluxKind<Value=K::Value> {
		type Output: FluxKind<Value=K::Value>;
		fn sub(self, kind: K) -> <Self as Sub<K>>::Output;
	}
	
	impl<A, B> Sub<B> for A
	where
		A: FluxKind + ops::Add<B>,
		B: FluxKind<Value=A::Value>,
		<A as ops::Add<B>>::Output: FluxKind<Value=A::Value>,
	{
		type Output = <A as ops::Add<B>>::Output;
		fn sub(self, kind: B) -> <A as ops::Add<B>>::Output {
			self + (kind * Scalar(-1.0))
		}
	}
	
	/// Squaring a kind of change.
	pub trait Sqr: FluxKind {
		type Output: FluxKind<Value=Self::Value>;
		fn sqr(self) -> <Self as Sqr>::Output;
	}
	
	impl<K: FluxKind> Sqr for K
	where
		K: ops::Mul,
		<K as ops::Mul>::Output: FluxKind<Value=K::Value>,
	{
		type Output = <K as ops::Mul>::Output;
		fn sqr(self) -> <Self as Sqr>::Output {
			self * self
		}
	}
}

/// Change accumulator.
/// 
/// Converts a discrete pattern of change into a desired form.
pub trait FluxAccum<'a, K: FluxKind> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self;
}

impl<K: FluxKind> FluxAccum<'_, K> for () {
	fn from_kind(_kind: FluxAccumKind<'_, K>) -> Self {}
}

/// General accumulator arguments.
#[non_exhaustive]
pub enum FluxAccumKind<'a, K: FluxKind> {
	Value {
		value: &'a mut K::Value,
		depth: usize,
		time: Time,
		base_time: Time,
	},
	Poly {
		poly: &'a mut Poly<K>,
		depth: usize,
		time: Time,
		base_time: Time,
	},
}

/// A polynomial wrapper for [`FluxKind`].
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<K>(K);

impl<K: FluxKind> Poly<K> {
	pub fn with_value(value: K::Value) -> Self {
		Poly(K::from(value))
	}
	
	pub fn sqr(self) -> Poly<<K as ops::Sqr>::Output>
	where
		K: ops::Sqr
	{
		(*self).sqr().into()
	}
	
	pub fn real_roots(self) -> Result<RootList, RootList>
	where
		K: Roots
	{
		//! Returns all real-valued roots of this polynomial in ascending order.
		//! If not all roots are known, `Err` is returned.
		
		let cleanup = |roots: RootList| {
			let mut roots = roots.into_vec();
			roots.retain(|r| !r.is_nan());
			roots.sort_unstable_by(|a,b| a.total_cmp(b));
			roots.into_boxed_slice()
		};
		
		self.roots()
			.map(cleanup)
			.map_err(cleanup)
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	pub(crate) fn when_sign(&self, order: Ordering, offset: Time) -> TimeRanges
	where
		K: Roots + PartialOrd,
		K::Value: PartialOrd,
	{
		let initial_order = self.initial_order();
		
		 // Convert Roots to Ranges:
		let range_list = match self.real_roots() {
			Ok(roots) => {
				let mut list = Vec::with_capacity(
					if order == Ordering::Equal {
						roots.len()
					} else {
						1 + (roots.len() / 2)
					}
				);
				let mut prev_point = if initial_order == Some(order) {
					Some(f64::NEG_INFINITY)
				} else {
					None
				};
				for &point in roots.iter() {
					if let Some(prev) = prev_point {
						if point != prev {
							list.push((prev, point));
						}
						prev_point = None;
					} else if order == Ordering::Equal {
						list.push((point, point));
					} else {
						prev_point = Some(point);
					}
				}
				if let Some(point) = prev_point {
					list.push((point, f64::INFINITY));
				}
				list.into_iter()
			},
			Err(_) => Default::default(),
		};
		
		 // Convert Ranges to Times:
		let time = offset.as_secs_f64();
		let max_time = Time::MAX.as_secs_f64();
		range_list
			.filter_map(|(a, b)| {
				debug_assert!(a <= b);
				let a = a + time;
				let b = b + time;
				if b < 0.0 || a >= max_time {
					None
				} else {
					Some((
						Time::try_from_secs_f64(a).unwrap_or(Time::ZERO),
						Time::try_from_secs_f64(b).unwrap_or(Time::MAX)
					))
				}
			})
			.collect()
	}
	
	/// Times when the value is equal to zero.
	pub(crate) fn when_zero(&self, offset: Time) -> Times
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		let real_roots = self.real_roots()
			.unwrap_or_default()
			.into_vec();
		
		 // Constant Equality:
		if real_roots.is_empty() && self.is_zero() {
			return [Time::ZERO]
				.into_iter().collect();
		}
		
		 // Convert Roots to Times:
		let time = offset.as_secs_f64();
		let max_time = Time::MAX.as_secs_f64();
		real_roots.into_iter()
			.filter_map(|t| {
				let t = t + time;
				if t < 0.0 || t >= max_time {
					None
				} else {
					Some(Time::from_secs_f64(t))
				}
			})
			.collect()
	}
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Self(K::zero())
	}
}

impl<K: FluxKind> From<K> for Poly<K> {
	fn from(value: K) -> Self {
		Poly(value)
	}
}

impl<K> Deref for Poly<K> {
	type Target = K;
	fn deref(&self) -> &Self::Target {
		&self.0
	}
}

impl<K> DerefMut for Poly<K> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.0
	}
}

impl<A: FluxKind, B: FluxKind> Add<Poly<B>> for Poly<A>
where
	A: ops::Add<B>
{
	type Output = Poly<<A as ops::Add<B>>::Output>;
	fn add(self, rhs: Poly<B>) -> Self::Output {
		Poly((*self).add(*rhs))
	}
}

impl<A: FluxKind, B: FluxKind> Sub<Poly<B>> for Poly<A>
where
	A: ops::Sub<B>
{
	type Output = Poly<<A as ops::Sub<B>>::Output>;
	fn sub(self, rhs: Poly<B>) -> Self::Output {
		Poly((*self).sub(*rhs))
	}
}

impl<K: FluxKind> Mul<Scalar> for Poly<K> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		*self = *self * rhs;
		self
	}
}

// !!! This should be an iterator at some point. Need to be able to produce a
// generating function for roots.
pub type RootList = Box<[f64]>;

/// Roots of a [`Poly`]nomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: FluxKind {
	fn roots(self) -> Result<RootList, RootList>;
}