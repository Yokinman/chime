//! Defining a *kind* of change over time.

use super::*;

use std::fmt::Debug;
use std::ops::{Add, Mul, Sub};

/// Defines a kind of change as the structure of a polynomial.
/// 
/// ??? Should `Poly` just be absorbed into this trait? It feels like it would
/// be convenient if the constant were managed by these types. 
pub trait FluxKind:
	'static + Copy + Clone + Debug
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Accum<'a>: FluxAccum<'a, Self>;
	
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

/// Change accumulator.
/// 
/// Converts a discrete pattern of change into a desired form.
pub trait FluxAccum<'a, K: FluxKind> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self;
}

/// General accumulator arguments.
#[non_exhaustive]
pub enum FluxAccumKind<'a, K: FluxKind> {
	Value {
		value: &'a mut K::Value,
		depth: usize,
		time: Time,
		offset: Time,
	},
	Poly {
		poly: &'a mut Poly<K>,
		depth: usize,
	},
}

/// A polynomial in standard form; e.g. `a + b x + c x^2 + d x^3`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<K: FluxKind>(pub K::Value, pub K); // !!! Use named fields.

impl<K: FluxKind> Poly<K> {
	pub fn constant(&self) -> K::Value {
		self.0
	}
	
	pub fn is_zero(&self) -> bool
	where
		K: PartialEq,
		K::Value: PartialEq,
	{
		self.0.is_zero() && self.1.is_zero()
	}
	
	// ??? Could instead make the polynomial evaluable at `-infinity`.
	pub fn initial_order(&self) -> Option<Ordering>
	where
		K::Value: PartialOrd
	{
		match self.1.initial_order() {
			Some(Ordering::Equal) => self.0.partial_cmp(&K::Value::zero()),
			order => order
		}
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
		
		K::roots(self)
			.map(cleanup)
			.map_err(cleanup)
	}
	
	pub fn when(&self, order: Ordering, offset: Time) -> TimeRanges
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
	
	pub fn when_eq(&self, offset: Time) -> Times
	where
		K: Roots + PartialEq,
		K::Value: PartialEq,
	{
		let mut real_roots = self.real_roots()
			.unwrap_or_default()
			.into_vec();
		
		 // Constant Equality:
		if real_roots.is_empty() && self.is_zero() {
			real_roots.push(0.0);
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
		Self(Linear::zero(), K::zero())
	}
}

impl<A, B> Add<Poly<B>> for Poly<A>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Value=A::Value>,
	<A as Add<B>>::Output: FluxKind<Value=A::Value>,
	A::Value: Add<B::Value, Output=A::Value>,
{
	type Output = Poly<<A as Add<B>>::Output>;
	fn add(self, rhs: Poly<B>) -> Self::Output {
		Poly(
			self.0 + rhs.0,
			self.1 + rhs.1,
		)
	}
}

impl<A, B> Sub<Poly<B>> for Poly<A>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Value=A::Value>,
	<A as Add<B>>::Output: FluxKind<Value=A::Value>,
	A::Value: Add<B::Value, Output=A::Value>,
{
	type Output = Poly<<A as Add<B>>::Output>;
	fn sub(self, rhs: Poly<B>) -> Self::Output {
		self + (rhs * Scalar(-1.0))
	}
}

impl<K: FluxKind> Mul<Scalar> for Poly<K> {
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.0 = self.0 * rhs;
		self.1 = self.1 * rhs;
		self
	}
}

/// Squaring [`Poly`]nomials.
pub trait PolySqr: private::Sealed {
	type Output: FluxKind;
	fn sqr(self) -> Poly<Self::Output>;
}

impl<K: FluxKind> private::Sealed for Poly<K> {}

impl<K: FluxKind> PolySqr for Poly<K>
where
	K: Mul + Mul<K::Value, Output=K> + Add<<K as Mul>::Output>,
	K::Value: Mul<Output=K::Value>,
	<K as Add<<K as Mul>::Output>>::Output: FluxKind<Value=K::Value>
{
	type Output = <K as Add<<K as Mul>::Output>>::Output;
	fn sqr(self) -> Poly<Self::Output> {
		Poly(
			self.0*self.0,
			self.1*self.0*Scalar(2.0) + self.1*self.1,
		)
	}
}

mod private {
	pub trait Sealed {}
}

// !!! This should be an iterator at some point. Need to be able to produce a
// generating function for roots.
pub type RootList = Box<[f64]>;

/// Roots of a [`Poly`]nomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: FluxKind {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList>;
}