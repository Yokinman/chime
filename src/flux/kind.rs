//! Defining a *kind* of change over time.

use super::*;

use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::ops::{Add, Index, IndexMut, Mul, Shr, Sub};
use std::slice::{Iter, IterMut};

/// Defines a kind of change as the structure of a polynomial.
pub trait FluxKind:
	'static + Copy + Clone + Debug
	+ Mul<Scalar, Output=Self>
{
	type Value: Linear;
	
	type Coeffs:
		Copy + Clone + Debug
		+ IntoIterator<Item=Self>
		+ AsRef<[Self]> + AsMut<[Self]>
		+ Borrow<[Self]> + BorrowMut<[Self]>
		+ Index<usize, Output=Self> + IndexMut<usize>;
	
	type Accum<'a>: FluxAccum<'a, Self>;
	
	const DEGREE: usize;
	
	fn zero() -> Self;
	
	fn zero_coeffs() -> Self::Coeffs;
	
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

/// Used with `Shr` and `Shl` for upgrading/downgrading a [`FluxKind`].
#[derive(Copy, Clone, Debug)]
pub struct DegShift;
// !!! assert { K::DEGREE + 1 } == <K as Shr<DegShift>>::Output::DEGREE ?

/// A polynomial in standard form; e.g. `a + b x + c x^2 + d x^3`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Poly<K: FluxKind>(pub K::Value, pub K::Coeffs);

impl<K: FluxKind> Poly<K> {
	pub fn constant(&self) -> K::Value {
		self.0
	}
	
	pub fn coeff(&self, coeff_index: usize) -> Option<K> {
		self.1.as_ref().get(coeff_index).copied()
	}
	
	pub fn coeff_iter(&self) -> Iter<K> {
		self.1.as_ref().iter()
	}
	
	pub fn coeff_iter_mut(&mut self) -> IterMut<K> {
		self.1.as_mut().iter_mut()
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
}

impl<K: FluxKind> Default for Poly<K> {
	fn default() -> Self {
		Self(Linear::zero(), K::zero_coeffs())
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
		let mut poly = Poly {
			0: self.constant() + rhs.constant(),
			..Default::default()
		};
		let mut a = self.coeff_iter();
		let mut b = rhs.coeff_iter();
		for i in 0..<<A as Add<B>>::Output as FluxKind>::DEGREE {
			let coeff = IndexMut::index_mut(&mut poly.1, i);
			*coeff = *a.next().unwrap_or(&A::zero())
				+ *b.next().unwrap_or(&B::zero())
		}
		poly
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
		self.0 = self.constant() * rhs;
		for term in self.coeff_iter_mut() {
			*term = *term * rhs
		}
		self
	}
}

impl<K: FluxKind> Shr<DegShift> for Poly<K>
where
	K: Shr<DegShift>,
	<K as Shr<DegShift>>::Output: FluxKind<Value=K::Value>,
	K: From<K::Value>,
{
	type Output = Poly<<K as Shr<DegShift>>::Output>;
	fn shr(self, rhs: DegShift) -> Self::Output {
		let mut poly = Poly::default();
		let coeff = IndexMut::index_mut(&mut poly.1, 0);
		*coeff = K::from(self.0) >> rhs;
		for d in 0..K::DEGREE {
			let coeff = IndexMut::index_mut(&mut poly.1, d + 1);
			*coeff = self.1[d] >> rhs;
		}
		poly
	}
}

pub type RootList = Box<[f64]>;

/// Roots of a [`Poly`]nomial.
pub trait Roots: FluxKind {
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList>;
}