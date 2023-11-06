//! Defining a *kind* of change over time.

use super::*;

use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::ops::{Index, IndexMut, Mul};

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

/// Used with `Shr` and `Shl` for upgrading/downgrading a [`FluxKind`].
#[derive(Copy, Clone, Debug)]
pub struct DegShift;
// !!! assert { K::DEGREE + 1 } == <K as Shr<DegShift>>::Output::DEGREE ?