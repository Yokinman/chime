//! ...

use std::fmt::Debug;
use std::ops::{Add, Mul};

use impl_op::impl_op;

use crate::kind::{Deg, FluxKind};
use crate::polynomial::{Poly, RootList, Roots};

/// A scalar value, used for multiplication with any [`LinearValue`].
#[derive(Copy, Clone, Debug)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) => Scalar(a.0 * b.0),
	(f64, Scalar)['commut] => a * b.0,
	(f32, Scalar)['commut] => a * (b.0 as f32),
	// (u64, Scalar)['commut] |
	// (u32, Scalar)['commut] |
	// (u16, Scalar)['commut] |
	// (u8,  Scalar)['commut] |
	// (i64, Scalar)['commut] |
	// (i32, Scalar)['commut] |
	// (i16, Scalar)['commut] |
	// (i8,  Scalar)['commut] => ((a as f64) * b.0).round() as Self,
}}

/// Any vector type that has addition and [`Scalar`] multiplication.
pub trait LinearValue:
	Sized + Copy + Clone + Debug
	+ Add<Output=Self>
	+ Mul<Scalar, Output=Self>
{
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

impl<T> LinearValue for T
where T:
	Sized + Copy + Clone + Debug
	+ Add<Output=T>
	+ Mul<Scalar, Output=T>
	+ Default
{
	fn zero() -> Self {
		Self::default()
	}
}

/// A mapping of a vector space that preserves addition & multiplication.
/// ??? Maybe rename to LinearMap, since it's not *really* isomorphic.
/// 
/// # Properties
/// 
/// - Generally isomorphic       - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Maps vector addition       - `map(A + B) = map(A) • map(B)`
/// - Maps scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso<Mapped>: LinearValue {
	fn map(self) -> Mapped;
	fn inv_map(value: Mapped) -> Self;
	fn set(&mut self, value: Mapped) {
		*self = Self::inv_map(value);
	}
}

impl<T: LinearValue> LinearIso<T> for T {
	fn map(self) -> T {
		self
	}
	fn inv_map(value: T) -> Self {
		value
	}
}

macro_rules! impl_iso_for_int {
	($a:ty: $($b:ty),+) => {$(
		impl LinearIso<$b> for $a {
			fn map(self) -> $b {
				self.round() as $b
			}
			fn inv_map(value: $b) -> Self {
				value as $a
			}
		}
	)+}
}
impl_iso_for_int!(f32: u8, u16, u32, u64, i8, i16, i32, i64);
impl_iso_for_int!(f64: u8, u16, u32, u64, i8, i16, i32, i64);

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exp(A * B)`
#[derive(Copy, Clone, Debug, PartialOrd, PartialEq)]
pub struct Exp<T: LinearValue>(T);

impl<T: LinearValue> Add for Exp<T> {
	type Output = Self;
	fn add(self, rhs: Self) -> Self {
		Self(self.0 + rhs.0)
	}
}

impl<T: LinearValue> Mul<Scalar> for Exp<T> {
	type Output = Self;
	fn mul(self, rhs: Scalar) -> Self {
		Self(self.0 * rhs)
	}
}

impl<T: LinearValue> Default for Exp<T> {
	fn default() -> Self {
		Self(T::zero())
	}
}

impl<T: LinearValue> From<T> for Exp<T> {
	fn from(value: T) -> Exp<T> {
		Exp(value)
	}
}

impl<T: LinearValue, const DEG: usize> Roots for Deg<Exp<T>, DEG>
where
	Deg<T, DEG>: FluxKind<Linear=T> + Roots
{
	fn roots(poly: Poly<Self>) -> Result<RootList, RootList> {
		let mut b_poly = Poly::<Deg<T, DEG>>::default();
		let mut b_coeff_iter = b_poly.coeff_iter_mut();
		for &Deg(Exp(coeff)) in poly.coeff_iter() {
			*b_coeff_iter.next().unwrap() = Deg(coeff);
		}
		b_poly.0 = poly.constant().0;
		Roots::roots(b_poly)
	}
}

impl LinearIso<f64> for Exp<f64> {
	fn map(self) -> f64 {
		self.0.exp()
	}
	fn inv_map(value: f64) -> Self {
		Self(value.ln())
	}
}

impl LinearIso<u64> for Exp<f64> {
	fn map(self) -> u64 {
		self.0.exp().round() as u64
	}
	fn inv_map(value: u64) -> Self {
		Self((value as f64).ln())
	}
}