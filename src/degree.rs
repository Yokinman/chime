//! Compile-time assured degree of change over time.
//! 
//! - `Deg<0>`: Constant
//! - `Deg<1>`: Linear
//! - `Deg<2>`: Quadratic
//! - `Deg<3>`: Cubic
//! - etc.

use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::ops::{Index, IndexMut};

mod private {
	pub trait Sealed {}
}

/// Degree type.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Deg<const DEG: usize>;
impl<const DEG: usize> private::Sealed for Deg<DEG> {}

/// Degree category.
pub trait IsDeg: private::Sealed + Copy + Clone + Debug {
	// ??? Could name this "Degree", but that seems confusing and ambiguous with
	// the `Deg` type. Similar to `Ord` vs `Ordering` or `Iterator` vs `Iter`. 
	
	const USIZE: usize;
	
	type Array<T: Copy + Clone + Debug>:
		Copy + Clone + Debug
		+ IntoIterator<Item=T>
		+ AsRef<[T]> + AsMut<[T]>
		+ Borrow<[T]> + BorrowMut<[T]>
		+ Index<usize, Output=T> + IndexMut<usize>;
	
	fn array_from<T: Copy + Clone + Debug>(value: T) -> Self::Array<T>;
}
impl<const DEG: usize> IsDeg for Deg<DEG> {
	const USIZE: usize = DEG;
	
	type Array<T: Copy + Clone + Debug> = [T; DEG];
	
	fn array_from<T: Copy + Clone + Debug>(value: T) -> Self::Array<T> {
		[value; DEG]
	}
}

/// Degree increment.
pub trait HasUpDeg: IsDeg {
	type Up: IsDeg;
}

/// Degree decrement.
pub trait HasDownDeg: IsDeg {
	type Down: IsDeg;
}

/// Degree sequential ordering.
macro_rules! impl_deg_order {
	(1  1  $($num:tt)*) => { impl_deg_order!(2  $($num)*); };
	(2  2  $($num:tt)*) => { impl_deg_order!(4  $($num)*); };
	(4  4  $($num:tt)*) => { impl_deg_order!(8  $($num)*); };
	(8  8  $($num:tt)*) => { impl_deg_order!(16 $($num)*); };
	(16 16 $($num:tt)*) => { impl_deg_order!(32 $($num)*); };
	(32 32 $($num:tt)*) => { impl_deg_order!(64 $($num)*); };
	(64) => {/* break */};
	($($num:tt)+) => {
		impl HasUpDeg for Deg<{ $($num +)+ 0 - 1 }> {
			type Up = Deg<{ $($num +)+ 0 }>;
		}
		impl HasDownDeg for Deg<{ $($num +)+ 0 }> {
			type Down = Deg<{ $($num +)+ 0 - 1 }>;
		}
		impl_deg_order!(1 $($num)+);
	};
}
impl_deg_order!(1);

/// Degree comparison.
/// 
/// # Examples
/// 
/// ```
/// use flux::degree::*;
/// fn three_below_seven() where Deg<3>: IsBelowDeg<Deg<7>> {}
/// fn fifteen_below_twenty_one() where Deg<15>: IsBelowDeg<Deg<21>> {}
/// ```
/// 
/// ```compile_fail
/// use flux::degree::*;
/// fn twenty_one_below_fifteen() where Deg<21>: IsBelowDeg<Deg<15>> {}
/// ```
pub trait IsBelowDeg<D: IsDeg>: IsDeg {}

 // 0 < D:
impl<D> IsBelowDeg<D> for Deg<0> where D: IsDeg + HasDownDeg {}

 // A < B, if A-1 < B-1: 
impl<A, B> IsBelowDeg<B> for A
where
	A: IsDeg + HasDownDeg,
	B: IsDeg + HasDownDeg,
	A::Down: IsBelowDeg<B::Down>, 
{}

/// Degree maximum.
/// 
/// # Examples
/// 
/// ```
/// use flux::degree::*;
/// use std::any::TypeId;
/// assert_eq!(
///     TypeId::of::<<Deg<3> as MaxDeg<Deg<7>>>::Max>(),
///     TypeId::of::<Deg<7>>()
/// );
/// assert_eq!(
///     TypeId::of::<<Deg<21> as MaxDeg<Deg<15>>>::Max>(),
///     TypeId::of::<Deg<21>>()
/// );
/// ```
pub trait MaxDeg<D: IsDeg, N: IsDeg = Deg<0>>: IsDeg {
	type Max: IsDeg;
}

 // max(0, 0) + N = N:
impl<N: IsDeg> MaxDeg<Deg<0>, N> for Deg<0> {
	type Max = N;
}

 // max(0, d) + N = max(0, d-1) + (N+1):
impl<D, N> MaxDeg<D, N> for Deg<0>
where
	D: IsDeg + HasDownDeg,
	N: IsDeg + HasUpDeg,
	Deg<0>: MaxDeg<D::Down, N::Up>,
{
	type Max = <Deg<0> as MaxDeg<D::Down, N::Up>>::Max;
}

 // max(D, 0) + N = max(D-1, 0) + (N+1):
impl<D, N> MaxDeg<Deg<0>, N> for D
where
	D: IsDeg + HasDownDeg,
	N: IsDeg + HasUpDeg,
	Deg<0>: MaxDeg<D::Down, N::Up>,
{
	type Max = <Deg<0> as MaxDeg<D::Down, N::Up>>::Max;
}

 // max(A, B) + N = max(A-1, B-1) + (N+1):
impl<A, B, N> MaxDeg<B, N> for A
where
	A: IsDeg + HasDownDeg,
	B: IsDeg + HasDownDeg,
	N: IsDeg + HasUpDeg,
    A::Down: MaxDeg<B::Down, N::Up>,
{
	type Max = <A::Down as MaxDeg<B::Down, N::Up>>::Max;
}

// /// Degree addition.
// trait AddDeg<D: IsDeg>: IsDeg {
// 	type Sum: IsDeg;
// }
// 
//  // D + 0 = D:
// impl<D: IsDeg + HasDownDeg> AddDeg<Deg<0>> for D {
// 	type Sum = D;
// }
// 
//  // A + B = (A+1) + (B-1):
// impl<A, B> AddDeg<B> for A
// where
// 	A: IsDeg + HasUpDeg,
// 	B: IsDeg + HasDownDeg,
// 	A::Up: AddDeg<B::Down>,
// {
// 	type Sum = <A::Up as AddDeg<B::Down>>::Sum;
// }