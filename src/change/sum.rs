//! Summation over time.

use crate::change::{ApplyChange, Change, PopChange};
use crate::linear::{Basis, Linear};
use crate::poly::Poly;

/// The pattern of no change.
pub type Nil = ();

/// The pattern of repeated addition.
///
/// Represents the addition of `T` (a [`Basis`] type) per unit of time, which
/// itself may have a rate of change in `C` (a [`Change`] type).
///
/// - `(Times<A>, Plus<B>)` ~ `(_ * a) + b`
/// - `Plus<A, Times<B>>` ~ `_ + (a * b)`
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Plus<T, C=Nil> {
	pub basis: T,
	pub change: C,
}

impl<T, C> Change<T> for Plus<T, C>
where
	T: Basis,
	C: Change<T>,
	Nil: ChangeIntoPoly<Nil, Self, T>,
{
	type Poly = <Nil as ChangeIntoPoly<Nil, Self, T>>::Output;
	fn into_poly(self, basis: T) -> Self::Poly {
		().change_into_poly((), self, basis)
	}
	fn scale(self, scalar: <T as Basis>::Inner) -> Self {
		Self {
			basis: self.basis.map_inner(|n| n.mul(scalar)),
			change: self.change.scale(scalar),
		}
	}
}

impl<T, C> PopChange for Plus<T, C> {
	type Inner = Nil;
	type Change = Self;
	fn pop_change(self) -> (Self::Inner, Self::Change) {
		((), self)
	}
}

impl<T, C> ApplyChange<Nil> for Plus<T, C> {
	type Output = Self;
	fn apply_change(self, _change: Nil) -> Self::Output {
		self
	}
}

impl<T, A, B> ApplyChange<Plus<T, B>> for Plus<T, A>
where
	T: Basis,
	A: ApplyChange<B>,
{
	type Output = Plus<T, A::Output>;
	fn apply_change(self, change: Plus<T, B>) -> Self::Output {
		Plus {
			basis: self.basis.zip_map_inner(change.basis, Linear::add),
			change: self.change.apply_change(change.change),
		}
	}
}

impl<T, A, B> ApplyChange<Times<T, B>> for Plus<T, A> {
	type Output = (Self, Times<T, B>);
	fn apply_change(self, change: Times<T, B>) -> Self::Output {
		(self, change)
	}
}

/// The pattern of repeated multiplication.
///
/// Represents the multiplication of `T` (a [`Basis`] type) per unit of time,
/// which itself may have a rate of change in `C` (a [`Change`] type).
///
/// - `(Plus<A>, Times<B>)` ~ `(_ + a) * b`
/// - `Times<A, Plus<B>>` ~ `_ * (a + b)`
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Times<T, C=Nil> {
	pub basis: T,
	pub change: C,
}

impl<T, C> Change<T> for Times<T, C>
where
	T: Basis,
	C: Change<T>,
	Nil: ChangeIntoPoly<Nil, Self, T>,
{
	type Poly = <Nil as ChangeIntoPoly<Nil, Self, T>>::Output;
	fn into_poly(self, basis: T) -> Self::Poly {
		().change_into_poly((), self, basis)
	}
	fn scale(self, scalar: <T as Basis>::Inner) -> Self {
		Self {
			basis: self.basis.map_inner(|n| n.mul(scalar)),
			change: self.change.scale(scalar),
		}
	}
}

impl<T, C> PopChange for Times<T, C> {
	type Inner = Nil;
	type Change = Self;
	fn pop_change(self) -> (Self::Inner, Self::Change) {
		((), self)
	}
}

impl<T, C> ApplyChange<Nil> for Times<T, C> {
	type Output = Self;
	fn apply_change(self, _change: Nil) -> Self::Output {
		self
	}
}

impl<T, A, B> ApplyChange<Times<T, B>> for Times<T, A>
where
	T: Basis,
	A: ApplyChange<B>,
{
	type Output = Times<T, A::Output>;
	fn apply_change(self, change: Times<T, B>) -> Self::Output {
		Times {
			basis: self.basis.zip_map_inner(change.basis, Linear::mul),
			change: self.change.apply_change(change.change),
		}
	}
}

impl<T, A, B> ApplyChange<Plus<T, B>> for Times<T, A> {
	type Output = (Self, Plus<T, B>);
	fn apply_change(self, change: Plus<T, B>) -> Self::Output {
		(self, change)
	}
}

/// ...
pub trait ChangeIntoPoly<A, B, T: Basis> {
	type Output: Poly<T>;
	fn change_into_poly(self, first: A, second: B, basis: T) -> Self::Output;
}

mod _impl_change_into_poly {
	use crate::change::{Change, PopChange, Nil, Plus, Times};
	use crate::linear::Basis;
	use crate::poly::Poly;
	use std::ops::{Add, Mul, Div};
	use super::ChangeIntoPoly;
	use symb_poly::{Integ, ProdInteg, SymbolizeInner, Unsymbolize};
	
	// A + Integ!(C)
	impl<T, A, C, O, P> ChangeIntoPoly<Nil, Plus<T, C>, T> for A
	where
		A: Change<T, Poly: Add<P, Output=O>>,
		C: Change<T, Poly: Integ<Output=P>>,
		O: Poly<T>,
		T: Basis,
	{
		type Output = O;
		fn change_into_poly(self, _first: Nil, second: Plus<T, C>, basis: T) -> Self::Output {
			self.into_poly(basis) + second.change.into_poly(second.basis).integ()
		}
	}
	
	// A * ProdInteg!(C)
	impl<T, A, C, O, P> ChangeIntoPoly<Nil, Times<T, C>, T> for A
	where
		A: Change<T, Poly: Mul<P, Output=O>>,
		C: Change<T, Poly: ProdInteg<Output=P>>,
		O: Poly<T>,
		T: Basis,
	{
		type Output = O;
		fn change_into_poly(self, _first: Nil, second: Times<T, C>, basis: T) -> Self::Output {
			self.into_poly(basis) * second.change.into_poly(second.basis).prod_integ()
		}
	}
	
	// (A + Integ!(B)) + Integ!(C)
	impl<T, A, B, C, O, P> ChangeIntoPoly<Plus<T, B>, Plus<T, C>, T> for A
	where
		A: PopChange<Inner: ChangeIntoPoly<A::Change, Plus<T, B>, T, Output: Add<P::Output, Output=O>>>,
		C: Change<T, Poly=P>,
		O: Poly<T>,
		P: Integ,
		T: Basis,
	{
		type Output = O;
		fn change_into_poly(self, first: Plus<T, B>, second: Plus<T, C>, basis: T) -> Self::Output {
			let (a_inner, a_outer) = self.pop_change();
			let p = second.change.into_poly(second.basis);
			a_inner.change_into_poly(a_outer, first, basis) + p.integ()
		}
	}
	
	// (A * ProdInteg!(B)) * ProdInteg!(C)
	impl<T, A, B, C, O, P> ChangeIntoPoly<Times<T, B>, Times<T, C>, T> for A
	where
		A: PopChange<Inner: ChangeIntoPoly<A::Change, Times<T, B>, T, Output: Mul<P::Output, Output=O>>>,
		C: Change<T, Poly=P>,
		O: Poly<T>,
		P: ProdInteg,
		T: Basis,
	{
		type Output = O;
		fn change_into_poly(self, first: Times<T, B>, second: Times<T, C>, basis: T) -> Self::Output {
			let (a_inner, a_outer) = self.pop_change();
			let p = second.change.into_poly(second.basis);
			a_inner.change_into_poly(a_outer, first, basis) * p.prod_integ()
		}
	}
	
	// (A * ProdInteg!(B)) + Integ!(C) -> (A + Integ!(C / ProdInteg!(B))) * ProdInteg!(B)
	impl<A, B, C, I, J, O, P, Q, R, T> ChangeIntoPoly<Times<T, B>, Plus<T, C>, T> for A
	where
		A: Change<T, Poly: SymbolizeInner<typenum::U0, Output=O, NextSymbol=I>>,
		B: Change<T, Poly: SymbolizeInner<I, Output=P, NextSymbol=J>>,
		C: Change<T, Poly: SymbolizeInner<J, Output=Q>>,
		O: Add<<P::Output as Integ>::Output, Output: Mul<Q::Output, Output=R>>,
		P: Div<Q::Output, Output: Integ>,
		Q: ProdInteg<Output: Clone>,
		R: Unsymbolize<Output: Poly<T>>,
		T: Basis,
	{
		type Output = R::Output;
		fn change_into_poly(self, b: Times<T, B>, c: Plus<T, C>, basis: T) -> Self::Output {
			let (o, i) = self.into_poly(basis).symbolize_inner(Default::default());
			let (p, j) = b.change.into_poly(b.basis).symbolize_inner(i);
			let (q, _) = c.change.into_poly(c.basis).symbolize_inner(j);
			let q_expr = q.prod_integ();
			let p_expr = (p / q_expr.clone()).integ();
			let o_expr = (o + p_expr) * q_expr;
			o_expr.unsymbolize()
		}
	}
	
	// (A + Integ!(B)) * ProdInteg!(C) -> (A + Integ!(B * C/ProdInteg!(C))) * ProdInteg!(C)
	impl<A, B, C, I, J, O, P, Q, R, T> ChangeIntoPoly<Plus<T, B>, Times<T, C>, T> for A
	where
		A: Change<T, Poly: SymbolizeInner<typenum::U0, Output=O, NextSymbol=I>>,
		B: Change<T, Poly: SymbolizeInner<I, Output=P, NextSymbol=J>>,
		C: Change<T, Poly: SymbolizeInner<J, Output=Q>>,
		O: Add<<<P::Output as Div<Q::Output>>::Output as Integ>::Output, Output: Mul<Q::Output, Output=R>>,
		P: Mul<Q, Output: Div<Q::Output, Output: Integ>>,
		Q: ProdInteg<Output: Clone> + Clone,
		R: Unsymbolize<Output: Poly<T>>,
		T: Basis,
	{
		type Output = R::Output;
		fn change_into_poly(self, b: Plus<T, B>, c: Times<T, C>, basis: T) -> Self::Output {
			let (o, i) = self.into_poly(basis).symbolize_inner(Default::default());
			let (p, j) = b.change.into_poly(b.basis).symbolize_inner(i);
			let (q, _) = c.change.into_poly(c.basis).symbolize_inner(j);
			let q_expr = q.clone().prod_integ();
			let p_expr = (p * q / q_expr.clone()).integ();
			let o_expr = (o + p_expr) * q_expr;
			o_expr.unsymbolize()
		}
	}
}

// #[cfg(feature = "glam")]
// impl<const D: usize> Roots for SumPoly<glam::DVec2, D>
// where
// 	SumPoly<f64, D>: Poly<Basis=f64> + Roots,
// 	<SumPoly<f64, D> as Roots>::Output: IntoIterator<Item=f64>,
// {
// 	type Output = [f64; D];
// 	fn roots(self) -> <Self as Roots>::Output {
// 		let mut root_list = [f64::NAN; D];
// 		let mut root_count = 0;
//
// 		let mut x_poly = SumPoly::zero();
// 		let mut y_poly = SumPoly::zero();
// 		for (index, coeff) in self.into_iter().enumerate() {
// 			x_poly[index] = coeff.x;
// 			y_poly[index] = coeff.y;
// 		}
//
// 		let mut x_roots = x_poly.roots().into_iter().collect::<Vec<_>>();
// 		let mut y_roots = y_poly.roots().into_iter().collect::<Vec<_>>();
// 		x_roots.sort_unstable_by(f64::total_cmp);
// 		y_roots.sort_unstable_by(f64::total_cmp);
//
// 		let mut x_iter = x_roots.into_iter();
// 		let mut y_iter = y_roots.into_iter();
// 		let mut x = x_iter.next();
// 		let mut y = y_iter.next();
// 		while let (Some(a), Some(b)) = (x.as_ref(), y.as_ref()) {
// 			match a.total_cmp(b) {
// 				Ordering::Less    => x = x_iter.next(),
// 				Ordering::Greater => y = y_iter.next(),
// 				Ordering::Equal   => {
// 					root_list[root_count] = *a;
// 					root_count += 1;
// 					x = x_iter.next();
// 					y = y_iter.next();
// 					// !!! What if the next value(s) are equal to the current?
// 					// Should only one of these be advanced?
// 				}
// 			}
// 		}
//
// 		root_list
// 	}
// }

#[cfg(test)]
mod tests {
	use crate::time::*;
	use crate::temporal::Temporal;
	use crate::pred::Prediction;
	use crate::*;

	#[cfg(feature = "glam")]
	#[test]
	fn vec() {
		use crate as chime;
		
		#[derive(PartialEq, Flux, ToMoment, ToMomentMut)]
		struct Pos {
			#[basis]
			value: Iso<glam::DVec2, glam::IVec2>,
			#[change(add_per(SEC))]
			spd: Spd,
		}
		
		#[derive(PartialEq, Flux, ToMoment, ToMomentMut)]
		struct Spd {
			#[basis]
			value: Iso<glam::DVec2, glam::IVec2>
		}
		
		let a_pos = Temporal::from(Pos {
			value: glam::IVec2::new(10, 10).into(),
			spd: Spd {
				value: glam::IVec2::new(6, 4).into()
			}
		});
		
		let b_pos = Temporal::from(Pos {
			value: glam::IVec2::new(14, 18).into(),
			spd: Spd {
				value: glam::IVec2::new(4, 0).into()
			}
		});
		
		// println!("{:?}", a_pos.poly() - b_pos.poly());

		todo!("re-add Roots for glam::* impls")
		// assert_eq!(
		// 	Vec::from_iter(a_pos.to_poly().when_eq(b_pos.to_poly())
		// 		.into_ranges(Time::ZERO)
		// 		.inclusive()),
		// 	[(2*SEC - 83333333*NANOSEC, 2*SEC + 83333333*NANOSEC)]
		// );
	}
	
	#[test]
	fn precise() {
		// https://www.desmos.com/calculator/1z97cqlopx
		
		 // Basic Check:
		let a = Temporal::new(poly!(|x| -193.99999999999997 + 4.481238217799146e-6*x - 500.*x^2), SEC);
		let b = Temporal::new(poly!(|_| -194.), SEC);
		assert_eq!(
			Vec::from_iter(a.when_eq(b)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(SEC-5*NANOSEC, SEC-3*NANOSEC),
				(SEC+12*NANOSEC, SEC+14*NANOSEC)
			]
		);
		
		 // Distance Check:
		let a = Temporal::new([
			poly!(|x| 0.0036784761334161292 + 1.1687626970174242e-7*x + 0.*x^2),
			poly!(|x| -182.00000057575835 - 7.537214753878195e-7*x - 500.*x^2)
		], SEC);
		let b = Temporal::new([
			poly!(|_| -3.8808053943969956e-5),
			poly!(|_| -193.99999999999997)
		], SEC);
		let dis = Temporal::new(poly!(|_| 12.), SEC);
		assert_eq!(
			Vec::from_iter(a.when_dis_eq(b, dis)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(780910981*NANOSEC, 780910981*NANOSEC),
				(SEC-13*NANOSEC, SEC-10*NANOSEC),
				(SEC+2*NANOSEC, SEC+8*NANOSEC),
				(1219089016*NANOSEC, 1219089016*NANOSEC),
			]
		);
	}
}