use std::cmp::Ordering;
use crate::linear::Basis;
use crate::time::Time;

pub use crate::change::constant::Constant;

/// An abstract description of change over time.
/// 
/// Used to define the standard representations of [`Flux`] types. In
/// other words, the layout of a polynomial.
pub trait Poly<B: Basis>: Clone {
	const DEGREE: usize;
	
	fn eval(&self, time: B::Inner) -> B;
	
	/// The order at or immediately preceding the value at a time.
	/// 
	/// This should be the first non-zero [`Poly::eval`] value of this kind
	/// or its derivatives; reversed for odd derivatives.
	fn initial_order(self, time: B::Inner) -> Option<Ordering>
	where
		Self: Deriv<B>,
		B: PartialOrd,
	{
		// !!! Alternative: Translate polynomial using `to_time` and then check
		// leading terms in order. Unknown which is more precise/faster.
		
		let order = self.eval(time)
			.partial_cmp(&Basis::zero());
		
		if order != Some(Ordering::Equal) || Self::DEGREE == 0 {
			return order
		}

		// Deriv::deriv(self).initial_order(time).map(Ordering::reverse)
		todo!()
	}
	
	fn is_zero(&self) -> bool;
}

mod _impl_poly {
	use crate::linear::Basis;
	use crate::poly::Constant;
	use super::Poly;
	use symb_poly::{Eval, InsertVal, Invar, InvarTerm, Limit, Func, Sum, Power, Prod, Var, Variable};
	
	impl<T, B, const SIZE: usize> Poly<[B; SIZE]> for [T; SIZE]
	where
		T: Poly<B>,
		B: Basis,
	{
		const DEGREE: usize = T::DEGREE;
		fn eval(&self, time: B::Inner) -> [B; SIZE] {
			self.each_ref().map(|x| T::eval(x, time))
		}
		fn is_zero(&self) -> bool {
			self.iter().all(T::is_zero)
		}
	}
	
	impl<T, B> Poly<B> for Invar<T>
	where
		T: InvarTerm + InvarPoly<B, T::Kind> + Clone,
		B: Basis,
	{
		const DEGREE: usize = 0;
		
		fn eval(&self, time: B::Inner) -> B {
			self.0.clone().invar_eval(time)
		}
		
		fn is_zero(&self) -> bool {
			InvarPoly::is_zero(&self.0)
		}
	}
	
	/// ...
	pub trait InvarPoly<B: Basis, Kind> {
		fn invar_eval(self, _time: B::Inner) -> B;
		fn is_zero(&self) -> bool;
	}
	
	mod _impl_invar_poly {
		use crate::linear::Basis;
		use crate::poly::Constant;
		use super::{InvarPoly, Poly};
		use symb_poly::{ExprKind, Finite, FromInvar, Invar, Known, NegInf, NumKind, Unknown, Unsymbolize, Zero};
		
		impl<T, B> InvarPoly<B, NumKind<Zero>> for T
		where
			B: Basis,
			Constant<B>: FromInvar<T>,
		{
			fn invar_eval(self, _time: B::Inner) -> B {
				<Constant<B>>::from_invar(self).0
			}
			fn is_zero(&self) -> bool {
				true
			}
		}
		
		impl<T, B> InvarPoly<B, NumKind<NegInf>> for T
		where
			B: Basis,
			Constant<B>: FromInvar<T>,
		{
			fn invar_eval(self, _time: B::Inner) -> B {
				<Constant<B>>::from_invar(self).0
			}
			fn is_zero(&self) -> bool {
				false
			}
		}
		
		impl<T, B, N> InvarPoly<B, NumKind<Finite<Known<N>>>> for T
		where
			B: Basis,
			Constant<B>: FromInvar<T>,
		{
			fn invar_eval(self, _time: B::Inner) -> B {
				<Constant<B>>::from_invar(self).0
			}
			fn is_zero(&self) -> bool {
				false
			}
		}
		
		impl<B> InvarPoly<B, NumKind<Finite<Unknown>>> for Constant<B>
		where
			B: Basis,
			// Constant<B>: FromInvar<T>,
		{
			fn invar_eval(self, _time: B::Inner) -> B {
				/*<Constant<B>>::from_invar*/(self).0
			}
			fn is_zero(&self) -> bool {
				self.0 == B::zero()
			}
		}
		
		impl<T, B> InvarPoly<B, ExprKind> for T
		where
			Invar<T>: Unsymbolize<Output: Poly<B>>,
			T: Clone,
			B: Basis,
		{
			fn invar_eval(self, time: B::Inner) -> B {
				Invar(self).unsymbolize().eval(time)
			}
			fn is_zero(&self) -> bool {
				Invar(self.clone()).unsymbolize().is_zero()
			}
		}
	}
	
	impl<T> Poly<T> for Var
	where
		T: Basis
	{
		const DEGREE: usize = 0;
		fn eval(&self, time: T::Inner) -> T {
			T::from_inner(time)
		}
		fn is_zero(&self) -> bool {
			false
		}
	}
	
	impl<A, B, T> Poly<T> for Func<Sum, (A, B)>
	where
		A: Poly<T>,
		B: Poly<T>,
		T: Basis,
		Self: InsertVal<Variable, Constant<T>, Output: Eval<Output=Constant<T>>>,
	{
		const DEGREE: usize = 0;
		fn eval(&self, time: T::Inner) -> T {
			self.clone()
				.insert_val(<Variable>::default(), Constant(T::from_inner(time)))
				.eval().0
		}
		fn is_zero(&self) -> bool {
			let (a, b) = self.args();
			a.is_zero() && b.is_zero()
		}
	}
	
	impl<A, B, T> Poly<T> for Func<Prod, (A, B)>
	where
		A: Poly<T>,
		B: Poly<T>,
		T: Basis,
		Self: InsertVal<Variable, Constant<T>, Output: Eval<Output=Constant<T>>>,
	{
		const DEGREE: usize = 0;
		fn eval(&self, time: T::Inner) -> T {
			self.clone()
				.insert_val(<Variable>::default(), Constant(T::from_inner(time)))
				.eval().0
		}
		fn is_zero(&self) -> bool {
			let (a, _b) = self.args();
			a.is_zero()// || b.is_zero()
		}
	}
	
	impl<A, B, T> Poly<T> for Func<Power, (A, B)>
	where
		A: Poly<T>,
		B: Poly<T>,
		T: Basis,
		Self: InsertVal<Variable, Constant<T>, Output: Eval<Output=Constant<T>>>,
	{
		const DEGREE: usize = 0;
		
		fn eval(&self, time: T::Inner) -> T {
			self.clone()
				.insert_val(<Variable>::default(), Constant(T::from_inner(time)))
				.eval().0
		}
		
		fn is_zero(&self) -> bool {
			let (a, _b) = self.args();
			a.is_zero()// && !b.is_zero()
		}
	}
	
	impl<A, I, R, T> Poly<T> for Limit<A, I, R>
	where
		T: Basis,
		Self: Clone + InsertVal<Variable, Constant<T>, Output: Eval<Output=Constant<T>>>,
	{
		const DEGREE: usize = 0;
		
		fn eval(&self, time: T::Inner) -> T {
			self.clone()
				.insert_val(<Variable>::default(), Constant(T::from_inner(time)))
				.eval().0
		}
		
		fn is_zero(&self) -> bool {
			// self.expr.is_zero()
			false
		}
	}
}

impl<T, B, const N: usize> Translate<[B; N]> for [T; N]
where
	T: Translate<B>,
	B: Basis,
{
	type Output = [T::Output; N];
	fn translate(self, amount: B::Inner) -> Self::Output {
		self.map(|x| x.translate(amount))
	}
}

impl<T, B, const N: usize> Deriv<[B; N]> for [T; N]
where
	T: Deriv<B>,
	B: Basis,
{
	type Deriv = [T::Deriv; N];
	fn deriv(self) -> Self::Deriv {
		self.map(T::deriv)
	}
}

/// ...
pub trait Deriv<B: Basis>: Poly<B> {
	type Deriv: Poly<B>;
	fn deriv(self) -> Self::Deriv;
}

mod _impl_deriv {
	use crate::linear::Basis;
	use super::{Deriv, Poly};

	impl<T, B> Deriv<B> for symb_poly::Invar<T>
	where
		Self: Poly<B> + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<B>>,
		B: Basis,
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
	
	impl<T, B> Deriv<B> for symb_poly::Var<T>
	where
		Self: Poly<B> + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<B>>,
		B: Basis,
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
	
	impl<F, A, B> Deriv<B> for symb_poly::Func<F, A>
	where
		Self: Poly<B> + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<B>>,
		B: Basis,
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
	
	impl<T, I, R, B> Deriv<B> for symb_poly::Limit<T, I, R>
	where
		Self: Poly<B> + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<B>>,
		B: Basis,
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
}

/// ...
pub trait Translate<B: Basis>: Poly<B> {
	type Output: Poly<B>;
	fn translate(self, amount: B::Inner) -> Self::Output;
}

mod _impl_translate {
	use super::Translate;
	use symb_poly::{Func, Sum, Invar, Limit, Var, Replace, Variable};
	use crate::linear::Basis;
	use crate::constant::Constant;
	use crate::poly::Poly;
	
	impl<A, B> Translate<B> for Invar<A>
	where
		Self: Poly<B>,
		B: Basis,
	{
		type Output = Self;
		fn translate(self, _amount: B::Inner) -> Self::Output {
			self
		}
	}
	
	impl<F, A, B> Translate<B> for Func<F, A>
	where
		Self: Poly<B> + Replace<Variable, Func<Sum, (Invar<Constant<B>>, Var)>, Output: Poly<B>>,
		B: Basis,
	{
		type Output = <Self as Replace<Variable, Func<Sum, (Invar<Constant<B>>, Var)>>>::Output;
		fn translate(self, amount: B::Inner) -> Self::Output {
			self.replace(Func::from_args((Invar(Constant(Basis::from_inner(amount))), Default::default())))
		}
	}
	
	impl<T, I, R, B> Translate<B> for Limit<T, I, R>
	where
		Self: Poly<B>,
		T: Replace<Variable, Func<Sum, (Invar<Constant<B>>, Var)>>,
		R: Replace<Variable, Func<Sum, (Invar<Constant<B>>, Var)>>,
		B: Basis,
		Limit<T::Output, I, R::Output>: Poly<B>,
	{
		type Output = Limit<T::Output, I, R::Output>;
		fn translate(self, amount: B::Inner) -> Self::Output {
			Limit {
				expr: self.expr.replace(Func::from_args((Invar(Constant(Basis::from_inner(amount))), Default::default()))),
				index: self.index,
				result: self.result.replace(Func::from_args((Invar(Constant(Basis::from_inner(amount))), Default::default()))),
			}
		}
	}
}

pub(crate) mod ops {
	use std::ops;
	
	/// Squaring a kind of change.
	pub trait Sqr {
		type Output;
		fn sqr(self) -> Self::Output;
	}
	
	impl<K> Sqr for K
	where
		K: ops::Mul + Clone
	{
		type Output = K::Output;
		fn sqr(self) -> Self::Output {
			self.clone() * self
		}
	}
}

/// Roots of a [`Poly`] polynomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots<B: Basis>: Poly<B> {
	type Output: IntoTimes;
	fn roots(self) -> <Self as Roots<B>>::Output;
}

mod _impl_roots {
	use crate::constant::Constant;
	use crate::poly::Poly;
	use std::cmp::Ordering;
	use num_traits::Pow;
	use super::Roots;
	use symb_poly::{Func, Invar, Num, Power, Prod, Sum, Var};
	use typenum::{P2, P3, P4};

	impl Roots<f64> for Invar<Constant<f64>> {
		type Output = [f64; 0];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			[0.; 0]
		}
	}

	impl Roots<f64> for Func<Prod, (Invar<Constant<f64>>, Var)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), _) = self.into_args();
			[0. / a; 1]
		}
	}

	impl Roots<f64> for Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), _) = self.into_args();
			[0. / a; 2]
		}
	}

	impl Roots<f64> for Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)> {
		type Output = [f64; 3];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), _) = self.into_args();
			[0. / a; 3]
		}
	}

	impl Roots<f64> for Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)> {
		type Output = [f64; 4];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), _) = self.into_args();
			[0. / a; 4]
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Prod, (Invar<Constant<f64>>, Var)>
	)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_prod) = self.into_args();
			let (Invar(Constant(b)), _) = b_prod.into_args();
			[-a / b; 1]
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_prod) = self.into_args();

			if a == 0. {
				return b_prod.roots()
			}

			let (Invar(Constant(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 2];
			roots[0] = (-a / b).powf((2 as f64).recip());
			roots[1] = -roots[0];
			roots
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_prod) = self.into_args();

			if a == 0. {
				return b_prod.roots()
			}

			let (Invar(Constant(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 3];
			roots[0] = (-a / b).powf((3 as f64).recip());
			roots
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_prod) = self.into_args();

			if a == 0. {
				return b_prod.roots()
			}

			let (Invar(Constant(b)), _) = b_prod.into_args();
			let mut roots = [f64::NAN; 4];
			roots[0] = (-a / b).powf((4 as f64).recip());
			roots[1] = -roots[0];
			roots
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>
		)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, c_prod) = b_sum.into_args();
			let (Invar(Constant(b)), _) = b_prod.into_args();
			let (Invar(Constant(c)), _) = c_prod.into_args();

			let mut n = -b / c;

			// Precision Breakdown (Discriminant Quotient):
			let x = -a / b;
			if x.is_nan() || n.is_nan() || (x / n).abs() < 1e-8 {
				return if x.is_nan() {
					[n, n]
				} else {
					[x, n-x]
				}
			}

			// General Quadratic:
			n /= 2.;
			let x = n.mul_add(n, -a/c).sqrt();
			[n-x, n+x]
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, d_prod) = b_sum.into_args();
			let (Invar(Constant(b)), _) = b_prod.into_args();
			let (Invar(Constant(d)), d_var) = d_prod.into_args();

			// Pseudo-linear:
			if b == 0. {
				return (Invar(Constant(a)) + Invar(Constant(d))*d_var).roots()
			}

			// Depressed Cubic:
			let p = -a / (2. * d);
			let q = -b / (3. * d);
			let r = p / q;
			let discriminant = if q < 0. {
				1.
			} else {
				let disc = f64::ln((r*r) / q.abs());
				if disc.abs() < 1e-15 {
					0.
				} else {
					disc
				}
			};
			match discriminant.partial_cmp(&0.) {
				// 3 Real Roots:
				Some(Ordering::Less) => {
					let sqrt_q = q.sqrt();
					let angle = f64::acos(r / sqrt_q);
					debug_assert!(!angle.is_nan(), "{:?}", self);
					use std::f64::consts::TAU;
					[
						2. * sqrt_q * f64::cos((angle         ) / 3.),
						2. * sqrt_q * f64::cos((angle +    TAU) / 3.),
						2. * sqrt_q * f64::cos((angle + 2.*TAU) / 3.),
					]
				},

				// 1 Real Root:
				Some(Ordering::Greater) => {
					let n = q.mul_add(-q*q, p*p).sqrt();
					debug_assert!(!n.is_nan());
					[(p + n).cbrt() + (p - n).cbrt(), f64::NAN, f64::NAN]
				},

				// 2 Real Roots:
				_ => [-r, -r, 2.*r]
			}
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			(self + (Invar(Constant(0.)) * <Var>::default())).roots()
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 3];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (c_prod, d_prod) = c_sum.into_args();
			let (Invar(Constant(b)), b_var) = b_prod.into_args();
			let (Invar(Constant(c)), c_var) = c_prod.into_args();
			let (Invar(Constant(d)), d_var) = d_prod.into_args();

			// Weak Constant:
			let x = -a / b;
			if x.is_nan() || (
				((x   * c) / b).abs() < 1e-5 && // ??? Adjust as needed
					((x*x * d) / b).abs() < 1e-16
			) {
				let [y, z] = (Invar(Constant(b)) + Invar(Constant(c))*b_var + Invar(Constant(d))*c_var).roots();
				return if x.is_nan() {
					[y, y, z]
				} else {
					[x, y-x, z]
				}
			}

			let mut n = -c / d;

			// Weak Leading Term:
			if n.is_nan() || (
				(b / (n   * c)).abs() < 1e-9 && // ??? Adjust as needed
					(a / (n*n * c)).abs() < 1e-13
			) {
				let [x, y] = (Invar(Constant(a)) + Invar(Constant(b))*b_var + Invar(Constant(c))*c_var).roots();
				return [x, y, n]
			}

			// Depressed Cubic:
			if c == 0. {
				return (Invar(Constant(a)) + Invar(Constant(b))*b_var + Invar(Constant(d))*d_var).roots()
			}

			// General Cubic:
			n /= 3.;
			let depressed_cubic = Invar(Constant(n.mul_add(n.mul_add(c * (2./3.), b), a)))
				+ Invar(Constant(n.mul_add(c, b)))*b_var
				+ Invar(Constant(d))*d_var;
			depressed_cubic.roots().map(|x| x + n)
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self
				+ Invar(Constant(0.))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant(0.))*<Var>::default().pow(Invar(Num::<P3>::default())))
				.roots()
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), c_sum) = self.into_args();
			let (c_prod, e_prod) = c_sum.into_args();
			let (Invar(Constant(c)), _) = c_prod.into_args();

			// Pseudo-linear:
			if c == 0. {
				return (Invar(Constant(a)) + e_prod).roots()
			}

			let (Invar(Constant(e)), _) = e_prod.into_args();

			// Biquadratic:
			let [x, y] = (Invar(Constant(a))
				+ Invar(Constant(c))*<Var>::default()
				+ Invar(Constant(e))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			let (x, y) = (x.sqrt(), y.sqrt());
			return [-x, x, -y, y];
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self
				+ Invar(Constant(0.))
				+ Invar(Constant(0.))*<Var>::default().pow(Invar(Num::<P2>::default())))
				.roots()
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (Invar(Constant(b)), _) = b_prod.into_args();

			// Biquadratic:
			if b == 0. {
				return (Invar(Constant(a)) + c_sum).roots()
			}

			let (c_prod, e_prod) = c_sum.into_args();
			let (Invar(Constant(c)), _) = c_prod.into_args();
			let (Invar(Constant(e)), _) = e_prod.into_args();

			// Depressed Quartic:
			let p = a / e;
			let q = b / (2. * e);
			let r = c / (2. * e);
			let resolvent_b = r.mul_add(r, -p);
			let resolvent_cubic = Invar(Constant(-(q * q) / 2.))
				+ Invar(Constant(resolvent_b))*<Var>::default()
				+ Invar(Constant(2. * r))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant(1.))*<Var>::default().pow(Invar(Num::<P3>::default()));
			let mut m = resolvent_cubic.roots().into_iter()
				.find(|&r| r > 0. && r.is_finite())
				.unwrap_or_else(|| panic!(
					"expected a positive root from: {:?}",
					resolvent_cubic
				));
			let y = resolvent_cubic.eval(m)
				/ m.mul_add(m.mul_add(3., 4.*r), resolvent_b);
			if m > y && y.is_finite() {
				m -= y; // Newton-Raphson step
			}
			let sqrt_2m = (2. * m).sqrt();
			let [x, y] = (Invar(Constant((q / sqrt_2m) + r + m))
				+ Invar(Constant(-sqrt_2m))*<Var>::default()
				+ Invar(Constant(1.))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			let [z, w] = (Invar(Constant(-(q / sqrt_2m) + r + m))
				+ Invar(Constant(sqrt_2m))*<Var>::default()
				+ Invar(Constant(1.))*<Var>::default().pow(Invar(Num::<P2>::default()))
			).roots();
			[x, y, z, w]
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self + Invar(Constant(0.))*<Var>::default().pow(Invar(Num::<P2>::default()))).roots()
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			(self + Invar(Constant(0.))*<Var>::default()).roots()
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Sum, (
				Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P2>>)>)>,
				Func<Sum, (
					Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P3>>)>)>,
					Func<Prod, (Invar<Constant<f64>>, Func<Power, (Var, Invar<Num<P4>>)>)>
				)>
			)>
		)>
	)> {
		type Output = [f64; 4];
		fn roots(self) -> Self::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, c_sum) = b_sum.into_args();
			let (c_prod, d_sum) = c_sum.into_args();
			let (d_prod, e_prod) = d_sum.into_args();
			let (Invar(Constant(b)), b_var) = b_prod.into_args();
			let (Invar(Constant(c)), c_var) = c_prod.into_args();
			let (Invar(Constant(d)), d_var) = d_prod.into_args();
			let (Invar(Constant(e)), e_var) = e_prod.into_args();

			// Weak Constant:
			let x = -a / b;
			if x.is_nan() || (
				((x     * c) / b).abs() < 1e-4  && // ??? Adjust as needed
					((x*x   * d) / b).abs() < 1e-12 &&
					((x*x*x * e) / b).abs() < 1e-20
			) {
				let [y, z, w] = (Invar(Constant(b)) + Invar(Constant(c))*b_var + Invar(Constant(d))*c_var + Invar(Constant(e))*d_var).roots();
				return if x.is_nan() {
					[y, y, z, w]
				} else {
					[x, y-x, z, w]
				}
			}

			let mut n = -d / e;

			// Weak Leading Term:
			if n.is_nan() || (
				(c / (n     * d)).abs() < 1e-7  && // ??? Adjust as needed
					(b / (n*n   * d)).abs() < 1e-10 &&
					(a / (n*n*n * d)).abs() < 1e-15
			) {
				let [x, y, z] = (Invar(Constant(a)) + Invar(Constant(b))*b_var + Invar(Constant(c))*c_var + Invar(Constant(d))*d_var).roots();
				return [x, y, z, n]
			}

			// Depressed Quartic:
			if d == 0. {
				return (Invar(Constant(a)) + Invar(Constant(b))*b_var + Invar(Constant(c))*c_var + Invar(Constant(e))*e_var).roots()
			}

			// General Quartic:
			n /= 4.;
			let depressed_quartic = Invar(Constant(n.mul_add(n.mul_add(n.mul_add(d * (3./4.), c), b), a)))
				+ Invar(Constant(n.mul_add(n.mul_add(d, c) * 2., b)))*<Var>::default()
				+ Invar(Constant(n.mul_add(d * (3./2.), c)))*<Var>::default().pow(Invar(Num::<P2>::default()))
				+ Invar(Constant(e))*<Var>::default().pow(Invar(Num::<P4>::default()));
			let [x, y, z, w] = depressed_quartic.roots();
			[x+n, y+n, z+n, w+n]
		}
	}

	impl<A, B, const M: usize, const N: usize> Roots<f64> for Func<Sum, (Func<Prod, (Invar<Constant<f64>>, A)>, B)>
	where
		Self: Poly<f64> + std::ops::Div<Var, Output: Roots<f64, Output = [f64; M]>>,
		B: Roots<f64, Output = [f64; N]>,
	{
		type Output = [f64; N];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let mut roots = [0.; N];
			roots[1..=M].copy_from_slice(&(self / Var::default()).roots());
			roots
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Prod, (
			Invar<Constant<f64>>,
			Func<Power, (Invar<Constant<f64>>, Var)>
		)>
	)> {
		type Output = [f64; 1];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_prod) = self.into_args();
			let (Invar(Constant(b)), c_power) = b_prod.into_args();
			let (Invar(Constant(c)), _) = c_power.into_args();
			if c == 1. {
				return [1. - a/b]
			}
			[(-a/b).log(c)]
		}
	}

	impl Roots<f64> for Func<Sum, (
		Invar<Constant<f64>>,
		Func<Sum, (
			Func<Prod, (Invar<Constant<f64>>, Var)>,
			Func<Prod, (
				Invar<Constant<f64>>,
				Func<Power, (Invar<Constant<f64>>, Var)>
			)>
		)>
	)> {
		type Output = [f64; 2];
		fn roots(self) -> <Self as Roots<f64>>::Output {
			let (Invar(Constant(a)), b_sum) = self.into_args();
			let (b_prod, c_prod) = b_sum.into_args();
			let (Invar(Constant(b)), _) = b_prod.into_args();
			let (Invar(Constant(c)), d_power) = c_prod.into_args();
			let (Invar(Constant(d)), _) = d_power.into_args();

			if d == 1. {
				todo!()
			}
			// Convert to `a + b*x + c^x` and solve using product log.
			let x = a / c;
			let y = b / c;
			let d_ln = d.ln();
			let z = d_ln / (y * d.powf(x / y));
			[
				-(lambert_w::lambert_wm1(z) / d_ln) - (x / y),
				-(lambert_w::lambert_w0(z)  / d_ln) - (x / y),
			]
		}
	}
	
	impl<T, I, R, B> Roots<B> for symb_poly::Limit<T, I, R>
	where
		Self: Poly<B>,
		B: crate::linear::Basis,
		R: Roots<B>,
	{
		type Output = R::Output;
		fn roots(self) -> Self::Output {
			self.result.roots()
		}
	}
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

/// ...
pub struct Exponent<const N: usize>;

impl<T, const N: usize> std::ops::BitXor<Exponent<N>> for symb_poly::Invar<T>
where
	typenum::Const<N>: typenum::ToUInt,
	typenum::U<N>: typenum::Unsigned + typenum::NonZero,
	symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>: Default,
	Self: num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>,
{
	type Output = <Self as num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>>::Output;
	fn bitxor(self, _rhs: Exponent<N>) -> Self::Output {
		num_traits::Pow::pow(self, Default::default())
	}
}

impl<T, const N: usize> std::ops::BitXor<Exponent<N>> for symb_poly::Var<T>
where
	typenum::Const<N>: typenum::ToUInt,
	typenum::U<N>: typenum::Unsigned + typenum::NonZero,
	symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>: Default,
	Self: num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>,
{
	type Output = <Self as num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>>::Output;
	fn bitxor(self, _rhs: Exponent<N>) -> Self::Output {
		num_traits::Pow::pow(self, Default::default())
	}
}

impl<F, A, const N: usize> std::ops::BitXor<Exponent<N>> for symb_poly::Func<F, A>
where
	typenum::Const<N>: typenum::ToUInt,
	typenum::U<N>: typenum::Unsigned + typenum::NonZero,
	symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>: Default,
	Self: num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>,
{
	type Output = <Self as num_traits::Pow<symb_poly::Invar<symb_poly::Num<typenum::PInt<typenum::U<N>>>>>>::Output;
	fn bitxor(self, _rhs: Exponent<N>) -> Self::Output {
		num_traits::Pow::pow(self, Default::default())
	}
}

/// ...
pub trait IntoTerm {
	type Term;
	fn into_term(self) -> Self::Term;
}

impl<T> IntoTerm for T
where
	T: Basis
{
	type Term = symb_poly::Invar<Constant<T>>;
	fn into_term(self) -> Self::Term {
		symb_poly::Invar(Constant(self))
	}
}

impl<T> IntoTerm for symb_poly::Invar<T> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

impl<T> IntoTerm for symb_poly::Var<T> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

impl<F, A> IntoTerm for symb_poly::Func<F, A> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

/// ...
#[macro_export]
macro_rules! poly {
	 // Callback Handling:
	(@ return ($($t:tt)*) -> {
		$call:ident
		$(($([$($lhs:tt)*]:)? in $(:[$($rhs:tt)*])? ! $(, $args:tt)*))?
		$(-> $($rest:tt)+)?
	}) => {
		$crate::poly!{@ $call($($($($lhs)*)?)? $($t)* $($($($rhs)*)?$(, $args)*)?) $(-> $($rest)+)?}
	};
	(@ $rule:ident $input:tt -> $call:ident $($rest:tt)*) => {
		$crate::poly!{@ $rule $input -> { $call $($rest)* } }
	};
	
	 // Parse Terms:
	(@ parse_expr([- $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_expr([$($input)*], [1])
			-> parse_binop([-]:in!, $weight)
			-> $ret }
	};
	(@ parse_expr([$x:ident $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_binop($crate::poly::IntoTerm::into_term($x), [$($input)*], $weight) -> $ret }
	};
	(@ parse_expr([$x:literal $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_binop($crate::poly::IntoTerm::into_term($x), [$($input)*], $weight) -> $ret }
	};
	(@ parse_expr([($($x:tt)*) $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_expr([$($x)*], [])
			-> take_first
			-> parse_binop(in!, [$($input)*], $weight)
			-> $ret }
	};
	(@ parse_expr([], $weight:tt) -> $ret:tt) => {
		compile_error!("expected expression, found end of input")
	};
	
	 // Parse Exponents:
	(@ parse_const_expr([$x:ident $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_binop(($crate::poly::Exponent::<$x>), [$($input)*], $weight) -> $ret }
	};
	(@ parse_const_expr([$x:literal $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_binop(($crate::poly::Exponent::<$x>), [$($input)*], $weight) -> $ret }
	};
	(@ parse_const_expr([{$($x:tt)*} $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ parse_binop(($crate::poly::Exponent::<{$($x)*}>), [$($input)*], $weight) -> $ret }
	};
	(@ parse_const_expr([$x:tt $($input:tt)*], $weight:tt) -> $ret:tt) => {
		compile_error!(concat!("expected const expression, found `", stringify!($x), "`"))
	};
	
	 // Parse Binary Operations:
	(@ parse_binop($lhs:expr, [$op:tt $($input:tt)*], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ weigh_binop($op, $weight)
			-> parse_binop_rhs(in!, [$($input)*], $weight, $lhs, $op)
			-> $ret }
	};
	(@ parse_binop($lhs:expr, [/*No Tokens*/], $weight:tt) -> $ret:tt) => {
		$crate::poly!{@ return ($lhs, []) -> $ret }
	};
	(@ parse_binop_rhs([], [$($input:tt)*], [], $lhs:expr, $op:tt) -> $ret:tt) => {
		 // Invalid Operator:
		compile_error!(concat!("expected arithmetic operator, found `", stringify!($op), "`"))
	};
	(@ parse_binop_rhs([], [$($input:tt)*], $weight:tt, $lhs:expr, $op:tt) -> $ret:tt) => {
		 // Found Operator of Lower/Same Precedence:
		$crate::poly!{@ return ($lhs, [$op $($input)*]) -> $ret }
	};
	(@ parse_binop_rhs($op_weight:tt, $input:tt, $weight:tt, $lhs:expr, ^) -> $ret:tt) => {
		 // Found Operator `^`, Parse Rhs as Const & Continue Search:
		$crate::poly!{@ parse_const_expr($input, $op_weight)
			-> parse_binop([$lhs ^]:in!, $weight)
			-> $ret }
	};
	(@ parse_binop_rhs($op_weight:tt, $input:tt, $weight:tt, $lhs:expr, $op:tt) -> $ret:tt) => {
		 // Found Operator of Higher Precedence, Parse Rhs & Continue Search:
		$crate::poly!{@ parse_expr($input, $op_weight)
			-> parse_binop([$lhs $op]:in!, $weight)
			-> $ret }
	};
	
	 // Operator Weight Comparison:
	(@ weigh_binop(+, [])           -> $ret:tt) => { $crate::poly!{@ return ([1])     -> $ret } };
	(@ weigh_binop(-, [])           -> $ret:tt) => { $crate::poly!{@ return ([1])     -> $ret } };
	(@ weigh_binop(*, [$(1)?])      -> $ret:tt) => { $crate::poly!{@ return ([1 1])   -> $ret } };
	(@ weigh_binop(/, [$(1)?])      -> $ret:tt) => { $crate::poly!{@ return ([1 1])   -> $ret } };
	(@ weigh_binop(^, [$(1$(1)?)?]) -> $ret:tt) => { $crate::poly!{@ return ([1 1 1]) -> $ret } };
	(@ weigh_binop($op:tt, $max:tt) -> $ret:tt) => { $crate::poly!{@ return ([])      -> $ret } };
	
	 // Token Handling:
	(@ take_first($first:tt $($rest:tt)*) -> $ret:tt) => {
		$crate::poly!{@ return ($first) -> $ret }
	};
	(@ expand($($t:tt)*)) => {
		$($t)*
	};
	
	 // Recursion Fallback:
	(@ $($t:tt)*) => {
		compile_error!(concat!("cannot find rule `", stringify!(@ $($t)*), "` in this macro"))
	};
	
	 // User Entry Point:
	(|$x:ident /*$(: $basis:ty)?*/| $($input:tt)+) => {
		{
			let $x = <::symb_poly::Var>::default();
			$crate::poly!(|_| $($input)+)
		}
	};
	(|_ /*$(: $basis:ty)?*/| $($input:tt)+) => {
		$crate::poly!{@ parse_expr([$($input)+], []) -> take_first -> expand }
	};
}

#[cfg(test)]
mod tests {
	use crate::time::*;
	use crate::*;
	use super::*;

	fn real_roots<K>(poly: Temporal<K>) -> impl Iterator<Item=f64>
	where
		K: Roots<f64>,
		<K as Roots<f64>>::Output: IntoIterator<Item=f64>,
	{
		poly.inner.roots().into_iter()
			.filter(|r| r.is_finite())
	}

	#[test]
	fn add() {
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(a + b, poly!(|x| 8.6 + 8.4*x + 6.300000000000001*x^2 + 4.5*x^3 + 5.7*x^4));
		assert_eq!(b + a, poly!(|x| 8.6 + 8.4*x + 6.300000000000001*x^2 + 4.5*x^3 + 5.7*x^4));
	}

	#[test]
	fn sub() {
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(
			a - b,
			poly!(|x| -5.6 - 3.4000000000000004*x + 0.10000000000000009*x^2 + 4.5*x^3 + 5.7*x^4)
		);
		assert_eq!(
			b - a,
			poly!(|x| 5.6 + 3.4000000000000004*x - 0.10000000000000009*x^2 - 4.5*x^3 - 5.7*x^4)
		);
	}

	#[test]
	fn scalar_mul() {
		let a = poly!(|x| 1.5 + 2.5*x + 3.2*x^2 + 4.5*x^3 + 5.7*x^4);
		let b = poly!(|x| 7.1 + 5.9*x + 3.1*x^2);
		assert_eq!(poly!(|_| 1.5 * a), poly!(|x| 2.25 + 3.75*x + 4.800000000000001*x^2 + 6.75*x^3 + 8.55*x^4));
		assert_eq!(poly!(|_| 1.5 * b), poly!(|x| 10.649999999999999 + 8.850000000000001*x + 4.65*x^2));
	}

	fn assert_roots<K: Poly<f64>>(p: K, expected_roots: &[f64])
	where
		K: Roots<f64>,
		<K as Roots<f64>>::Output: IntoIterator<Item=f64>,
	{
		let mut r = real_roots(Temporal::from(p)).into_iter()
			.collect::<Vec<f64>>();
		let mut expected_roots = expected_roots.into_iter()
			.copied()
			.collect::<Vec<f64>>();
		r.sort_unstable_by(f64::total_cmp);
		expected_roots.sort_unstable_by(f64::total_cmp);
		assert_eq!(
			r.len(), expected_roots.len(),
			"{:?} vs {:?}",
			r, expected_roots
		);
		for (x, y) in r.into_iter().zip(expected_roots.into_iter()) {
			assert!(
				(x - y).abs() < 1e-9 * 10_f64.powf(1. + y.abs().log10().ceil().max(0.)),
				"{:?} vs {:?}",
				x, y
			);
		}
	}

	#[test]
	fn constant() {
		assert_roots(poly!(|_| 2.), &[]);
		assert_roots(poly!(|_| 0.), &[]);
	}

	#[test]
	fn linear() {
		assert_roots(poly!(|x| 20. + -4.*x), &[5.]);
		assert_roots(poly!(|x| 20. + 4.*x), &[-5.]);
		assert_roots(poly!(|x| 4./3.*x), &[0.]);
	}

	#[test]
	fn quadratic() {
		assert_roots(poly!(|x| 20. + 4.*x + 7.*x^2), &[]);
		assert_roots(
			poly!(|x| 20. + 4.*x - 7.*x^2),
			&[-10./7., 2.]
		);
		let r = real_roots(Temporal::from(poly!(|x| -40./3. - 2./3.*x + 17./100.*x^2)))
			.into_iter().collect::<Vec<_>>();
		assert_roots(
			poly!(|x| 40./3. + 2./3.*x - 17./100.*x^2),
			r.as_slice()
		);
		assert_roots(
			poly!(|x| 4./6.*x - 17./100.*x^2),
			&[0., 200./51.]
		);
	}

	#[test]
	fn cubic() {
		assert_roots(
			poly!(|x| 6. - 2077.5*x - 17000./77.*x^2 + 6712./70.*x^3),
			&[-3.64550618348, 0.00288720188, 5.94514363216]
		);
		assert_roots(
			poly!(|x| -3.6350710512584225e-33
				+ 0.10240000000000019*x
				+ -0.5455003017809628*x^2
				+ 1.*x^3),
			&[3.54987407349455e-32]
		);
		assert_eq!(
			real_roots(Temporal::from(poly!(|x| -236263115684.8131
				- 9476965815.566229*x
				- 95034754.784949*x^2
				+ 1.0*x^3
			))).count(),
			3
		);

		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64) -> impl Roots<f64, Output: IntoIterator<Item=f64>> {
			poly!(|x| a
				+ (b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s))*x
				+ (c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + d/(2.*3.*s*s*s))*x^2
				+ (d/(2.*3.*s*s*s))*x^3)
		}
		assert_roots(
			sum_poly(1000000000., 1., 1., 1., -1.),
			&[4591405718.8520711602007]
		);
		assert_roots(
			sum_poly(1000000000.*60., 1., 1., 1., -1.),
			&[275484343229.41354766068234]
		);
		assert_roots(
			sum_poly(1000000000.*60.*60., 1., 1., 1., -1.),
			&[16529060593863.10213769318]
		);
		assert_roots(
			sum_poly(1000000000.*60.*60., 1000., 300., 10., -4.),
			&[
				-55432951711728.79099027553,
				-13201315814560.382043758431976512,
				95634267526286.173034033963943,
			]
		);
	}

	#[test]
	fn quartic() {
		fn sum_poly(s: f64, a: f64, b: f64, c: f64, d: f64, e: f64) -> impl Roots<f64, Output: IntoIterator<Item=f64>> {
			poly!(|x| a
				+ (b/s + c/(2.*s*s) + (2.*d)/(2.*3.*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s))*x
				+ (c/(2.*s*s) + d/(2.*3.*s*s*s) + (2.*d)/(2.*3.*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s) + (2.*3.*e)/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s))*x^2
				+ (d/(2.*3.*s*s*s) + e/(2.*3.*4.*s*s*s*s) + (2.*e)/(2.*3.*4.*s*s*s*s) + (3.*e)/(2.*3.*4.*s*s*s*s))*x^3
				+ (e/(2.*3.*4.*s*s*s*s))*x^4)
		}
		assert_roots(
			sum_poly(1000000000.*60.*60., 7.5, -6.4, -10., 7.8, 2.7),
			&[
				-51671989204288.265625,
				-5830778969957.388671875,
				2846556446843.169921875,
				13056211727396.474609375,
			]
		);
		let t = 1000000000.*60.*60.;
		assert_roots(
			poly!(|x| 7500. - 1000./(t*t)*x^2 + 27./(t*t*t*t)*x^4),
			&[
				-2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
				-2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				2000000000000. * f64::sqrt(60. - 6.*f64::sqrt(19.)),
				2000000000000. * f64::sqrt(60. + 6.*f64::sqrt(19.)),
			]
		);
		assert_roots(
			poly!(|x| 4.906800313619897e-6
				+ 15145.164260286278*x
				+ 23999.999769993723*x^2
				- 236643.191*x^3
				+ 250000.0*x^4
			),
			&[
				-0.19230902079054427,
				-3.23984644667874e-10,
				0.4732863823239847,
				0.6655954027905443
			]
		);
		assert_roots(
			poly!(|x| -35.99999999999999 + 8.046918796812349e-6*x + 2999.99988981303*x^2 - 54772.25541522836*x^3 + 250000.0*x^4),
			&[-0.067702232, 0.177246743]
		);
		assert_eq!(Temporal::new(poly!(|x| -181.99999999999994 - 7.289202428347473e-6*x - 500.*x^2), 1209618449*NANOSEC).eval(1209618450*NANOSEC), -181.99999999999994);
		assert_roots(
			poly!(|x| 9.094947017729282e-13 - 1.7967326421342023e-5*x - 8983.663173028655*x^2 + 997.5710159206409*x^3 + 250000.*x^4),
			&[-0.191570016, -1.11113e-8, 9.11132e-9, 0.187579734]
		);
	}
}