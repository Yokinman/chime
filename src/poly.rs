use std::cmp::Ordering;
use crate::linear::Basis;
use crate::time::Time;

pub use crate::change::constant::Constant2;

/// An abstract description of change over time.
/// 
/// Used to define the standard representations of [`Flux`] types. In
/// other words, the layout of a polynomial.
pub trait Poly: Clone {
	const DEGREE: usize;
	
	type Basis: Basis;
	
	fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis;
	
	/// The order at or immediately preceding the value at a time.
	/// 
	/// This should be the first non-zero [`Poly::eval`] value of this kind
	/// or its derivatives; reversed for odd derivatives.
	fn initial_order(self, time: <Self::Basis as Basis>::Inner) -> Option<Ordering>
	where
		Self: Deriv,
		Self::Basis: PartialOrd,
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
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

mod _impl_poly {
	use crate::constant::Constant2;
	use crate::linear::{Basis, Linear};
	use super::Poly;
	use symb_poly::{Invar, Func, Sum, Power, Prod};
	
	impl<T: Poly, const SIZE: usize> Poly for [T; SIZE] {
		const DEGREE: usize = T::DEGREE;
		type Basis = [T::Basis; SIZE];
		fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			self.each_ref().map(|x| T::eval(x, time))
		}
		fn zero() -> Self {
			std::array::from_fn(|_| T::zero())
		}
	}
	
	impl Poly for Invar<symb_poly::Num<typenum::Z0>> {
		const DEGREE: usize = 0;
		
		type Basis = f64;
		
		fn eval(&self, _time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			0.
		}
		
		fn zero() -> Self {
			Self::default()
		}
	}
	
	impl<T> Poly for Invar<Constant2<T>>
	where
		T: Basis
	{
		const DEGREE: usize = 0;
		
		type Basis = T;
		
		fn eval(&self, _time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			let Invar(Constant2(basis)) = self.clone();
			basis
		}
		
		fn zero() -> Self {
			Invar(Constant2(T::zero()))
		}
	}
	
	impl<A, B> Poly for Func<Sum, (A, B)>
	where
		A: Poly,
		B: Poly<Basis = A::Basis>,
	{
		const DEGREE: usize = 0;
		
		type Basis = A::Basis;
		
		fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			let (a, b) = self.args();
			a.eval(time).zip_map_inner(b.eval(time), Linear::add)
		}
		
		fn zero() -> Self {
			Func::from_args((A::zero(), B::zero()))
		}
	}
	
	impl<A, B> Poly for Func<Prod, (A, B)>
	where
		A: Poly,
		B: Clone + Default + symb_poly::Replace<symb_poly::Variable, Invar<Constant2<A::Basis>>,
			Output = Invar<Constant2<A::Basis>>>,
	{
		const DEGREE: usize = 0;
		
		type Basis = A::Basis;
		
		fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			let (a, b) = self.args();
			let Invar(Constant2(b_eval)) = b.clone().replace(Invar(Constant2(A::Basis::from_inner(time))));
			a.eval(time).zip_map_inner(b_eval, Linear::mul)
		}
		
		fn zero() -> Self {
			Func::from_args((A::zero(), B::default()))
		}
	}
	
	impl<A, B> Poly for Func<Power, (A, B)>
	where
		A: Poly,
		B: Clone + Default + symb_poly::Replace<symb_poly::Variable, Invar<Constant2<A::Basis>>,
			Output = Invar<Constant2<A::Basis>>>,
	{
		const DEGREE: usize = 0;
		
		type Basis = A::Basis;
		
		fn eval(&self, time: <Self::Basis as Basis>::Inner) -> Self::Basis {
			let (a, b) = self.args();
			let Invar(Constant2(b_eval)) = b.clone().replace(Invar(Constant2(A::Basis::from_inner(time))));
			a.eval(time).zip_map_inner(b_eval, Linear::pow)
		}
		
		fn zero() -> Self {
			Func::from_args((A::zero(), B::default()))
		}
	}
}

impl<T, const N: usize> Translate for [T; N]
where
	T: Translate
{
	type Output = [T::Output; N];
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
		self.map(|x| x.translate(amount))
	}
}

impl<T, const N: usize> Deriv for [T; N]
where
	T: Deriv
{
	type Deriv = [T::Deriv; N];
	fn deriv(self) -> Self::Deriv {
		self.map(<T as Deriv>::deriv)
	}
}

/// ...
pub trait Deriv: Poly {
	type Deriv: Poly<Basis = Self::Basis>;
	fn deriv(self) -> Self::Deriv;
}

mod _impl_deriv {
	use super::{Deriv, Poly};

	impl<T> Deriv for symb_poly::Invar<T>
	where
		Self: Poly + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<Basis = <Self as Poly>::Basis>>
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
	
	impl<T> Deriv for symb_poly::Var<T>
	where
		Self: Poly + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<Basis = <Self as Poly>::Basis>>
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
	
	impl<F, A> Deriv for symb_poly::Func<F, A>
	where
		Self: Poly + symb_poly::DerivExpr<symb_poly::Variable, Output: Poly<Basis = <Self as Poly>::Basis>>
	{
		type Deriv = <Self as symb_poly::DerivExpr<symb_poly::Variable>>::Output;
		fn deriv(self) -> Self::Deriv {
			symb_poly::DerivExpr::deriv_expr(self, Default::default())
		}
	}
}

/// ...
pub trait Translate: Poly {
	type Output: Poly<Basis = Self::Basis>;
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output;
}

mod _impl_translate {
	use super::Translate;
	use symb_poly::{Func, Sum, Invar, Var, Replace, Variable};
	use crate::linear::Basis;
	use crate::constant::Constant2;
	use crate::poly::Poly;
	
	impl<A> Translate for Invar<A>
	where
		Self: Poly
	{
		type Output = Self;
		fn translate(self, _amount: <Self::Basis as Basis>::Inner) -> Self::Output {
			self
		}
	}
	
	impl<F, A> Translate for Func<F, A>
	where
		Self: Poly + Replace<Variable, Func<Sum, (Invar<Constant2<<Self as Poly>::Basis>>, Var)>, Output: Poly<Basis = <Self as Poly>::Basis>>
	{
		type Output = <Self as Replace<Variable, Func<Sum, (Invar<Constant2<<Self as Poly>::Basis>>, Var)>>>::Output;
		fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output {
			self.replace(Func::from_args((Invar(Constant2(Basis::from_inner(amount))), Default::default())))
		}
	}
}

/// Combining [`Poly`] types.
/// 
/// Primarily this serves as a way to put two kinds of change-over-time into
/// the same space, for combination or comparison purposes.
pub mod ops {
	use std::ops;
	use crate::linear::Basis;
	use super::Poly;
	
	/// Squaring a kind of change.
	pub trait Sqr: Poly {
		type Output: Poly<Basis: Basis<Inner = <Self::Basis as Basis>::Inner>>;
		fn sqr(self) -> <Self as Sqr>::Output;
	}
	
	impl<K: Poly> Sqr for K
	where
		K: ops::Mul,
		<K as ops::Mul>::Output: Poly<Basis: Basis<Inner = <K::Basis as Basis>::Inner>>,
	{
		type Output = <K as ops::Mul>::Output;
		fn sqr(self) -> <Self as Sqr>::Output {
			self.clone() * self
		}
	}
}

/// Roots of a [`Poly`] polynomial.
/// 
/// For discontinuous change-over-time, roots should also include any moments
/// where the polynomial discontinuously "teleports" across 0.
pub trait Roots: Poly {
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
	type Term = symb_poly::Invar<Constant2<T>>;
	fn into_term(self) -> Self::Term {
		symb_poly::Invar(Constant2(self))
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