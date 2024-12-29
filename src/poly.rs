use std::cmp::Ordering;
use crate::linear::Basis;
use crate::time::Time;

pub use crate::change::constant::Constant;
pub use crate::change::sum::SumPoly;
pub use crate::change::prod::ProdPoly;
use crate::change::sum::Monomial;
pub use crate::change::sumprod::SumProdPoly;

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
		
		Deriv::deriv(self).initial_order(time).map(Ordering::reverse)
	}
	
	fn zero() -> Self;
	
	fn is_zero(&self) -> bool
	where
		Self: PartialEq
	{
		self.eq(&Self::zero())
	}
}

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
	type Deriv: Deriv<Basis = Self::Basis>;
	fn deriv(self) -> Self::Deriv;
}

/// ...
pub trait Translate: Poly {
	type Output: Poly<Basis = Self::Basis>;
	fn translate(self, amount: <Self::Basis as Basis>::Inner) -> Self::Output;
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

/// ...
#[derive(Copy, Clone, Debug)]
pub struct Variable<const D: usize>;

impl<const T: usize> IntoTerm for Variable<T> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

impl<const N: usize> std::ops::BitXor<Exponent<N>> for Variable<1> {
	type Output = Variable<N>;
	fn bitxor(self, _rhs: Exponent<N>) -> Self::Output {
		Variable::<N>
	}
}

impl<T, const D: usize> std::ops::Mul<Variable<D>> for Constant<T> {
	type Output = Monomial<T, D>;
	fn mul(self, _rhs: Variable<D>) -> Self::Output {
		Monomial(self.0)
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
	type Term = Constant<T>;
	fn into_term(self) -> Self::Term {
		Constant(self)
	}
}

impl<T> IntoTerm for Constant<T> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

impl<T, const D: usize> IntoTerm for Monomial<T, D> {
	type Term = Self;
	fn into_term(self) -> Self::Term {
		self
	}
}

impl<A, B> IntoTerm for crate::change::sum::Binomial<A, B> {
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
			let $x = $crate::poly::Variable::<1>;
			$crate::poly!(|_| $($input)+)
		}
	};
	(|_ /*$(: $basis:ty)?*/| $($input:tt)+) => {
		$crate::poly!{@ parse_expr([$($input)+], []) -> take_first -> expand }
	};
}