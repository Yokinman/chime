//! Change describing utilities ...

use std::cmp::Ordering;
use std::fmt::Debug;
use std::iter::Sum;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub, Div};

use impl_op::impl_op;

use time::{Time, TimeUnit};

use crate::polynomial::Polynomial;
use crate::degree::{Degree, Deg, IsDeg, IsBelowDeg, HasUpDeg, MaxDeg};

/// A value that changes over time (linearly, exponentially, etc.).
#[derive(Debug, Clone)] // ??? Copy may be worth deriving
pub struct Value<Iso: LinearIso, D: IsDeg> {
	value: LinearFluxValue<Iso::Linear>,
	degree: PhantomData<D>,
}

impl<Iso: LinearIso> Value<Iso, Deg<0>> {
	pub fn new(initial_value: impl Into<Iso>) -> Self {
		Self {
			value: LinearFluxValue::new(initial_value.into().inv_map()),
			degree: PhantomData,
		}
	}
}

impl<Iso: LinearIso, D: IsDeg> Value<Iso, D> {
	fn at(&self, time: Time) -> Iso {
		Iso::map(self.value.at(time /* - self.time_offset */))
	}
}

impl<Iso, D> Value<Iso, D> 
where
	Iso: LinearIso,
	D: IsDeg,
{
	pub fn add_change<T>(mut self, rate: impl Into<Value<Iso, T>>, unit: TimeUnit)
		-> Value::<Iso, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<Iso, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData
		}
	}
	
	pub fn sub_change<T>(mut self, rate: impl Into<Value<Iso, T>>, unit: TimeUnit)
		-> Value::<Iso, <D as AddChange<T>>::Sum>
	where
		T: IsDeg,
		D: AddChange<T>,
	{
		let change = Change::new(Scalar(-1.0), rate.into().value, unit);
		self.value.add_change(change);
		Value::<Iso, <D as AddChange<T>>::Sum> {
			value: self.value,
			degree: PhantomData
		}
	}
}

impl<Iso, D> Value<Iso, D>
where
	Iso: LinearIso,
	Iso::Linear: PartialEq,
	D: IsDeg + IsBelowDeg<Deg<5>>,
{
	fn when_eq<T>(&self, eq_value: &Value<Iso, T>) -> Vec<Time>
	where
		T: IsDeg + IsBelowDeg<Deg<5>>
	{
		self.value.when_eq(&eq_value.value)
	}
}

impl<Iso, D> Value<Iso, D>
where
	Iso: LinearIso,
	Iso::Linear: PartialOrd,
	D: IsDeg + IsBelowDeg<Deg<5>>,
{
	fn when<T>(&self, cmp_order: Ordering, cmp_value: &Value<Iso, T>) -> Vec<(Time, Time)>
	where
		T: IsDeg + IsBelowDeg<Deg<5>>
	{
		self.value.when(cmp_order, &cmp_value.value)
	}
}

// impl<Iso: LinearIso> AsRef<Value<Iso::Linear>> for FluxValue<Iso> {
// 	fn as_ref(&self) -> &Value<Iso::Linear> {
// 		&self.value
// 	}
// }

/// A value that linearly changes over time.
#[derive(Debug, Clone)] // ??? Copy may be worth deriving
struct LinearFluxValue<T: LinearValue> {
	initial_value: T,
	change_list: Vec<Change<T>>,
}

impl<T: LinearValue> LinearFluxValue<T> {
	fn new(initial_value: T) -> Self {
		Self {
			initial_value,
			change_list: vec![],
		}
	}
	
	fn degree(&self) -> Degree {
		self.change_list
			.iter()
			.map(|change| change.rate.degree() + 1)
			.max()
			.unwrap_or(0)
	}
	
	fn calc(&self, t: f64, d: f64) -> T {
		let term_iter = self.change_list.iter().map(|change| {
			let v = change.rate.calc(t, d + 1.0);
			let u = t / ((change.time_unit >> TimeUnit::Nanosecs) as f64);
			v * Scalar((u + d) / (d + 1.0)) * change.scalar
		});
		self.initial_value + term_iter.sum()
	}
	
	fn at(&self, time: Time) -> T {
		self.calc((time >> TimeUnit::Nanosecs) as f64, 0.0)
	}
	
	fn coeff_calc(&self, d: f64, d_min: f64) -> T {
		let coeff_iter = self.change_list.iter().map(|change| {
			let coeff = change.rate.coeff_calc(d + 1.0, d_min);
			let scalar = change.scalar;
			if d < d_min {
				let coeff = coeff * Scalar(1.0 / ((change.time_unit >> TimeUnit::Nanosecs) as f64));
				if d == 0.0 {
					coeff * scalar
				} else {
					(coeff + (change.rate.coeff_calc(d + 1.0, d_min + 1.0) * Scalar(d))) * scalar
				}
			} else {
				coeff * Scalar(d) * scalar
			}
		});
		if d == 0.0 && d_min == 0.0 {
			self.initial_value
		} else {
			let coeff_sum = coeff_iter.sum::<T>() * Scalar(1.0 / (d + 1.0));
			if d >= d_min {
				self.initial_value + coeff_sum
			} else {
				coeff_sum
			}
		}
		// https://www.desmos.com/calculator/tn0udn7swy
	}
	
	fn polynomial(&self) -> Polynomial<T> {
		match self.degree() {
			0 => Polynomial::Constant ([self.coeff_calc(0.0, 0.0)]),
			1 => Polynomial::Linear   ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0)]),
			2 => Polynomial::Quadratic([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0)]),
			3 => Polynomial::Cubic    ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0), self.coeff_calc(0.0, 3.0)]),
			4 => Polynomial::Quartic  ([self.coeff_calc(0.0, 0.0), self.coeff_calc(0.0, 1.0), self.coeff_calc(0.0, 2.0), self.coeff_calc(0.0, 3.0), self.coeff_calc(0.0, 4.0)]),
			_ => unreachable!()
		}
	}
	
	fn roots(&self) -> Result<Vec<Time>, Vec<Time>> {
		let offset = 0.0; // (time_offset >> TimeUnit::Nanosecs) as f64;
		
		let into_time = |roots: Vec<f64>| -> Vec<Time> {
			 roots.into_iter()
				.filter_map(|r| {
					debug_assert!(!r.is_nan());
					let r = r - offset;
					if r < 0.0 || r > (u64::MAX as f64) {
						None
					} else {
						Some((r as u64)*TimeUnit::Nanosecs)
					}
				})
				.collect()
		};
		
		self.polynomial().real_roots()
			.map(into_time)
			.map_err(into_time)
	}
	
	fn update(&mut self, time: Time) {
		self.initial_value = self.at(time);
		for change in &mut self.change_list {
			change.rate.update(time);
		}
		// ??? Could probably reuse most of the initial_value calculation for
		// updating each sub-change's initial_value.
	}
	
	fn add_change(&mut self, change: Change<T>) {
		// ??? Maybe use independent time offsets instead of a shared time
		// offset, so less has to be updated when a change is made to a value.
		// Downside - this may complicate value calculation & root finding.
		self.change_list.push(change);
	}
}

impl<T: LinearValue> From<T> for LinearFluxValue<T> {
	fn from(value: T) -> Self {
		LinearFluxValue::new(value)
	}
}

impl<T> LinearFluxValue<T>
where
	T: LinearValue + PartialEq,
{
	fn when_eq(&self, eq_value: &Self) -> Vec<Time> {
		//! Will panic if used for a value over degree 4.
		
		let polynomial = self.polynomial().clone()
			- eq_value.polynomial().clone();
		
		let mut real_roots = polynomial.real_roots().unwrap_or(vec![]);
		
		 // Constant Equality:
		if real_roots.is_empty() {
			if polynomial.iter().all(|&term| term == term*Scalar(0.0)) {
				real_roots.push(0.0)
			} else {
				return vec![]
			}
		}
		
		real_roots.into_iter()
			.filter_map(|t| {
				if t < 0.0 || t > (u64::MAX as f64) {
					None
				} else {
					Some((t as u64)*TimeUnit::Nanosecs)
				}
			})
			.collect()
	}
}

impl<T> LinearFluxValue<T>
where
	T: LinearValue + PartialOrd,
{
	fn when(&self, cmp_order: Ordering, cmp_value: &Self) -> Vec<(Time, Time)> {
		//! Will panic if used for a value over degree 4.
		
		let polynomial = self.polynomial().clone()
			- cmp_value.polynomial().clone();
		
		 // Find Initial Order:
		let mut degree = polynomial.degree();
		let mut initial_order = loop {
			let coeff = *polynomial.term(degree).unwrap();
			let order = coeff.partial_cmp(&(coeff * Scalar(0.0)));
			if degree == 0 || order != Some(Ordering::Equal) {
				break if degree % 2 == 0 {
					order
				} else {
					order.map(|o| o.reverse())
				}
			}
			degree -= 1;
		};
		
		 // Convert Roots to Ranges:
		let range_list = match polynomial.real_roots() {
			Ok(roots) => {
				let mut list = Vec::with_capacity(
					if cmp_order == Ordering::Equal {
						roots.len()
					} else {
						1 + (roots.len() / 2)
					}
				);
				let mut prev_point = if initial_order == Some(cmp_order) {
					Some(f64::NEG_INFINITY)
				} else {
					None
				};
				for point in roots {
					if let Some(prev) = prev_point {
						if prev != point {
							list.push((prev, point));
						}
						prev_point = None;
					} else if cmp_order == Ordering::Equal {
						list.push((point, point));
					} else {
						prev_point = Some(point);
					}
				}
				if let Some(point) = prev_point {
					list.push((point, f64::INFINITY));
				}
				list
			},
			Err(_) => vec![],
		};
		
		 // Convert Ranges to Time Ranges:
		range_list.into_iter()
			.filter_map(|(a, b)| {
				assert!(a <= b);
				if b < 0.0 || a > (u64::MAX as f64) {
					None
				} else {
					Some((
						(a as u64)*TimeUnit::Nanosecs,
						(b as u64)*TimeUnit::Nanosecs
					))
				}
			})
			.collect()
	}
}

/// A linear change over time.
#[derive(Debug, Clone)]
struct Change<T: LinearValue> {
	scalar: Scalar,
	rate: LinearFluxValue<T>,
	time_unit: TimeUnit,
}

impl<T: LinearValue> Change<T> {
	pub fn new(scalar: Scalar, rate: LinearFluxValue<T>, time_unit: TimeUnit) -> Self {
		Self {
			scalar,
			rate,
			time_unit,
		}
	}
}

/// Produces the degree of the resulting change, à la: `max(A, B+1)`. 
pub trait AddChange<D: IsDeg>: IsDeg {
	type Sum: IsDeg;
}

impl<A, B> AddChange<B> for A
where
	A: IsDeg,
	B: IsDeg + HasUpDeg,
	B::Up: MaxDeg<A>,
{
	type Sum = <B::Up as MaxDeg<A>>::Max;
}

/// A scalar value, used for multiplication with any [`LinearValue`].
#[derive(Debug, Copy, Clone)]
pub struct Scalar(pub f64);

impl_op!{ a * b {
	(Scalar, Scalar) as (Scalar(a), Scalar(b)) => {
		if a == 0. || b == 0. {
			Scalar(0.)
		} else {
			Scalar(a * b)
		}
	},
	(f64, Scalar)['commut] => if b.0 == 0. { 0. } else { a * b.0 },
	(f32, Scalar)['commut] => if b.0 == 0. { 0. } else { a * (b.0 as f32) },
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
pub trait LinearValue
where
	Self:
		Sized + Copy + Clone
		+ Add<Output=Self> + Mul<Scalar, Output=Self>
		// + Sub<Output=Self> + Div<Scalar, Output=Self>
		+ Sum
		+ PartialEq + PartialOrd // !!! Temporary, wouldn't work for something like vectors
		+ Debug // ??? Could be optional
	// Scalar:
	// 	Mul<Self, Output=Self>,
{
	fn roots(polynomial: &Polynomial<Self>) -> Result<Vec<f64>, Vec<f64>> {
		//! Returns all real-valued roots of the given polynomial. If not all
		//! roots are known, `Err` should be returned.
		//! 
		//! [`Polynomial::real_roots`] should generally be used instead.
		
		Err(vec![])
	}
}

/// The binary operators for addition & subtraction.
#[derive(Debug, Copy, Clone)]
pub enum LinearOp {
	Add,
	Sub,
}

impl Operator for LinearOp {
	fn scalar(&self) -> Scalar {
		match self {
			Self::Add => Scalar(1.0),
			Self::Sub => Scalar(-1.0),
		}
	}
}

/// The binary operators for multiplication & division.
#[derive(Debug, Copy, Clone)]
pub enum ExponentialOp {
	Mul,
	Div,
}

impl Operator for ExponentialOp {
	fn scalar(&self) -> Scalar {
		match self {
			Self::Mul => Scalar(1.0),
			Self::Div => Scalar(-1.0),
		}
	}
}

/// A binary operator described by a [`LinearIso`].
pub trait Operator: Debug + Copy + Clone {
	fn scalar(&self) -> Scalar;
}

/// A mapping of a vector space that preserves addition & multiplication.
/// 
/// # Properties
/// 
/// - Generally isomorphic  - `inv_map(map(T)) = T`, `map(inv_map(U)) = U`
/// - Vector addition       - `map(A + B) = map(A) • map(B)`
/// - Scalar multiplication - `map(A * S) = map(A) ^ S`
pub trait LinearIso {
	type Linear: LinearValue;
	type Op: Operator;
	fn map(value: Self::Linear) -> Self;
	fn inv_map(self) -> Self::Linear;
	
	// /// Combine terms using an operator.
	// fn sum<T: LinearValue>(op: Self::Op, augend: Self, addend: Self) -> Self {
	// 	let augend = Self::inv_map(augend);
	// 	let addend = Self::inv_map(addend) * op.scalar();
	// 	Self::map(augend + addend)
	// }
}

/// Used to construct a [`Flux`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(TimeUnit::Secs)` 
pub trait Per: Sized {
	fn per(self, time_unit: TimeUnit) -> Flux<Self> {
		Flux {
			rate: self,
			time_unit
		}
	}
}
impl Per for f64 {}
impl Per for u64 {}
impl Per for i64 {}
impl<Iso: LinearIso, D: IsDeg> Per for Value<Iso, D> {}

/// A change over time used with arithmetic operators to construct [`Value`]s.
pub struct Flux<T> {
	rate: T,
	time_unit: TimeUnit,
}

macro_rules! impl_flux_ops {
	($op_trait:ident::$op_fn:ident, $change_fn:ident, $change_type:ident<$t:ty>) => {
		 // Num + Num
		impl $op_trait<Flux<$t>> for $t {
			type Output = Value<$change_type<$t>, Deg<1>>;
			fn $op_fn(self, rhs: Flux<$t>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.time_unit)
			}
		}
		
		 // Num + Value
		impl<D> $op_trait<Flux<Value<$change_type<$t>, D>>> for $t
		where
			D: IsDeg,
			Deg<0>: AddChange<D>,
		{
			type Output = Value<$change_type<$t>, <Deg<0> as AddChange<D>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$change_type<$t>, D>>) -> Self::Output {
				Value::new(self).$change_fn(rhs.rate, rhs.time_unit)
			}
		}
		
		 // Value + Num
		impl<D> $op_trait<Flux<$t>> for Value<$change_type<$t>, D>
		where
			D: IsDeg + AddChange<Deg<0>>,
		{
			type Output = Value<$change_type<$t>, <D as AddChange<Deg<0>>>::Sum>;
			fn $op_fn(self, rhs: Flux<$t>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.time_unit)
			}
		}
		
		 // Value + Value
		impl<D, T> $op_trait<Flux<Value<$change_type<$t>, T>>> for Value<$change_type<$t>, D>
		where
			D: IsDeg + AddChange<T>,
			T: IsDeg,
		{
			type Output = Value<$change_type<$t>, <D as AddChange<T>>::Sum>;
			fn $op_fn(self, rhs: Flux<Value<$change_type<$t>, T>>) -> Self::Output {
				self.$change_fn(rhs.rate, rhs.time_unit)
			}
		}
	}
}
impl_flux_ops!(Add::add, add_change, Linear<f64>);
impl_flux_ops!(Add::add, add_change, Linear<u64>);
impl_flux_ops!(Add::add, add_change, Linear<i64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<f64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<u64>);
impl_flux_ops!(Sub::sub, sub_change, Linear<i64>);
impl_flux_ops!(Mul::mul, add_change, Exponential<f64>);
impl_flux_ops!(Mul::mul, add_change, Exponential<u64>);
impl_flux_ops!(Div::div, sub_change, Exponential<f64>);
impl_flux_ops!(Div::div, sub_change, Exponential<u64>);

/// A linear map that performs no change - pure addition.
#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub struct Linear<T>(T);

impl<T> From<T> for Linear<T> {
	fn from(value: T) -> Linear<T> {
		Linear(value)
	}
}

impl<T> From<T> for Value<Linear<T>, Deg<0>>
where
	Linear<T>: LinearIso
{
	fn from(value: T) -> Self {
		Self {
			value: LinearFluxValue::new(Linear(value).inv_map()),
			degree: PhantomData,
		}
	}
}

impl LinearIso for Linear<f64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value) }
	fn inv_map(self) -> Self::Linear { self.0 }
}

impl LinearIso for Linear<u64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value.round() as u64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

impl LinearIso for Linear<i64> {
	type Linear = f64;
	type Op = LinearOp;
	fn map(value: Self::Linear) -> Self { Self(value.round() as i64) }
	fn inv_map(self) -> Self::Linear { self.0 as Self::Linear }
}

/// A linear map that translates between addition and multiplication.
/// 
/// `map(inv_map(A) + inv_map(B)) <=> Exponential(A * B)`
#[derive(Debug, Copy, Clone, PartialOrd, PartialEq)]
pub struct Exponential<T>(T);

impl<T> From<T> for Exponential<T> {
	fn from(value: T) -> Exponential<T> {
		Exponential(value)
	}
}

impl<T> From<T> for Value<Exponential<T>, Deg<0>>
where
	Exponential<T>: LinearIso
{
	fn from(value: T) -> Self {
		Self {
			value: LinearFluxValue::new(Exponential(value).inv_map()),
			degree: PhantomData,
		}
	}
}

impl LinearIso for Exponential<f64> {
	type Linear = f64;
	type Op = ExponentialOp;
	fn map(value: Self::Linear) -> Self { Self(value.exp()) }
	fn inv_map(self) -> Self::Linear { self.0.ln() }
}

impl LinearIso for Exponential<u64> {
	type Linear = f64;
	type Op = ExponentialOp;
	fn map(value: Self::Linear) -> Self { Self(value.exp().round() as u64) }
	fn inv_map(self) -> Self::Linear { (self.0 as Self::Linear).ln() }
}

#[cfg(test)]
mod value_tests {
	use super::*;
	use TimeUnit::*;
	
	fn linear() -> (
		Value<Linear<i64>, Deg<1>>,
		Value<Linear<i64>, Deg<2>>,
		Value<Linear<i64>, Deg<3>>,
		Value<Linear<i64>, Deg<4>>,
	) {
		let v0 = 3_i64 - 5.per(Nanosecs);
		let v1 = 20 + v0.clone().per(Nanosecs) + v0.clone().per(Nanosecs);
		let v2 = 52 + v1.clone().per(Nanosecs);
		let v3 = 150 + v2.clone().per(Nanosecs) + v0.clone().per(Nanosecs);
		(v0, v1, v2, v3)
	}
	
	// #[test]
	// fn linear_polynomial() {
	// 	let v = linear();
	// 	assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([3.0, -5.0]));
	// 	assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([20.0, 6.0, -10.0]));
	// 	assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([52.0, 20.0, 6.0, -10.0]));
	// 	assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([150.0, 55.0, 15.0, 6.0, -10.0]));
	// }
	
	#[test]
	fn linear_value() {
		let v = linear();
		let n0 = [3, -2, -7, -12, -17, -22, -27, -32, -37, -42];
		let n1 = [20, 16, 2, -22, -56, -100, -154, -218, -292, -376];
		let n2 = [52, 68, 70, 48, -8, -108, -262, -480, -772, -1148];
		let n3 = [150, 216, 279, 315, 290, 160, -129, -641, -1450, -2640];
		for t in 0..10_u64 {
			assert_eq!(v.0.at(t*Nanosecs), n0[t as usize].into());
			assert_eq!(v.1.at(t*Nanosecs), n1[t as usize].into());
			assert_eq!(v.2.at(t*Nanosecs), n2[t as usize].into());
			assert_eq!(v.3.at(t*Nanosecs), n3[t as usize].into());
		}
	}
	
	#[test]
	fn calc_test() {
		let d = 1.0 + 2.0.per(Mins);
		let c = 1.0 + d.per(Mins) + 7.0.per(Secs);
		let c1 = 3.0 + 9.0.per(Millisecs);
		let b = 10.0 + c.per(Microsecs) + c1.per(Mins);
		let a = 30.0 + b.per(Microsecs);
		
		assert_eq!(a.value.polynomial(), Polynomial::Quartic([
			30.0,
			0.018166666666666664,
			3.668167918055556e-6,
			1.1764138889120372e-15,
			2.314814814814815e-29
		]));
	}
	
	#[test]
	fn linear_roots() {
		let (v0, v1, v2, v3) = linear();
		assert_eq!(v0.when_eq(&Value::new(0)), vec![0*Nanosecs]); // 0.6
		assert_eq!(v1.when_eq(&Value::new(0)), vec![2*Nanosecs]); // -1.902, 2.102
		assert_eq!(v2.when_eq(&Value::new(0)), vec![3*Nanosecs]); // 3.892
		assert_eq!(v3.when_eq(&Value::new(0)), vec![5*Nanosecs]); // -3.687, 5.631
	}
	
	fn exponential() -> (
		Value<Exponential<f64>, Deg<1>>,
		Value<Exponential<f64>, Deg<2>>,
		Value<Exponential<f64>, Deg<3>>,
		Value<Exponential<f64>, Deg<4>>,
	) {
		let v0 = 3.2 / 1.1.per(Nanosecs);
		let v1 = 14.515 * v0.clone().per(Nanosecs) * v0.clone().per(Nanosecs);
		let v2 = 30.4 * v1.clone().per(Nanosecs);
		let v3 = 300.0 * v2.clone().per(Nanosecs) * v0.clone().per(Nanosecs);
		(v0, v1, v2, v3)
	}
	
	// #[test]
	// fn exponential_polynomial() {
	// 	let v = exponential();
	// 	fn ln(n: f64) -> f64 { n.ln() }
	// 	assert_eq!(v[0].polynomial(), FluxPolynomial::Linear([ln(3.2), -ln(1.1)]));
	// 	assert_eq!(v[1].polynomial(), FluxPolynomial::Quadratic([ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// 	assert_eq!(v[2].polynomial(), FluxPolynomial::Cubic([ln(30.4), ln(14.515), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// 	assert_eq!(v[3].polynomial(), FluxPolynomial::Quartic([ln(300.0), ln(30.4) - ln(3.2), ln(14.515) + ln(1.1), 2.0*ln(3.2), -2.0*ln(1.1)]));
	// }
	
	#[test]
	fn exponential_value() {
		let (v0, v1, v2, v3) = exponential();
		let n0 = [3.2, 2.9090909, 2.6446281, 2.4042074, 2.1856431];
		let n1 = [14.515, 122.83769, 859.13387, 4965.97682, 23722.647933];
		let n2 = [30.4, 3734.26565, 3208234.11489, 15932016253.16265, 377949612447409.6];
		for t in 0..5_u64 {
			assert!((v0.at(t*Nanosecs).0 - n0[t as usize]).abs() < 0.00001);
			assert!((v1.at(t*Nanosecs).0 - n1[t as usize]).abs() < 0.00001);
			assert!((v2.at(t*Nanosecs).0 - n2[t as usize]).abs() < 0.00001);
		}
		assert!(v3.at(100*Nanosecs).0.abs() < 0.00001);
	}
	
	#[test]
	fn exponential_roots() {
		let (v0, v1, v2, v3) = exponential();
		assert_eq!(v0.when_eq(&Value::new(1.0)), vec![12*Nanosecs]); // 12.204
		assert_eq!(v1.when_eq(&Value::new(1.0)), vec![24*Nanosecs]); // -1.143, 24.551
		assert_eq!(v2.when_eq(&Value::new(1.0)), vec![36*Nanosecs]); // 36.91 
		assert_eq!(v3.when_eq(&Value::new(1.0)), vec![49*Nanosecs]); // -4.752, 49.329
	}
	
	#[test]
	fn exponential_non_pos() {
		let v = 0.0 * 2.0.per(Nanosecs);
		
		// assert_eq!(v.polynomial(), FluxPolynomial::Linear([f64::NEG_INFINITY, f64::ln(2.0)]));
		assert_eq!(v.at(0*Nanosecs).0, 0.0);
		assert_eq!(v.at(u64::MAX*Nanosecs).0, 0.0);
		
		let v = -1.5 * 2.0.per(Nanosecs);
		
		// match v.polynomial() {
		// 	FluxPolynomial::Linear([a, _]) => assert!(a.is_nan()),
		// 	_ => panic!()
		// }
		assert!(v.at(0*Nanosecs).0.is_nan());
		assert!(v.at(u64::MAX*Nanosecs).0.is_nan());
	}
	
	#[test]
	fn long_term() {
		let h = 10.0 + 1.0.per(Hours);
		let n = 10.0 + 1.0.per(Nanosecs);
		assert_eq!(h.at(1*Hours), n.at(1*Nanosecs));
		assert_eq!(h.at(10*Hours), n.at(10*Nanosecs));
		
		let h1 = 30.0 + h.per(Hours);
		let n1 = 30.0 + n.per(Nanosecs);
		assert_eq!(h1.at(1*Hours), n1.at(1*Nanosecs));
		assert_eq!(h1.at(10*Hours), n1.at(10*Nanosecs));
	}
	
	#[test]
	fn discrete_update() {
		let (mut v0, _, _, mut v3) = linear();
		
		v0.value.update(10*Nanosecs);
		let new_v0 = v0.clone() + 7.per(Nanosecs);
		assert_eq!(v0.at(0*Nanosecs), new_v0.at(0*Nanosecs));
		assert_eq!(new_v0.when_eq(&Value::new(0)), vec![23*Nanosecs]); // 23.5
		
		v3.value.update(6*Nanosecs);
		let new_v3 = v3.clone() + 1000.per(Nanosecs);
		assert_eq!(v3.at(0*Nanosecs), new_v3.at(0*Nanosecs));
		assert_eq!(new_v3.when_eq(&Value::new(0)), vec![0*Nanosecs, 3*Nanosecs]); // 0.22, 3.684
	}
	
	#[test]
	fn when() {
		use Ordering::*;
		let (v0, v1, v2, v3) = linear();
		assert_eq!(v0.when(Equal,   &Value::new(-5)),  vec![(1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v0.when(Less,    &Value::new(-5)),  vec![(1*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v0.when(Greater, &Value::new(-5)),  vec![(0*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v1.when(Equal,   &Value::new(-5)),  vec![(2*Nanosecs, 2*Nanosecs)]);
		assert_eq!(v1.when(Less,    &Value::new(-5)),  vec![(2*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v1.when(Greater, &Value::new(-25)), vec![(0*Nanosecs, 3*Nanosecs)]);
		assert_eq!(v2.when(Equal,   &Value::new(71)),  vec![(1*Nanosecs, 1*Nanosecs), (1*Nanosecs, 1*Nanosecs)]);
		assert_eq!(v2.when(Less,    &Value::new(20)),  vec![(3*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v2.when(Greater, &Value::new(69)),  vec![(1*Nanosecs, 2*Nanosecs)]);
		assert_eq!(v3.when(Equal,   &Value::new(140)), vec![(5*Nanosecs, 5*Nanosecs)]);
		assert_eq!(v3.when(Less,    &Value::new(160)), vec![(0*Nanosecs, 0*Nanosecs), (5*Nanosecs, u64::MAX*Nanosecs)]);
		assert_eq!(v3.when(Greater, &Value::new(280)), vec![(2*Nanosecs, 4*Nanosecs)]);
	}
	
	#[test]
	fn when_eq() {
		let (v0, v1, v2, v3) = linear();
		let v: [Value<Linear<i64>, Deg<4>>; 4] = [
			Value { value: v0.value, degree: PhantomData },
			Value { value: v1.value, degree: PhantomData },
			Value { value: v2.value, degree: PhantomData },
			v3
		];
		for value in v {
			assert_eq!(
				value.when_eq(&Value::new(3)),
				value.when(Ordering::Equal, &Value::new(3)).into_iter()
					.map(|(a, _)| a)
					.collect::<Vec<Time>>()
			)
		}
	}
}