//! ...

use std::cmp::Ordering;
use std::ops::{Add, Mul};

use crate::linear::{Linear, LinearIso, LinearIsoVec, Scalar};
use crate::time;
use crate::time::{Time, InclusiveTimeRanges};
use crate::kind::*;

/// Function that converts a root value to a Time, or ignores it.
pub(crate) trait TimeFilterMap: Clone {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time>;
}

impl<T: Fn(Time, bool) -> Option<Time> + Clone> TimeFilterMap for T {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time> {
		self(time, is_end)
	}
}

/// ...
#[derive(Clone)]
pub struct PreTimeFilterMap;

impl TimeFilterMap for PreTimeFilterMap {
	fn cool(&self, time: Time, is_end: bool) -> Option<Time> {
		if is_end {
			Some(time)
		} else {
			time.checked_sub(time::NANOSEC)
		}
	}
}

/// ...
pub struct DiffTimeFilterMap<A, B, D, I, J, L> {
	a_poly: Poly<A, I>,
	b_poly: Poly<B, J>,
	diff_poly: Poly<D, L>,
}

impl<A, B, D, I, J, L> Clone for DiffTimeFilterMap<A, B, D, I, J, L>
where
	A: Clone,
	B: Clone,
	D: Clone,
{
	fn clone(&self) -> Self {
		Self {
			a_poly: self.a_poly.clone(),
			b_poly: self.b_poly.clone(),
			diff_poly: self.diff_poly.clone(),
		}
	}
}

impl<A, B, D, I, J, L> TimeFilterMap for DiffTimeFilterMap<A, B, D, I, J, L>
where
	A: FluxKind,
	B: FluxKind<Value=A::Value>,
	D: FluxKind<Value=A::Value>,
	A::Value: PartialEq,
	I: LinearIso<A::Value>,
	J: LinearIso<A::Value>,
	L: LinearIso<A::Value>,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		
		let Self { a_poly, b_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).sign();
		
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = if is_end {
				time.checked_add(inc_time)
			} else {
				time.checked_sub(inc_time)
			} {
				 // Stop Before Rate Reverses:
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				if I::linear_id(a_poly.at(next_time)) != J::linear_id(b_poly.at(next_time)) {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub struct DisTimeFilterMap<const SIZE: usize, A, B, D, E, F, I, J, L>
where
	E: FluxKind,
	F: FluxKind,
{
	a_pos: PolyVec<SIZE, A, I>,
	b_pos: PolyVec<SIZE, B, J>,
	dis_poly: Poly<D, L>,
	pos_poly: Poly<E, E::Value>,
	diff_poly: Poly<F, F::Value>,
}

impl<const SIZE: usize, A, B, D, E, F, I, J, L> Clone
	for DisTimeFilterMap<SIZE, A, B, D, E, F, I, J, L>
where
	A: Clone,
	B: Clone,
	D: Clone,
	E: FluxKind,
	F: FluxKind,
{
	fn clone(&self) -> Self {
		Self {
			a_pos: self.a_pos.clone(),
			b_pos: self.b_pos.clone(),
			dis_poly: self.dis_poly.clone(),
			pos_poly: self.pos_poly.clone(),
			diff_poly: self.diff_poly.clone(),
		}
	}
}

impl<const SIZE: usize, A, B, D, E, F, I, J, L> TimeFilterMap
	for DisTimeFilterMap<SIZE, A, B, D, E, F, I, J, L>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	D: FluxKind,
	D::Value: PartialEq + Mul<Output=D::Value>,
	A::Kind: FluxKind<Value=D::Value>,
	B::Kind: FluxKind<Value=D::Value>,
	E: FluxKind<Value=D::Value>,
	F: FluxKind<Value=D::Value>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	L: LinearIso<D::Value>,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		// To handle rounding, the lower bound of equality is undershot.
		// For example, a pair of IVec2 points can round towards each other up
		// to `0.5` along each axis, or `sqrt(n)` in n-dimensional distance. 
		
		let Self { a_pos, b_pos, dis_poly, pos_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).sign();
		
		 // Rounding Buffer:
		if !is_end {
			let round_factor = Scalar(0.5 / (SIZE as f64).sqrt());
			loop {
				let mut inc_time = time::NANOSEC;
				while let Some(next_time) = time.checked_sub(inc_time) {
					 // Stop Before Rate Reverses:
					let rate = diff_poly.rate_at(next_time);
					if sign != rate.sign() && !rate.is_zero() {
						if inc_time == time::NANOSEC {
							return Some(time)
						}
						break
					}
					
					 // Calculate Actual Distances:
					let dis = dis_poly.at(next_time);
					let mut a_dis = D::Value::zero();
					let mut b_dis = D::Value::zero();
					let mut real_diff = D::Value::zero();
					for i in 0..SIZE {
						let a = a_pos.index_poly(i).at(next_time);
						let b = b_pos.index_poly(i).at(next_time);
						a_dis = a_dis + a*a;
						b_dis = b_dis + b*b;
						let x = a - b;
						real_diff = real_diff + x*x;
					}
					a_dis = I::Value::linear_id(a_dis.sqrt());
					b_dis = J::Value::linear_id(b_dis.sqrt());
					real_diff = Mul::<Scalar>::mul(real_diff.sqrt() - dis, round_factor);
					let c_dis = L::linear_id(dis);
					
					 // Undershoot Actual Distances:
					if
						a_dis != I::Value::linear_id(a_dis + real_diff) &&
						b_dis != J::Value::linear_id(b_dis + real_diff) &&
						c_dis != L::linear_id(c_dis + real_diff)
					{
						 // Undershoot Predicted Distances:
						let pred_diff = Mul::<Scalar>::mul(
							pos_poly.at(next_time).sqrt() - dis,
							round_factor
						);
						if
							a_dis != I::Value::linear_id(a_dis + pred_diff) &&
							b_dis != J::Value::linear_id(b_dis + pred_diff) &&
							c_dis != L::linear_id(c_dis + pred_diff)
						{
							break
						}
					}
					
					time = next_time;
					inc_time += inc_time;
				}
				if inc_time == time::NANOSEC {
					break
				}
			}
			
			return Some(time)
		}
		
		 // Fully Cover Zero Range:
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = time.checked_add(inc_time) {
				 // Stop Before Rate Reverses:
				let rate = diff_poly.rate_at(next_time);
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				let mut pos = D::Value::zero();
				for i in 0..SIZE {
					let x = I::Value::linear_id(a_pos.index_poly(i).at(next_time))
						- J::Value::linear_id(b_pos.index_poly(i).at(next_time));
					pos = pos + x*x;
				}
				let dis = L::linear_id(dis_poly.at(next_time));
				if pos != dis*dis {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub trait Prediction {
	type TimeRanges: time::TimeRangeIter;
	fn into_time_ranges(self) -> Self::TimeRanges;
	
	fn into_inclusive_time_ranges(self) -> InclusiveTimeRanges<Self::TimeRanges>
	where
		Self: Sized
	{
		InclusiveTimeRanges {
			times: self.into_time_ranges()
		}
	}
	
	/// Decrements the lower bound of each range by 1 nanosecond.  
	fn pre(self) -> PredFilter<Self, PreTimeFilterMap>
	where
		Self: Sized
	{
		PredFilter {
			pred: self,
			filter: PreTimeFilterMap,
		}
	}
}

macro_rules! impl_prediction {
	(for<$($param:ident),*> $pred:ty) => {
		impl<$($param),*> Prediction for $pred
		where
			Self: IntoIterator<Item = time::TimeRange>
		{
			type TimeRanges = <Self as IntoIterator>::IntoIter;
			fn into_time_ranges(self) -> Self::TimeRanges {
				self.into_iter()
			}
		}
		
		impl<$($param,)* T> std::ops::BitAnd<T> for $pred {
			type Output = PredInter<Self, T>;
			fn bitand(self, b_pred: T) -> Self::Output {
				PredInter {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param,)* T> std::ops::BitOr<T> for $pred {
			type Output = PredUnion<Self, T>;
			fn bitor(self, b_pred: T) -> Self::Output {
				PredUnion {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param,)* T> std::ops::BitXor<T> for $pred {
			type Output = PredSymDiff<Self, T>;
			fn bitxor(self, b_pred: T) -> Self::Output {
				PredSymDiff {
					a_pred: self,
					b_pred,
				}
			}
		}
		
		impl<$($param),*> std::ops::Not for $pred {
			type Output = PredInv<Self>;
			fn not(self) -> Self::Output {
				PredInv {
					pred: self,
				}
			}
		}
	};
	($(for<$($param:ident),*> $pred:ty;)+) => {
		$(impl_prediction!{for<$($param),*> $pred})+
	};
}

impl_prediction!{
	for<K, I> Pred<K, I>;
	for<K, I> PredEq<K, I>;
	for<P, F> PredFilter<P, F>;
	for<A, B> PredInter<A, B>;
	for<A, B> PredUnion<A, B>;
	for<A, B> PredSymDiff<A, B>;
	for<P> PredInv<P>;
	for<> DynPred;
}

impl Prediction for Time {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self..=self))
	}
}

impl Prediction for std::ops::Range<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeInclusive<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeTo<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeToInclusive<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeFrom<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeFull {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl<P: Prediction> Prediction for Option<P> {
	type TimeRanges = time::OptionTimeRanges<P::TimeRanges>;
	fn into_time_ranges(self) -> Self::TimeRanges {
		time::OptionTimeRanges::new(self.map(|x| x.into_time_ranges()))
	}
}

/// ...
pub struct Pred<K, I> {
	pub(crate) poly: Poly<K, I>,
	pub(crate) order: Ordering,
}

impl<K, I> Clone for Pred<K, I>
where
	K: Clone,
{
	fn clone(&self) -> Self {
		Self {
			poly: self.poly.clone(),
			order: self.order,
		}
	}
}

impl<K, I> IntoIterator for Pred<K, I>
where
	K: Roots + PartialOrd,
	K::Value: PartialOrd,
	I: LinearIso<K::Value>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_iter(self) -> Self::IntoIter {
		let basis = self.poly.time();
		let basis_order = self.poly
			.initial_order(Time::ZERO)
			.unwrap_or(Ordering::Equal);
		let times = RootFilterMap {
			times: self.poly.into_inner().roots().into_times(),
			basis,
		};
		time::TimeRangeBuilder::new(times, basis_order, self.order)
	}
}

/// ...
pub struct PredEq<K, I> {
	pub(crate) poly: Poly<K, I>,
}

impl<K, I> Clone for PredEq<K, I>
where
	K: Clone,
{
	fn clone(&self) -> Self {
		Self {
			poly: self.poly.clone(),
		}
	}
}

impl<K, I> IntoIterator for PredEq<K, I>
where
	K: Roots + PartialEq,
	K::Value: PartialEq,
	I: LinearIso<K::Value>,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_iter(self) -> Self::IntoIter {
		let basis = self.poly.time();
		let basis_order = if self.poly.into_inner().is_zero() {
			Ordering::Equal
		} else {
			Ordering::Greater
		};
		let times = RootFilterMap {
			times: self.poly.into_inner().roots().into_times(),
			basis,
		};
		time::TimeRangeBuilder::new(times, basis_order, Ordering::Equal)
	}
}

/// ...
pub struct DynPred {
	inner: time::DynTimeRanges,
}

impl DynPred {
	pub fn new<T>(pred: T) -> Self
	where
		T: Prediction,
		T::TimeRanges: Send + Sync + 'static,
	{
		Self {
			inner: time::DynTimeRanges::new(pred.into_time_ranges())
		}
	}
}

impl IntoIterator for DynPred {
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::DynTimeRanges;
	fn into_iter(self) -> Self::IntoIter {
		self.inner
	}
}

/// ...
#[derive(Clone)]
pub struct PredFilter<P, F> {
	pub(crate) pred: P,
	pub(crate) filter: F,
}

impl<P, F> IntoIterator for PredFilter<P, F>
where
	P: Prediction,
	F: TimeFilterMap,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeFilter<P::TimeRanges, F>;
	fn into_iter(self) -> Self::IntoIter {
		time::TimeFilter::new(self.pred.into_time_ranges(), self.filter)
	}
}

/// ...
#[derive(Clone)]
pub struct PredInter<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> IntoIterator for PredInter<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangesInter<A::TimeRanges, B::TimeRanges>;
	fn into_iter(self) -> Self::IntoIter {
		time::TimeRangesInter::new(
			self.a_pred.into_time_ranges(),
			self.b_pred.into_time_ranges(),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredUnion<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> IntoIterator for PredUnion<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangesUnion<A::TimeRanges, B::TimeRanges>;
	fn into_iter(self) -> Self::IntoIter {
		time::TimeRangesUnion::new(
			self.a_pred.into_time_ranges(),
			self.b_pred.into_time_ranges(),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredSymDiff<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> IntoIterator for PredSymDiff<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangesSymDiff<A::TimeRanges, B::TimeRanges>;
	fn into_iter(self) -> Self::IntoIter {
		time::TimeRangesSymDiff::new(
			self.a_pred.into_time_ranges(),
			self.b_pred.into_time_ranges(),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredInv<P> {
	pred: P,
}

impl<P> IntoIterator for PredInv<P>
where
	P: Prediction,
{
	type Item = <Self::IntoIter as Iterator>::Item;
	type IntoIter = time::TimeRangesInv<P::TimeRanges>;
	fn into_iter(self) -> Self::IntoIter {
		time::TimeRangesInv::new(self.pred.into_time_ranges())
	}
}

/// [`crate::Flux::when`] predictive comparison.
pub trait When<B, J> {
	type Pred: Prediction;
	fn when(self, order: Ordering, poly: Poly<B, J>) -> Self::Pred;
}

impl<A, B, I, J> When<B, J> for Poly<A, I>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind,
	I: LinearIso<A::Value>,
	J: LinearIso<B::Value>,
	<A as ops::Sub<B>>::Output: Roots + PartialOrd,
	A::Value: PartialOrd,
{
	type Pred = PredFilter<
		Pred<<A as ops::Sub<B>>::Output, I>,
		DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output, I, J, I>,
	>;
	fn when(self, order: Ordering, poly: Poly<B, J>) -> Self::Pred {
		let diff_poly = self - poly;
		diff_poly
			.when_sign(order, DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B, J> {
	type Pred: Prediction;
	fn when_eq(self, poly: Poly<B, J>) -> Self::Pred;
}

impl<A, B, I, J> WhenEq<B, J> for Poly<A, I>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind,
	I: LinearIso<A::Value>,
	J: LinearIso<B::Value>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	A::Value: PartialEq,
{
	type Pred = PredFilter<PredEq<<A as ops::Sub<B>>::Output, I>, DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output, I, J, I>>;
	fn when_eq(self, poly: Poly<B, J>) -> Self::Pred {
		let diff_poly = self - poly;
		diff_poly
			.when_zero(DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::FluxVec::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B, D, J, L> {
	type Pred: Prediction;
	fn when_dis(
		self,
		poly: PolyVec<SIZE, B, J>,
		order: Ordering,
		dis: Poly<D, L>,
	) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D, I, J, L> WhenDis<SIZE, B, D, J, L>
	for PolyVec<SIZE, A, I>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	A::Kind: ops::Sub<B::Kind>,
	<A::Kind as ops::Sub<B::Kind>>::Output: ops::Sqr,
	<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialOrd,
	<A::Kind as FluxKind>::Value:
		Mul<Output = <A::Kind as FluxKind>::Value> + PartialOrd,
	D: FluxKind<Value = <A::Kind as FluxKind>::Value> + ops::Sqr,
	L: LinearIso<D::Value>,
{
	type Pred = PredFilter<
		Pred<
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output as FluxKind>::Value,
		>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			I, J, L
		>
	>;
	fn when_dis(
		self,
		poly: PolyVec<SIZE, B, J>,
		order: Ordering,
		dis: Poly<D, L>,
	) -> Self::Pred {
		use ops::*;
		
		let basis = self.time();
		
		let mut sum = <<A::Kind as Sub<B::Kind>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).into_inner()
				.sub(poly.index_poly(i).to_time(basis).into_inner())
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_sign(order, DisTimeFilterMap {
				a_pos: self,
				b_pos: poly,
				dis_poly: dis,
				pos_poly: sum,
				diff_poly,
			})
	}
}

/// [`crate::FluxVec::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D, J, L> {
	type Pred: Prediction;
	fn when_dis_eq(
		self,
		poly: PolyVec<SIZE, B, J>,
		dis: Poly<D, L>,
	) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D, I, J, L> WhenDisEq<SIZE, B, D, J, L> for PolyVec<SIZE, A, I>
where
	A: FluxKindVec<SIZE>,
	B: FluxKindVec<SIZE>,
	I: LinearIsoVec<SIZE, A::Value>,
	J: LinearIsoVec<SIZE, B::Value>,
	A::Kind: ops::Sub<B::Kind>,
	<A::Kind as ops::Sub<B::Kind>>::Output: ops::Sqr,
	<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	<A::Kind as FluxKind>::Value:
		Mul<Output = <A::Kind as FluxKind>::Value> + PartialEq,
	D: FluxKind<Value = <A::Kind as FluxKind>::Value> + ops::Sqr,
	L: LinearIso<D::Value>,
{
	type Pred = PredFilter<
		PredEq<
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output as FluxKind>::Value
		>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			<<A::Kind as ops::Sub<B::Kind>>::Output as ops::Sqr>::Output,
			I, J, L
		>
	>;
	fn when_dis_eq(
		self,
		poly: PolyVec<SIZE, B, J>,
		dis: Poly<D, L>,
	) -> Self::Pred {
		use ops::*;
		
		let basis = self.time();
		
		let mut sum = <<A::Kind as Sub<B::Kind>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index_poly(i).into_inner()
				.sub(poly.index_poly(i).to_time(basis).into_inner())
				.sqr();
		}
		
		let sum = Poly::new(sum, basis);
		let diff_poly = sum - dis.sqr();
		
		diff_poly
			.when_zero(DisTimeFilterMap {
				a_pos: self,
				b_pos: poly,
				dis_poly: dis,
				pos_poly: sum,
				diff_poly,
			})
	}
}

#[test]
fn consistent_sign_pred() {
	use crate::sum::Sum;
	fn toast(time: Time) -> Vec<(Time, Time)> {
		let poly = PolyVec::new([
			Sum::new(-2., [5., -2.]),
			Sum::zero()
		], time);
		poly.when_dis(
			PolyVec::new([Sum::<f64, 2>::zero(); 2], time),
			Ordering::Greater,
			Poly::new(crate::Constant::from(1.), time),
		).into_inclusive_time_ranges()
			.map(|(a, b)| (
				a.saturating_sub(time),
				if b == Time::MAX {
					b
				} else {
					b.saturating_sub(time)
				})
			)
			.collect::<Vec<_>>()
	}
	for t in 0..6 {
		assert_eq!(
			toast(Time::from_secs_f64((t as f64) * 0.5)),
			toast(Time::from_secs_f64(((t+1) as f64) * 0.5))
		);
	}
}