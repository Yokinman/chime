//! ...

use std::cmp::Ordering;
use std::ops::{Add, Sub};

use crate::linear::{Linear, Basis, Vector};
use crate::time;
use crate::time::Time;
use crate::temporal::Temporal;
use crate::poly::{IntoTimes, Poly, Roots, RootFilterMap, ops::Sqr, Deriv};

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
pub struct DiffTimeFilterMap<A, B, D> {
	a_poly: Temporal<A>,
	b_poly: Temporal<B>,
	diff_poly: Temporal<D>,
}

impl<A, B, D> Clone for DiffTimeFilterMap<A, B, D>
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

impl<A, B, D> TimeFilterMap for DiffTimeFilterMap<A, B, D>
where
	A: Poly<Basis: PartialEq>,
	B: Poly<Basis = A::Basis>,
	D: Poly<Basis = A::Basis> + Deriv,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		
		let Self { a_poly, b_poly, diff_poly } = self;
		let diff_rate = diff_poly.clone().deriv();
		let sign = diff_rate.eval(time).map_inner(|x| x.sign());
		
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = if is_end {
				time.checked_add(inc_time)
			} else {
				time.checked_sub(inc_time)
			} {
				 // Stop Before Rate Reverses:
				let rate = diff_rate.eval(next_time);
				if rate != Basis::zero() && rate.map_inner(|x| x.sign()) != sign {
					break
				}
				
				 // Stop Before Inequality:
				if a_poly.eval(next_time).inner_id().zip_map_inner(
					b_poly.eval(next_time).inner_id(),
					Linear::sub
				) != Basis::zero() {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				if is_end {
					time.checked_add(inc_time)?;
				} else {
					time.checked_sub(inc_time)?;
				}
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub struct DisTimeFilterMap<const SIZE: usize, A, B, D, E, F> {
	a_pos: Temporal<A>,
	b_pos: Temporal<B>,
	dis_poly: Temporal<D>,
	pos_poly: Temporal<E>,
	diff_poly: Temporal<F>,
}

impl<const SIZE: usize, A, B, D, E, F> Clone
	for DisTimeFilterMap<SIZE, A, B, D, E, F>
where
	A: Clone,
	B: Clone,
	D: Clone,
	E: Clone,
	F: Clone,
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

impl<const SIZE: usize, A, B, D, E, F> TimeFilterMap
	for DisTimeFilterMap<SIZE, A, B, D, E, F>
where
	A: Vector<SIZE, Output: Poly<Basis = D::Basis>> + Clone,
	B: Vector<SIZE, Output: Poly<Basis = D::Basis>> + Clone,
	D: Poly<Basis: PartialEq>,
	E: Poly<Basis = D::Basis>,
	F: Poly<Basis = D::Basis> + Deriv,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		// To handle rounding, the lower bound of equality is undershot.
		// For example, a pair of IVec2 points can round towards each other up
		// to `0.5` along each axis, or `sqrt(n)` in n-dimensional distance. 
		
		let Self { a_pos, b_pos, dis_poly, pos_poly, diff_poly } = self;
		let diff_rate = diff_poly.clone().deriv();
		let sign = diff_rate.eval(time).map_inner(|x| x.sign());
		
		 // Rounding Buffer:
		if !is_end {
			let round_factor = Linear::from_f64(0.5 / (SIZE as f64).sqrt());
			loop {
				let mut inc_time = time::NANOSEC;
				while let Some(next_time) = time.checked_sub(inc_time) {
					 // Stop Before Rate Reverses:
					let rate = diff_rate.eval(next_time);
					if rate != Basis::zero() && rate.map_inner(|x| x.sign()) != sign {
						if inc_time == time::NANOSEC {
							return Some(time)
						}
						break
					}
					
					 // Calculate Actual Distances:
					let dis = dis_poly.eval(next_time);
					let mut a_dis = D::Basis::zero();
					let mut b_dis = D::Basis::zero();
					let mut real_diff = D::Basis::zero();
					for i in 0..SIZE {
						let a = a_pos.index(i).eval(next_time);
						let b = b_pos.index(i).eval(next_time);
						a_dis = a_dis.zip_map_inner(a.clone(), |a, b| a.add(b.sqr()));
						b_dis = b_dis.zip_map_inner(b.clone(), |a, b| a.add(b.sqr()));
						real_diff = Basis::each_map_inner(
							[real_diff, a, b],
							|[real_diff, a, b]| real_diff.add(a.sub(b).sqr())
						);
					}
					a_dis = a_dis.map_inner(Linear::sqrt);
					b_dis = b_dis.map_inner(Linear::sqrt);
					real_diff = real_diff.zip_map_inner(
						dis.clone(),
						|x, dis| x.sqrt().sub(dis).mul(round_factor)
					);
					
					 // Undershoot Actual Distances:
					let diff_not_zero = |mut dis: D::Basis, diff: D::Basis| -> bool {
						dis = dis.inner_id();
						dis.clone().zip_map_inner(
							dis.zip_map_inner(diff, Linear::add).inner_id(),
							Linear::sub
						) != Basis::zero()
					};
					if
						diff_not_zero(a_dis.clone(), real_diff.clone()) &&
						diff_not_zero(b_dis.clone(), real_diff.clone()) &&
						diff_not_zero(dis.clone(), real_diff)
					{
						 // Undershoot Predicted Distances:
						let pred_diff = pos_poly.eval(next_time).zip_map_inner(
							dis.clone(),
							|x, dis| x.sqrt().sub(dis).mul(round_factor)
						);
						if
							diff_not_zero(a_dis, pred_diff.clone()) &&
							diff_not_zero(b_dis, pred_diff.clone()) &&
							diff_not_zero(dis, pred_diff)
						{
							break
						}
					}
					
					time = next_time;
					inc_time += inc_time;
				}
				if inc_time == time::NANOSEC {
					time.checked_sub(inc_time)?;
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
				let rate = diff_rate.eval(next_time);
				if rate != Basis::zero() && rate.map_inner(|x| x.sign()) != sign {
					break
				}
				
				 // Stop Before Inequality:
				let mut pos = D::Basis::zero();
				for i in 0..SIZE {
					pos = Basis::each_map_inner(
						[
							pos,
							a_pos.index(i).eval(next_time).inner_id(),
							b_pos.index(i).eval(next_time).inner_id(),
						],
						|[pos, a, b]| {
							pos.add(a.sub(b).sqr())
						}
					);
				}
				let dis = dis_poly.eval(next_time).inner_id();
				if pos.zip_map_inner(dis, |a, b| a.sub(b.sqr())) != Basis::zero() {
					break
				}
				
				time = next_time;
				inc_time += inc_time;
			}
			if inc_time == time::NANOSEC {
				time.checked_add(inc_time)?;
				break
			}
		}
		
		Some(time)
	}
}

/// ...
pub trait Prediction {
	type TimeRanges: time::TimeRanges;
	fn into_ranges(self, time: Time) -> Self::TimeRanges;
	
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
	for<I> time::TimeRangeBuilder<I>;
	for<K> Pred<K>;
	for<K> PredEq<K>;
	for<P, F> PredFilter<P, F>;
	for<A, B> PredInter<A, B>;
	for<A, B> PredUnion<A, B>;
	for<A, B> PredSymDiff<A, B>;
	for<P> PredInv<P>;
	for<> DynPred;
}

impl Prediction for () {
	type TimeRanges = time::TimeRangeBuilder<std::iter::Empty<Time>>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		time::TimeRangeBuilder::Empty
	}
}

impl<I> Prediction for time::TimeRangeBuilder<I>
where
	Self: IntoIterator<IntoIter: time::TimeRanges>,
{
	type TimeRanges = <Self as IntoIterator>::IntoIter;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		self.into_iter()
	}
}

impl Prediction for Time {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self..=self))
	}
}

impl Prediction for std::ops::Range<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeInclusive<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeTo<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeToInclusive<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeFrom<Time> {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl Prediction for std::ops::RangeFull {
	type TimeRanges = std::iter::Once<time::TimeRange>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		std::iter::once(time::TimeRange::from_range(self))
	}
}

impl<P: Prediction> Prediction for Option<P> {
	type TimeRanges = time::OptionTimeRanges<P::TimeRanges>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::OptionTimeRanges::new(self.map(|x| x.into_ranges(time)))
	}
}

/// ...
pub struct Pred<K> {
	pub(crate) poly: Temporal<K>,
	pub(crate) order: Ordering,
}

impl<K> Clone for Pred<K>
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

impl<K> Prediction for Pred<K>
where
	K: Roots + PartialEq + Deriv,
	K::Basis: PartialOrd,
{
	type TimeRanges = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		let basis_order = self.poly.clone()
			.initial_order(Time::ZERO)
			.unwrap_or(Ordering::Equal);
		let times = RootFilterMap {
			times: self.poly.inner.roots().into_times(),
			basis: self.poly.time,
			prev_time: Time::ZERO,
		};
		time::TimeRangeBuilder::new(times, basis_order, self.order)
	}
}

/// ...
pub struct PredEq<K> {
	pub(crate) poly: Temporal<K>,
}

impl<K> Clone for PredEq<K>
where
	K: Clone,
{
	fn clone(&self) -> Self {
		Self {
			poly: self.poly.clone(),
		}
	}
}

impl<K> Prediction for PredEq<K>
where
	K: Roots + PartialEq,
	K::Basis: PartialEq,
{
	type TimeRanges = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		let basis_order = if self.poly.inner.is_zero() {
			Ordering::Equal
		} else {
			Ordering::Greater
		};
		let times = RootFilterMap {
			times: self.poly.inner.roots().into_times(),
			basis: self.poly.time,
			prev_time: Time::ZERO,
		};
		time::TimeRangeBuilder::new(times, basis_order, Ordering::Equal)
	}
}

/// For defining `dyn Prediction` types with an unspecified `TimeRanges` type.
trait DynPrediction {
	#[allow(clippy::wrong_self_convention)]
	fn into_ranges(&mut self, time: Time) -> time::DynTimeRanges;
}

impl<T> DynPrediction for Option<T>
where
	T: Prediction,
	T::TimeRanges: 'static,
{
	fn into_ranges(&mut self, time: Time) -> time::DynTimeRanges {
		time::DynTimeRanges::new(self.take()
			.expect("should exist")
			.into_ranges(time))
	}
}

/// ...
pub struct DynPred {
	inner: Box<dyn DynPrediction>,
}

impl DynPred {
	pub fn new(pred: impl Prediction + 'static) -> Self {
		Self {
			inner: Box::new(Some(pred))
		}
	}
}

impl Prediction for DynPred {
	type TimeRanges = time::DynTimeRanges;
	fn into_ranges(mut self, time: Time) -> Self::TimeRanges {
		self.inner.into_ranges(time)
	}
}

/// ...
#[derive(Clone)]
pub struct PredFilter<P, F> {
	pub(crate) pred: P,
	pub(crate) filter: F,
}

impl<P, F> Prediction for PredFilter<P, F>
where
	P: Prediction,
	F: TimeFilterMap,
{
	type TimeRanges = time::TimeFilter<P::TimeRanges, F>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::TimeFilter::new(self.pred.into_ranges(time), self.filter)
	}
}

/// ...
#[derive(Clone)]
pub struct PredInter<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> Prediction for PredInter<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type TimeRanges = time::TimeRangesInter<A::TimeRanges, B::TimeRanges>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::TimeRangesInter::new(
			self.a_pred.into_ranges(time),
			self.b_pred.into_ranges(time),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredUnion<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> Prediction for PredUnion<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type TimeRanges = time::TimeRangesUnion<A::TimeRanges, B::TimeRanges>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::TimeRangesUnion::new(
			self.a_pred.into_ranges(time),
			self.b_pred.into_ranges(time),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredSymDiff<A, B> {
	a_pred: A,
	b_pred: B,
}

impl<A, B> Prediction for PredSymDiff<A, B>
where
	A: Prediction,
	B: Prediction,
{
	type TimeRanges = time::TimeRangesSymDiff<A::TimeRanges, B::TimeRanges>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::TimeRangesSymDiff::new(
			self.a_pred.into_ranges(time),
			self.b_pred.into_ranges(time),
		)
	}
}

/// ...
#[derive(Clone)]
pub struct PredInv<P> {
	pred: P,
}

impl<P> Prediction for PredInv<P>
where
	P: Prediction,
{
	type TimeRanges = time::InvTimeRanges<P::TimeRanges>;
	fn into_ranges(self, time: Time) -> Self::TimeRanges {
		time::InvTimeRanges::new(self.pred.into_ranges(time))
	}
}

/// [`crate::Temporal::when`] predictive comparison.
pub trait When<B: Poly>: Poly {
	type Pred: Prediction;
	fn when(
		a_poly: Temporal<Self>,
		cmp: Ordering,
		b_poly: Temporal<B>,
	) -> Self::Pred;
}

impl<A, B> When<B> for A
where
	A: Poly + Sub<B, Output: Roots + PartialEq + Poly<Basis = A::Basis> + Deriv>,
	B: Poly<Basis = A::Basis>,
	A::Basis: PartialOrd,
{
	type Pred = PredFilter<
		Pred<<A as Sub<B>>::Output>,
		DiffTimeFilterMap<A, B, <A as Sub<B>>::Output>,
	>;
	fn when(
		a_poly: Temporal<A>,
		cmp: Ordering,
		b_poly: Temporal<B>,
	) -> Self::Pred {
		let diff_poly = Temporal {
			inner: a_poly.inner.clone() - b_poly.inner.clone(),
			time: a_poly.time,
		};
		diff_poly.clone()
			.when_sign(cmp, DiffTimeFilterMap { a_poly, b_poly, diff_poly })
	}
}

/// [`crate::Temporal::when_eq`] predictive comparison.
pub trait WhenEq<B: Poly>: Poly {
	type Pred: Prediction;
	fn when_eq(
		a_poly: Temporal<Self>,
		b_poly: Temporal<B>,
	) -> Self::Pred;
	// ??? Add `diff_poly` method, make `When` a sub-trait, elide the returned
	// `DiffTimeFilterMap` type, and provide `when(_eq)` by default.
}

impl<A, B> WhenEq<B> for A
where
	A: Poly + Sub<B, Output: Roots + PartialEq + Poly<Basis = A::Basis> + Deriv>,
	B: Poly<Basis = A::Basis>,
	A::Basis: PartialEq,
{
	type Pred = PredFilter<
		PredEq<<A as Sub<B>>::Output>,
		DiffTimeFilterMap<A, B, <A as Sub<B>>::Output>,
	>;
	fn when_eq(
		a_poly: Temporal<A>,
		b_poly: Temporal<B>,
	) -> Self::Pred {
		let diff_poly = Temporal {
			inner: a_poly.inner.clone() - b_poly.inner.clone(),
			time: a_poly.time,
		};
		diff_poly.clone()
			.when_zero(DiffTimeFilterMap { a_poly, b_poly, diff_poly })
	}
}

/// [`crate::FluxVector::when_dis`] predictive distance comparison.
pub trait WhenDis<B: Poly, D: Poly, const SIZE: usize>: Poly {
	type Pred: Prediction;
	fn when_dis(
		a_poly: Temporal<Self>,
		b_poly: Temporal<B>,
		cmp: Ordering,
		dis: Temporal<D>,
	) -> Self::Pred;
}

impl<A, B, D, const SIZE: usize> WhenDis<B, D, SIZE> for A
where
	A: Poly + Vector<SIZE, Output: Poly<Basis = D::Basis>>,
	B: Poly + Vector<SIZE, Output: Poly<Basis = D::Basis>>,
	D: Poly<Basis: PartialOrd> + Sqr,
	A::Output: Sub<
		B::Output,
		Output: Sqr<Output:
			Add<Output = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output>
			+ Sub<<D as Sqr>::Output,
				Output = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output>
			+ Roots
			+ PartialEq
			+ Poly<Basis = D::Basis>
			+ Deriv
		>,
	>,
{
	type Pred = PredFilter<
		Pred<<<A::Output as Sub<B::Output>>::Output as Sqr>::Output>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Output as Sub<B::Output>>::Output as Sqr>::Output,
			<<A::Output as Sub<B::Output>>::Output as Sqr>::Output,
		>
	>;
	fn when_dis(
		a_pos: Temporal<A>,
		b_pos: Temporal<B>,
		cmp: Ordering,
		dis_poly: Temporal<D>,
	) -> Self::Pred {
		let mut sum = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + (a_pos.index(i).inner - b_pos.index(i).inner).sqr();
		}
		
		let pos_poly = Temporal::new(sum, a_pos.time);
		let diff_poly = Temporal {
			inner: pos_poly.inner.clone() - dis_poly.clone().sqr().inner,
			time: a_pos.time,
		};
		
		diff_poly.clone().when_sign(
			cmp,
			DisTimeFilterMap { a_pos, b_pos, dis_poly, pos_poly, diff_poly },
		)
	}
}

/// [`crate::FluxVector::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<B: Poly, D: Poly, const SIZE: usize>: Poly {
	type Pred: Prediction;
	fn when_dis_eq(
		a_pos: Temporal<Self>,
		b_pos: Temporal<B>,
		dis: Temporal<D>,
	) -> Self::Pred;
}

impl<A, B, D, const SIZE: usize> WhenDisEq<B, D, SIZE> for A
where
	A: Poly + Vector<SIZE, Output: Poly<Basis = D::Basis>>,
	B: Poly + Vector<SIZE, Output: Poly<Basis = D::Basis>>,
	D: Poly<Basis: PartialEq> + Sqr,
	A::Output: Sub<B::Output,
		Output: Sqr<Output:
			Add<Output = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output>
			+ Sub<<D as Sqr>::Output,
				Output = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output>
			+ Roots
			+ PartialEq
			+ Poly<Basis = D::Basis>
			+ Deriv
		>,
	>,
{
	type Pred = PredFilter<
		PredEq<<<A::Output as Sub<B::Output>>::Output as Sqr>::Output>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Output as Sub<B::Output>>::Output as Sqr>::Output,
			<<A::Output as Sub<B::Output>>::Output as Sqr>::Output,
		>
	>;
	fn when_dis_eq(
		a_pos: Temporal<A>,
		b_pos: Temporal<B>,
		dis_poly: Temporal<D>,
	) -> Self::Pred {
		let mut sum = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + (a_pos.index(i).inner - b_pos.index(i).inner).sqr();
		}
		
		let pos_poly = Temporal::new(sum, a_pos.time);
		let diff_poly = Temporal {
			inner: pos_poly.inner.clone() - dis_poly.clone().sqr().inner,
			time: a_pos.time,
		};
		
		diff_poly.clone()
			.when_zero(DisTimeFilterMap { a_pos, b_pos, dis_poly, pos_poly, diff_poly })
	}
}

#[test]
fn consistent_sign_pred() {
	use crate::time::TimeRanges;
	fn toast(time: Time) -> Vec<(Time, Time)> {
		let poly = Temporal::new([
			symb_poly::Invar(crate::constant::Constant2(-2.))
				+ symb_poly::Invar(crate::constant::Constant2(5.)) * <symb_poly::Var>::default()
				+ symb_poly::Invar(crate::constant::Constant2(-2.)) * <symb_poly::Var>::default() * <symb_poly::Var>::default(),
			symb_poly::Invar(crate::constant::Constant2(0.))
				+ symb_poly::Invar(crate::constant::Constant2(0.)) * <symb_poly::Var>::default()
				+ symb_poly::Invar(crate::constant::Constant2(0.)) * <symb_poly::Var>::default() * <symb_poly::Var>::default(),
		], time);
		poly.when_dis(
			Temporal::new([symb_poly::Invar(crate::constant::Constant2(0.)); 2], time),
			Ordering::Greater,
			Temporal::new(1., time).to_poly(),
		).into_ranges(Time::ZERO)
			.inclusive()
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