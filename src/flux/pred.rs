//! ...

use std::cmp::Ordering;
use std::ops::Add;

use crate::linear::{Linear, LinearIso, Basis, Scalar, Vector};
use crate::time;
use crate::time::Time;
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
pub struct DiffTimeFilterMap<A, B, D> {
	a_poly: Poly<A>,
	b_poly: Poly<B>,
	diff_poly: Poly<D>,
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
	A: FluxKind,
	B: FluxKind<Basis: Basis<Inner = KindLinear<A>>>,
	D: FluxKind<Basis: Basis<Inner = KindLinear<A>>>,
	KindLinear<A>: PartialEq,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		
		let Self { a_poly, b_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).into_inner().sign();
		
		loop {
			let mut inc_time = time::NANOSEC;
			while let Some(next_time) = if is_end {
				time.checked_add(inc_time)
			} else {
				time.checked_sub(inc_time)
			} {
				 // Stop Before Rate Reverses:
				let rate = diff_poly.rate_at(next_time).into_inner();
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				if A::Basis::inner_id(a_poly.at(next_time).into_inner())
					!= B::Basis::inner_id(b_poly.at(next_time).into_inner())
				{
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
pub struct DisTimeFilterMap<const SIZE: usize, A, B, D, E, F> {
	a_pos: Poly<A>,
	b_pos: Poly<B>,
	dis_poly: Poly<D>,
	pos_poly: Poly<E>,
	diff_poly: Poly<F>,
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
	A: Vector<SIZE, Output: FluxKind> + Clone,
	B: Vector<SIZE, Output: FluxKind> + Clone,
	D: FluxKind,
	KindLinear<D>: PartialEq,
	A::Output: FluxKind<Basis: Basis<Inner = KindLinear<D>>>,
	B::Output: FluxKind<Basis: Basis<Inner = KindLinear<D>>>,
	E: FluxKind<Basis: Basis<Inner = KindLinear<D>>>,
	F: FluxKind<Basis: Basis<Inner = KindLinear<D>>>,
{
	fn cool(&self, mut time: Time, is_end: bool) -> Option<Time> {
		// Covers the range of equality, but stops where the trend reverses.
		// To handle rounding, the lower bound of equality is undershot.
		// For example, a pair of IVec2 points can round towards each other up
		// to `0.5` along each axis, or `sqrt(n)` in n-dimensional distance. 
		
		let Self { a_pos, b_pos, dis_poly, pos_poly, diff_poly } = self;
		let sign = diff_poly.rate_at(time).into_inner().sign();
		
		 // Rounding Buffer:
		if !is_end {
			let round_factor = Scalar::from(0.5 / (SIZE as f64).sqrt());
			loop {
				let mut inc_time = time::NANOSEC;
				while let Some(next_time) = time.checked_sub(inc_time) {
					 // Stop Before Rate Reverses:
					let rate = diff_poly.rate_at(next_time).into_inner();
					if sign != rate.sign() && !rate.is_zero() {
						if inc_time == time::NANOSEC {
							return Some(time)
						}
						break
					}
					
					 // Calculate Actual Distances:
					let dis = dis_poly.at(next_time).into_inner();
					let mut a_dis = <KindLinear<D> as Linear>::zero();
					let mut b_dis = <KindLinear<D> as Linear>::zero();
					let mut real_diff = <KindLinear<D> as Linear>::zero();
					for i in 0..SIZE {
						let a = a_pos.index(i).at(next_time).into_inner();
						let b = b_pos.index(i).at(next_time).into_inner();
						a_dis = a_dis.add(a.clone().sqr());
						b_dis = b_dis.add(b.clone().sqr());
						real_diff = real_diff.add(a.sub(b).sqr());
					}
					a_dis = <A::Output as FluxKind>::Basis::inner_id(a_dis.sqrt());
					b_dis = <B::Output as FluxKind>::Basis::inner_id(b_dis.sqrt());
					real_diff = Linear::mul_scalar(real_diff.sqrt().sub(dis.clone()), round_factor);
					let c_dis = D::Basis::inner_id(dis.clone());
					
					 // Undershoot Actual Distances:
					if
						a_dis != <A::Output as FluxKind>::Basis::inner_id(a_dis.clone().add(real_diff.clone())) &&
						b_dis != <B::Output as FluxKind>::Basis::inner_id(b_dis.clone().add(real_diff.clone())) &&
						c_dis != D::Basis::inner_id(c_dis.clone().add(real_diff))
					{
						 // Undershoot Predicted Distances:
						let pred_diff = Linear::mul_scalar(
							pos_poly.at(next_time).into_inner().sqrt().sub(dis),
							round_factor
						);
						if
							a_dis != <A::Output as FluxKind>::Basis::inner_id(a_dis.clone().add(pred_diff.clone())) &&
							b_dis != <B::Output as FluxKind>::Basis::inner_id(b_dis.clone().add(pred_diff.clone())) &&
							c_dis != D::Basis::inner_id(c_dis.clone().add(pred_diff))
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
				let rate = diff_poly.rate_at(next_time).into_inner();
				if sign != rate.sign() && !rate.is_zero() {
					break
				}
				
				 // Stop Before Inequality:
				let mut pos = <KindLinear<D> as Linear>::zero();
				for i in 0..SIZE {
					let x = <A::Output as FluxKind>::Basis::inner_id(a_pos.index(i).at(next_time).into_inner())
						.sub(<B::Output as FluxKind>::Basis::inner_id(b_pos.index(i).at(next_time).into_inner()));
					pos = pos.add(x.sqr());
				}
				let dis = D::Basis::inner_id(dis_poly.at(next_time).into_inner());
				if pos != dis.sqr() {
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
	pub(crate) poly: Poly<K>,
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
	K: Roots + PartialEq,
	KindLinear<K>: PartialOrd,
{
	type TimeRanges = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		let basis_order = self.poly
			.initial_order(Time::ZERO)
			.unwrap_or(Ordering::Equal);
		let times = RootFilterMap {
			times: self.poly.kind.roots().into_times(),
			basis: self.poly.time,
			prev_time: Time::ZERO,
		};
		time::TimeRangeBuilder::new(times, basis_order, self.order)
	}
}

/// ...
pub struct PredEq<K> {
	pub(crate) poly: Poly<K>,
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
	KindLinear<K>: PartialEq,
{
	type TimeRanges = time::TimeRangeBuilder<RootFilterMap<<<K as Roots>::Output as IntoTimes>::TimeIter>>;
	fn into_ranges(self, _time: Time) -> Self::TimeRanges {
		let basis_order = if self.poly.kind.is_zero() {
			Ordering::Equal
		} else {
			Ordering::Greater
		};
		let times = RootFilterMap {
			times: self.poly.kind.roots().into_times(),
			basis: self.poly.time,
			prev_time: Time::ZERO,
		};
		time::TimeRangeBuilder::new(times, basis_order, Ordering::Equal)
	}
}

/// For defining `dyn Prediction` types with an unspecified `TimeRanges` type.
pub trait DynPrediction {
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

/// [`crate::Flux::when`] predictive comparison.
pub trait When<B> {
	type Pred: Prediction;
	fn when(self, order: Ordering, poly: Poly<B>) -> Self::Pred;
}

impl<A, B> When<B> for Poly<A>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind<Basis: Basis<Inner = KindLinear<A>>>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	KindLinear<A>: PartialOrd,
{
	type Pred = PredFilter<
		Pred<<A as ops::Sub<B>>::Output>,
		DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output>,
	>;
	fn when(self, order: Ordering, poly: Poly<B>) -> Self::Pred {
		let diff_poly = self.clone() - poly.clone();
		diff_poly.clone()
			.when_sign(order, DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::Flux::when_eq`] predictive comparison.
pub trait WhenEq<B> {
	type Pred: Prediction;
	fn when_eq(self, poly: Poly<B>) -> Self::Pred;
	// ??? Add `diff_poly` method, make `When` a sub-trait, elide the returned
	// `DiffTimeFilterMap` type, and provide `when(_eq)` by default.
}

impl<A, B> WhenEq<B> for Poly<A>
where
	A: FluxKind + ops::Sub<B>,
	B: FluxKind<Basis: Basis<Inner = KindLinear<A>>>,
	<A as ops::Sub<B>>::Output: Roots + PartialEq,
	KindLinear<A>: PartialEq,
{
	type Pred = PredFilter<PredEq<<A as ops::Sub<B>>::Output>, DiffTimeFilterMap<A, B, <A as ops::Sub<B>>::Output>>;
	fn when_eq(self, poly: Poly<B>) -> Self::Pred {
		let diff_poly = self.clone() - poly.clone();
		diff_poly.clone()
			.when_zero(DiffTimeFilterMap {
				a_poly: self,
				b_poly: poly,
				diff_poly
			})
	}
}

/// [`crate::FluxVector::when_dis`] predictive distance comparison.
pub trait WhenDis<const SIZE: usize, B, D> {
	type Pred: Prediction;
	fn when_dis(self, poly: Poly<B>, order: Ordering, dis: Poly<D>) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D> WhenDis<SIZE, B, D> for Poly<A>
where
	A: FluxKind + Vector<SIZE, Output: FluxKind> + Clone,
	B: FluxKind + Vector<SIZE, Output: FluxKind<Basis: Basis<Inner = KindLinear<A::Output>>>> + Clone,
	A::Output: ops::Sub<B::Output>,
	<A::Output as ops::Sub<B::Output>>::Output: ops::Sqr,
	<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	KindLinear<A::Output>: PartialOrd,
	D: FluxKind<Basis: Basis<Inner = KindLinear<A::Output>>> + ops::Sqr,
{
	type Pred = PredFilter<
		Pred<<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output,
			<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output,
		>
	>;
	fn when_dis(self, poly: Poly<B>, order: Ordering, dis: Poly<D>) -> Self::Pred {
		use ops::*;
		
		let mut sum = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index(i).kind
				.sub(poly.index(i).to_time(self.time).kind)
				.sqr();
		}
		
		let sum = Poly::new(sum, self.time);
		let diff_poly = sum.clone() - dis.clone().sqr();
		
		diff_poly.clone()
			.when_sign(order, DisTimeFilterMap {
				a_pos: self,
				b_pos: poly,
				dis_poly: dis,
				pos_poly: sum,
				diff_poly,
			})
	}
}

/// [`crate::FluxVector::when_dis_eq`] predictive distance comparison.
pub trait WhenDisEq<const SIZE: usize, B, D> {
	type Pred: Prediction;
	fn when_dis_eq(self, poly: Poly<B>, dis: Poly<D>) -> Self::Pred;
}

impl<const SIZE: usize, A, B, D> WhenDisEq<SIZE, B, D> for Poly<A>
where
	A: FluxKind + Vector<SIZE, Output: FluxKind> + Clone,
	B: FluxKind + Vector<SIZE, Output: FluxKind<Basis: Basis<Inner = KindLinear<A::Output>>>> + Clone,
	A::Output: ops::Sub<B::Output>,
	<A::Output as ops::Sub<B::Output>>::Output: ops::Sqr,
	<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output:
		Add<Output = <<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>
		+ ops::Sub<<D as ops::Sqr>::Output,
			Output = <<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>
		+ Roots
		+ PartialEq,
	KindLinear<A::Output>: PartialEq,
	D: FluxKind<Basis: Basis<Inner = KindLinear<A::Output>>> + ops::Sqr,
{
	type Pred = PredFilter<
		PredEq<<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output>,
		DisTimeFilterMap<
			SIZE,
			A, B, D,
			<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output,
			<<A::Output as ops::Sub<B::Output>>::Output as ops::Sqr>::Output,
		>
	>;
	fn when_dis_eq(self, poly: Poly<B>, dis: Poly<D>) -> Self::Pred {
		use ops::*;
		
		let mut sum = <<A::Output as Sub<B::Output>>::Output as Sqr>::Output::zero();
		for i in 0..SIZE {
			sum = sum + self.index(i).kind
				.sub(poly.index(i).to_time(self.time).kind)
				.sqr();
		}
		
		let sum = Poly::new(sum, self.time);
		let diff_poly = sum.clone() - dis.clone().sqr();
		
		diff_poly.clone()
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
	use crate::time::TimeRanges;
	fn toast(time: Time) -> Vec<(Time, Time)> {
		let poly = Poly::new([
			Sum::new(-2., [5., -2.]),
			Sum::zero()
		], time);
		poly.when_dis(
			Poly::new([Sum::<f64, 2>::zero(); 2], time),
			Ordering::Greater,
			Poly::new(crate::Constant::from(1.), time),
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