//! Working with time.

use std::ops::RangeBounds;
use std::cmp::Ordering;

/// An amount of time.
pub type Time = std::time::Duration;

/// Units of time.
mod units {
	use super::Time;
	
	/// 1 nanosecond (ns).
	pub const NANOSEC: Time = Time::from_nanos(1);
	/// 1 microsecond (μs).
	pub const MICROSEC: Time = Time::from_micros(1);
	/// 1 millisecond (ms).
	pub const MILLISEC: Time = Time::from_millis(1);
	/// 1 second (s).
	pub const SEC: Time = Time::from_secs(1);
	/// 1 minute (m).
	pub const MINUTE: Time = Time::from_secs(60);
	/// 1 hour (h).
	pub const HOUR: Time = Time::from_secs(60 * 60);
}
pub use units::*;

/// Iterators over an ordered sequence of [`TimeRange`] values.
pub trait TimeRanges: Iterator<Item=TimeRange> {
	fn inclusive(self) -> InclusiveTimeRanges<Self>
	where
		Self: Sized
	{
		InclusiveTimeRanges::new(self)
	}
}

macro_rules! impl_time_ranges {
	(for<$($param:ident),*> $iter:ty) => {
		impl<$($param),*> TimeRanges for $iter
		where
			Self: Iterator<Item=TimeRange>,
		{}
	};
	($(for<$($param:ident),*> $iter:ty;)+) => {
		$(impl_time_ranges!{for<$($param),*> $iter})+
	};
}

impl_time_ranges!{
	for<> std::iter::Once<TimeRange>;
	for<T> TimeRangeBuilder<T>;
	for<I, F> TimeFilter<I, F>;
	for<A, B> TimeRangesInter<A, B>;
	for<A, B> TimeRangesUnion<A, B>;
	for<A, B> TimeRangesSymDiff<A, B>;
	for<T> InvTimeRanges<T>;
	for<T> OptionTimeRanges<T>;
	for<> DynTimeRanges;
}

/// Upper or lower bound for a [`TimeRange`].
type TimeBound = std::ops::Bound<Time>;

fn reverse_bound(bound: TimeBound, is_end: bool) -> TimeBound {
	match bound {
		TimeBound::Included(x) => TimeBound::Excluded(x),
		TimeBound::Excluded(x) => TimeBound::Included(x),
		TimeBound::Unbounded => if is_end {
			TimeBound::Excluded(Time::MAX)
		} else {
			TimeBound::Excluded(Time::ZERO)
		},
	}
}

/// Range of time. A <= B.
/// 
/// Inclusivity examples:
/// ```text
///                                     :|0  |1  |2 :
/// TimeRange(Included(0), Included(2)) :[_______]  :
/// TimeRange(Included(0), Excluded(2)) :[_______)  :
/// TimeRange(Excluded(0), Included(1)) :(___]      :
/// TimeRange(Excluded(0), Excluded(1)) :(___)      :
/// TimeRange(Excluded(0), Included(0)) :|          : <- Zero-sized range
/// TimeRange(Unbounded,   Included(1)) :____]      :
/// TimeRange(Unbounded,   Unbounded  ) :___________:
/// ```
#[derive(Copy, Clone, Debug)]
pub struct TimeRange(TimeBound, TimeBound);

impl TimeRange {
	pub fn from_range(range: impl RangeBounds<Time>) -> Self {
		Self::try_from_range(range)
			.expect("must be true: lower bound <= upper bound")
	}
	
	pub fn try_from_range(range: impl RangeBounds<Time>) -> Option<Self> {
		let range = TimeRange(
			range.start_bound().map(|x| *x),
			range.end_bound().map(|x| *x),
		);
		if range.cmp_lower_to_upper(&range).is_gt() {
			return None
		}
		Some(range)
	}
	
	/// Compares the lower bounds of two ranges.
	fn cmp_lower(&self, other: &Self) -> Ordering {
		match (&self.0, &other.0) {
			(TimeBound::Included(a), TimeBound::Included(b)) |
			(TimeBound::Excluded(a), TimeBound::Excluded(b)) => {
				a.cmp(b)
			},
			(TimeBound::Included(a), TimeBound::Excluded(b)) => {
				a.cmp(b).then(Ordering::Less)
			},
			(TimeBound::Excluded(a), TimeBound::Included(b)) => {
				a.cmp(b).then(Ordering::Greater)
			},
			(TimeBound::Unbounded, TimeBound::Unbounded) => {
				Ordering::Equal
			},
			(_, TimeBound::Unbounded) => Ordering::Greater,
			(TimeBound::Unbounded, _) => Ordering::Less,
		}
	}
	
	/// Compares the upper bounds of two ranges.
	fn cmp_upper(&self, other: &Self) -> Ordering {
		match (&self.1, &other.1) {
			(TimeBound::Included(a), TimeBound::Included(b)) |
			(TimeBound::Excluded(a), TimeBound::Excluded(b)) => {
				a.cmp(b)
			},
			(TimeBound::Included(a), TimeBound::Excluded(b)) => {
				a.cmp(b).then(Ordering::Greater)
			},
			(TimeBound::Excluded(a), TimeBound::Included(b)) => {
				a.cmp(b).then(Ordering::Less)
			},
			(TimeBound::Unbounded, TimeBound::Unbounded) => {
				Ordering::Equal
			},
			(_, TimeBound::Unbounded) => Ordering::Less,
			(TimeBound::Unbounded, _) => Ordering::Greater,
		}
	}
	
	/// Compares the lower bound of self to the upper bound of other.
	fn cmp_lower_to_upper(&self, other: &TimeRange) -> Ordering {
		match (&self.0, &other.1) {
			(TimeBound::Included(a), TimeBound::Included(b)) => {
				a.cmp(b).then(Ordering::Less)
			},
			(TimeBound::Excluded(a), TimeBound::Excluded(b)) => {
				a.cmp(b).then(Ordering::Greater)
			},
			(TimeBound::Included(a), TimeBound::Excluded(b)) |
			(TimeBound::Excluded(a), TimeBound::Included(b)) => {
				a.cmp(b)
			},
			_ => Ordering::Less,
		}
	}
	
	/// State of overlap between two ranges. Outside, inside, or adjacent.
	fn overlap(&self, other: &TimeRange) -> Overlap {
		match self.cmp_lower_to_upper(other) {
			Ordering::Less => match other.cmp_lower_to_upper(self) {
				Ordering::Less => Overlap::Inside,
				Ordering::Greater => Overlap::Outside,
				Ordering::Equal => Overlap::Adjacent,
			},
			Ordering::Greater => Overlap::Outside,
			Ordering::Equal => Overlap::Adjacent,
		}
	}
	
	/// Intersection of two ranges, if overlapping or adjacent.
	fn and(self, other: TimeRange) -> Option<TimeRange> {
		match self.overlap(&other) {
			Overlap::Outside => None,
			Overlap::Inside | Overlap::Adjacent => {
				let a = if self.cmp_lower(&other).is_lt() {
					other.0
				} else {
					self.0
				};
				let b = if self.cmp_upper(&other).is_lt() {
					self.1
				} else {
					other.1
				};
				Some(TimeRange(a, b))
			},
		}
	}
	
	/// Union of two ranges, if overlapping or adjacent.
	fn or(self, other: TimeRange) -> Result<TimeRange, (TimeRange, TimeRange)> {
		match self.overlap(&other) {
			Overlap::Outside => Err((self, other)),
			Overlap::Inside | Overlap::Adjacent => {
				let a = if self.cmp_lower(&other).is_lt() {
					self.0
				} else {
					other.0
				};
				let b = if self.cmp_upper(&other).is_lt() {
					other.1
				} else {
					self.1
				};
				Ok(TimeRange(a, b))
			},
		}
	}
	
	/// Symmetric difference of two ranges, or union if adjacent.
	fn xor(self, other: TimeRange) -> Result<(TimeRange, TimeRange), TimeRange> {
		match self.overlap(&other) {
			Overlap::Outside => Ok((self, other)),
			Overlap::Inside => { 
				let (x1, x2) = if self.cmp_lower(&other).is_lt() {
					(self.0, other.0)
				} else {
					(other.0, self.0)
				};
				let (y1, y2) = if self.cmp_upper(&other).is_lt() {
					(self.1, other.1)
				} else {
					(other.1, self.1)
				};
				Ok((
					TimeRange(x1, reverse_bound(x2, false)),
					TimeRange(reverse_bound(y1, true), y2)
				))
			},
			Overlap::Adjacent => Err(self.or(other)
				.expect("should always work")),
		}
	}
}

/// Overlap between two ranges.
enum Overlap {
	Outside,
	Inside,
	Adjacent,
}

/// Converts an iterator of [`Time`]s into an iterator of [`TimeRange`]s.
/// 
/// !!! Seal this later.
pub enum TimeRangeBuilder<I> {
	Empty,
	Points(Time, I),
	Exclusive(Time, I),
	Unbounded(I),
}

impl<I> TimeRangeBuilder<I>
where
	I: Iterator<Item = Time>
{
	pub(crate) fn new(
		iter: impl IntoIterator<IntoIter = I>,
		initial_order: Ordering,
		order: Ordering,
	) -> Self {
		let mut iter = iter.into_iter();
		
		 // Starts Unbounded:
		if order == initial_order {
			return TimeRangeBuilder::Unbounded(iter)
		}
		
		 // Bounded by Initial Time:
		if let Some(time) = iter.next() {
			return if order.is_eq() {
				TimeRangeBuilder::Points(time, iter)
			} else {
				TimeRangeBuilder::Exclusive(time, iter)
			}
		}
		
		TimeRangeBuilder::Empty
	}
}

impl<I> Iterator for TimeRangeBuilder<I>
where
	I: Iterator<Item = Time>
{
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		match std::mem::replace(self, TimeRangeBuilder::Empty) {
			TimeRangeBuilder::Empty => None,
			TimeRangeBuilder::Points(time, mut iter) => {
				let t = TimeBound::Included(time);
				
				 // Increase Time by Next Duration:
				while let Some(add_time) = iter.next() {
					if add_time != Time::ZERO {
						*self = TimeRangeBuilder::Points(time + add_time, iter);
						break
					}
				}
				
				Some(TimeRange(t, t))
			},
			TimeRangeBuilder::Exclusive(mut time, mut iter) => {
				while let Some(add_time) = iter.next() {
					if add_time != Time::ZERO {
						let a = TimeBound::Excluded(time);
						time += add_time;
						let b = TimeBound::Excluded(time);
						if let Some(add_time) = iter.next() {
							*self = TimeRangeBuilder::Exclusive(time + add_time, iter);
						}
						return Some(TimeRange(a, b))
					}
					if let Some(add_time) = iter.next() {
						time += add_time;
					}
				}
				Some(TimeRange(TimeBound::Excluded(time), TimeBound::Unbounded))
			},
			TimeRangeBuilder::Unbounded(mut iter) => {
				let first = if let Some(time) = iter.next() {
					if let Some(add_time) = iter.next() {
						*self = TimeRangeBuilder::Exclusive(time + add_time, iter);
					}
					TimeBound::Excluded(time)
				} else {
					TimeBound::Unbounded
				};
				Some(TimeRange(TimeBound::Unbounded, first))
			},
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		match self {
			TimeRangeBuilder::Empty => (0, Some(0)),
			TimeRangeBuilder::Points(_, iter) => {
				let (_, max) = iter.size_hint();
				(1, max.and_then(|x| x.checked_add(1)))
			},
			TimeRangeBuilder::Exclusive(_, iter) |
			TimeRangeBuilder::Unbounded(iter) => {
				let (_, max) = iter.size_hint();
				(0, max.map(|x| 1 + x/2))
			},
		}
	}
}

/// Filters and/or maps the bounds of a [`TimeRanges`].
/// 
/// !!! Seal this later.
pub struct TimeFilter<I, F> {
	times: I,
	filter: F,
}

impl<I, F> TimeFilter<I, F> {
	pub(crate) fn new(times: I, filter: F) -> Self {
		Self { times, filter }
	}
}

impl<I: TimeRanges, F> Iterator for TimeFilter<I, F>
where
	F: crate::pred::TimeFilterMap
{
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(TimeRange(a, b)) = self.times.next() {
			Some(TimeRange(
				match a {
					TimeBound::Included(t) => self.filter.cool(t, false).map(TimeBound::Included),
					TimeBound::Excluded(t) => self.filter.cool(t, false).map(TimeBound::Excluded),
					TimeBound::Unbounded => None,
				}.unwrap_or(TimeBound::Unbounded),
				match b {
					TimeBound::Included(t) => self.filter.cool(t, true).map(TimeBound::Included),
					TimeBound::Excluded(t) => self.filter.cool(t, true).map(TimeBound::Excluded),
					TimeBound::Unbounded => None,
				}.unwrap_or(TimeBound::Unbounded)
			))
		} else {
			None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.times.size_hint()
	}
}

/// Converts an iterator of [`TimeRange`]s into 2-tuples of inclusive time.
pub struct InclusiveTimeRanges<I> {
	times: I
}

impl<I> InclusiveTimeRanges<I> {
	pub(crate) fn new(iter: impl IntoIterator<IntoIter = I>) -> Self {
		let times = iter.into_iter();
		InclusiveTimeRanges { times }
	}
}

impl<I: TimeRanges> Iterator for InclusiveTimeRanges<I> {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		match self.times.next()? {
			TimeRange(TimeBound::Included(a), TimeBound::Included(b)) => {
				debug_assert!(a <= b);
				Some((a, b))
			},
			TimeRange(TimeBound::Included(a), TimeBound::Excluded(b)) => {
				if a == b {
					self.next()
				} else {
					debug_assert!(a < b);
					Some((a, b - NANOSEC))
				}
			},
			TimeRange(TimeBound::Excluded(a), TimeBound::Included(b)) => {
				if a == b {
					self.next()
				} else {
					debug_assert!(a < b);
					Some((a + NANOSEC, b))
				}
			},
			TimeRange(TimeBound::Excluded(a), TimeBound::Excluded(b)) => {
				debug_assert!(a < b);
				if a == b - NANOSEC {
					self.next()
				} else {
					Some((a + NANOSEC, b - NANOSEC))
				}
			},
			TimeRange(TimeBound::Included(a), TimeBound::Unbounded) => {
				Some((a, Time::MAX))
			},
			TimeRange(TimeBound::Excluded(a), TimeBound::Unbounded) => {
				if a == Time::MAX {
					self.next()
				} else {
					Some((a + NANOSEC, Time::MAX))
				}
			},
			TimeRange(TimeBound::Unbounded, TimeBound::Included(b)) => {
				Some((Time::ZERO, b))
			},
			TimeRange(TimeBound::Unbounded, TimeBound::Excluded(b)) => {
				if b == Time::ZERO {
					self.next()
				} else {
					Some((Time::ZERO, b - NANOSEC))
				}
			},
			TimeRange(TimeBound::Unbounded, TimeBound::Unbounded) => {
				Some((Time::ZERO, Time::MAX))
			},
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, upper) = self.times.size_hint();
		(0, upper)
	}
}

/// Inverted [`TimeRanges`].
pub struct InvTimeRanges<I> {
	iter: I,
	prev: Option<TimeBound>,
}

impl<I> InvTimeRanges<I> {
	pub(crate) fn new(iter: I) -> Self {
		Self {
			iter,
			prev: None,
		}
	}
}

impl<I: TimeRanges> Iterator for InvTimeRanges<I> {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(TimeRange(mut a, mut b)) = self.iter.next() {
			a = reverse_bound(a, false);
			b = reverse_bound(b, true);
			let prev = self.prev.replace(b).unwrap_or(TimeBound::Unbounded);
			return if a == prev {
				self.next()
			} else {
				Some(TimeRange(prev, a))
			}
		}
		if let Some(prev) = self.prev.replace(TimeBound::Unbounded) {
			return if prev == TimeBound::Unbounded {
				None
			} else {
				Some(TimeRange(prev, TimeBound::Unbounded))
			}
		}
		Some(TimeRange(TimeBound::Unbounded, TimeBound::Unbounded))
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (min, max) = self.iter.size_hint();
		(min, max.and_then(|x| x.checked_add(1)))
	}
}

/// Intersection of two [`TimeRanges`].
pub struct TimeRangesInter<A, B> {
	iter: OrdTimeRanges<A, B>,
}

impl<A, B> TimeRangesInter<A, B>
where
	A: TimeRanges,
	B: TimeRanges,
{
	pub(crate) fn new(a: A, b: B) -> Self {
		Self {
			iter: OrdTimeRanges::new(a, b)
		}
	}
}

impl<A: TimeRanges, B: TimeRanges> Iterator for TimeRangesInter<A, B> {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		while let Some((a, Some(b))) = self.iter.next() {
			if let Some(range) = a.and(b) {
				return Some(range)
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.iter.size_hint().1)
	}
}

/// Union of two [`TimeRanges`].
pub struct TimeRangesUnion<A, B> {
	iter: OrdTimeRanges<A, B>,
}

impl<A, B> TimeRangesUnion<A, B>
where
	A: TimeRanges,
	B: TimeRanges,
{
	pub(crate) fn new(a: A, b: B) -> Self {
		Self {
			iter: OrdTimeRanges::new(a, b)
		}
	}
}

impl<A: TimeRanges, B: TimeRanges> Iterator for TimeRangesUnion<A, B> {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		match self.iter.next() {
			Some((a, Some(b))) => {
				match a.or(b) {
					Ok(mut range) => {
						while let Some(r) = self.iter.peek()
							.and_then(|(x, _)| range.or(x).ok())
						{
							range = r;
							self.iter.next();
						}
						Some(range)
					},
					Err((a, _)) => {
						Some(a)
					}
				}
			},
			Some((a, None)) => Some(a),
			None => None,
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.iter.size_hint().1)
	}
}

/// Symmetric difference of two [`TimeRanges`].
pub struct TimeRangesSymDiff<A, B> {
	iter: OrdTimeRanges<A, B>,
	range: Option<TimeRange>,
}

impl<A, B> TimeRangesSymDiff<A, B>
where
	A: TimeRanges,
	B: TimeRanges,
{
	pub(crate) fn new(a: A, b: B) -> Self {
		Self {
			iter: OrdTimeRanges::new(a, b),
			range: None,
		}
	}
}

impl<A: TimeRanges, B: TimeRanges> Iterator for TimeRangesSymDiff<A, B> {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		match self.iter.peek() {
			Some((mut a, Some(mut b))) => {
				if let Some(range) = self.range.take() {
					if a.1 == range.1 {
						a = range;
					} else if b.1 == range.1 {
						b = range;
					} else {
						return Some(range)
					}
				}
				
				self.iter.next();
				
				match a.xor(b) {
					Ok((x, y)) => {
						self.range = Some(y);
						Some(x)
					},
					Err(range) => {
						self.range = Some(range);
						self.next()
					},
				}
			},
			Some((a, None)) => {
				self.iter.next();
				self.range.take().or(Some(a))
			},
			None => self.range.take(),
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.iter.size_hint().1)
	}
}

/// Orders two [`InclusiveTimeRanges`] iterators in parallel.
struct OrdTimeRanges<A, B> {
	a_iter: A,
	a_next: Option<TimeRange>,
	b_iter: B,
	b_next: Option<TimeRange>,
}

impl<A: TimeRanges, B: TimeRanges> OrdTimeRanges<A, B> {
	fn new(mut a_iter: A, mut b_iter: B) -> Self {
		let a_next = a_iter.next();
		let b_next = b_iter.next();
		Self {
			a_iter,
			a_next,
			b_iter,
			b_next,
		}
	}
	
	fn peek(&self) -> Option<(TimeRange, Option<TimeRange>)> {
		match (self.a_next, self.b_next) {
			(Some(a), Some(b)) => {
				if a.cmp_upper(&b).is_le() {
					Some((a, Some(b)))
				} else {
					Some((b, Some(a)))
				}
			},
			(Some(a), _) => Some((a, None)),
			(_, Some(b)) => Some((b, None)),
			_ => None,
		}
	}
}

impl<A: TimeRanges, B: TimeRanges> Iterator for OrdTimeRanges<A, B> {
	type Item = (TimeRange, Option<TimeRange>);
	fn next(&mut self) -> Option<Self::Item> {
		match (self.a_next, self.b_next) {
			(Some(a), Some(b)) => match a.cmp_upper(&b) {
				Ordering::Less => {
					self.a_next = self.a_iter.next();
					Some((a, Some(b)))
				},
				Ordering::Greater => {
					self.b_next = self.b_iter.next();
					Some((b, Some(a)))
				},
				Ordering::Equal => {
					self.a_next = self.a_iter.next();
					self.b_next = self.b_iter.next();
					Some((a, Some(b)))
				},
			},
			(Some(a), None) => {
				self.a_next = self.a_iter.next();
				Some((a, None))
			},
			(None, Some(b)) => {
				self.b_next = self.b_iter.next();
				Some((b, None))
			},
			(None, None) => None,
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (a_min, a_max) = self.a_iter.size_hint();
		let (b_min, b_max) = self.b_iter.size_hint();
		let extra_num = (self.a_next.is_some() as usize)
			+ (self.b_next.is_some() as usize);
		(
			a_min.max(b_min) + extra_num,
			if let (Some(a_max), Some(b_max)) = (a_max, b_max) {
				Some(a_max + b_max + extra_num)
			} else {
				None
			}
		)
	}
}

/// ...
pub struct OptionTimeRanges<T> {
	times: Option<T>,
}

impl<T> OptionTimeRanges<T> {
	pub(crate) fn new(times: Option<T>) -> Self {
		Self { times }
	}
}

impl<T> Iterator for OptionTimeRanges<T>
where
	T: TimeRanges
{
	type Item = T::Item;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(times) = self.times.as_mut() {
			times.next()
		} else {
			None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		if let Some(times) = self.times.as_ref() {
			times.size_hint()
		} else {
			(0, Some(0))
		}
	}
}

/// ...
pub struct DynTimeRanges {
	inner: Box<dyn TimeRanges>,
}

impl DynTimeRanges {
	pub(crate) fn new(inner: impl TimeRanges + 'static) -> Self {
		Self { inner: Box::new(inner) }
	}
}

impl Iterator for DynTimeRanges {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		self.inner.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.inner.size_hint()
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	
	fn inter<A, B>(a: InclusiveTimeRanges<A>, b: InclusiveTimeRanges<B>)
		-> InclusiveTimeRanges<TimeRangesInter<A, B>>
	where
		A: TimeRanges,
		B: TimeRanges,
	{
		InclusiveTimeRanges::new(TimeRangesInter::new(a.times, b.times))
	}
	
	fn union<A, B>(a: InclusiveTimeRanges<A>, b: InclusiveTimeRanges<B>)
		-> InclusiveTimeRanges<TimeRangesUnion<A, B>>
	where
		A: TimeRanges,
		B: TimeRanges,
	{
		InclusiveTimeRanges::new(TimeRangesUnion::new(a.times, b.times))
	}
	
	fn sym_diff<A, B>(a: InclusiveTimeRanges<A>, b: InclusiveTimeRanges<B>)
		-> InclusiveTimeRanges<TimeRangesSymDiff<A, B>>
	where
		A: TimeRanges,
		B: TimeRanges,
	{
		InclusiveTimeRanges::new(TimeRangesSymDiff::new(a.times, b.times))
	}
	
	fn inv<T>(ranges: InclusiveTimeRanges<T>) -> InclusiveTimeRanges<InvTimeRanges<T>>
	where
		T: TimeRanges,
	{
		InclusiveTimeRanges::new(InvTimeRanges::new(ranges.times))
	}
	
	#[test]
	fn range_logic() {
		assert!(TimeRange::try_from_range(5*NANOSEC..4*NANOSEC).is_none());
		assert!(TimeRange::try_from_range(5*NANOSEC..5*NANOSEC).is_some());
		let t = SEC;
		let a = || InclusiveTimeRanges::new(TimeRangeBuilder::Unbounded([2*t, 1*t, 7*t, 10*t, 20*t, 0*t, 0*t, NANOSEC].into_iter()));
		let b = || InclusiveTimeRanges::new(TimeRangeBuilder::Exclusive(2*t, [3*t, 15*t, 20*t, 10*t, 0*t, 0*t].into_iter()));
		let c = || InclusiveTimeRanges::new(TimeRangeBuilder::Points(0*t, [7*t, 13*t, 0*t, 30*t - NANOSEC, NANOSEC, NANOSEC].into_iter()));
		assert_eq!(Vec::from_iter(a()), [
			(Time::ZERO, 2*t - NANOSEC),
			(3*t + NANOSEC, 10*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(b()), [
			(2*t + NANOSEC, 5*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(50*t + NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(c()), [
			(0*t, 0*t),
			(7*t, 7*t),
			(20*t, 20*t),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, 50*t),
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(inter(a(), b())), [
			(3*t + NANOSEC, 5*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(50*t + NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(union(a(), b())), [
			(0*t, 2*t - NANOSEC),
			(2*t + NANOSEC, 10*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX)
		]);
		assert_eq!(Vec::from_iter(sym_diff(a(), b())), [
			(0*t, 2*t - NANOSEC),
			(2*t + NANOSEC, 3*t),
			(5*t, 10*t - NANOSEC),
			(40*t + 2*NANOSEC, 50*t)
		]);
		assert_eq!(Vec::from_iter(inter(a(), a())), Vec::from_iter(a()));
		assert_eq!(Vec::from_iter(inter(b(), b())), Vec::from_iter(b()));
		assert_eq!(Vec::from_iter(inter(c(), c())), Vec::from_iter(c()));
		assert_eq!(Vec::from_iter(union(a(), a())), Vec::from_iter(a()));
		assert_eq!(Vec::from_iter(union(b(), b())), Vec::from_iter(b()));
		assert_eq!(Vec::from_iter(union(c(), c())), Vec::from_iter(c()));
		assert_eq!(Vec::from_iter(sym_diff(a(), a())), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(sym_diff(b(), b())), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(sym_diff(c(), c())), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(inv(inv(inv(a())))), [
			(2*t, 3*t),
			(10*t, 20*t),
			(40*t, 40*t + NANOSEC)
		]);
		assert_eq!(Vec::from_iter(inv(b())), [
			(Time::ZERO, 2*t),
			(5*t, 20*t),
			(40*t, 50*t)
		]);
		assert_eq!(Vec::from_iter(inter(c(), a())), [
			(0*t, 0*t),
			(7*t, 7*t),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, 50*t),
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(inter(b(), c())), [
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(union(c(), a())), [
			(Time::ZERO, 2*t - NANOSEC),
			(3*t + NANOSEC, 10*t - NANOSEC),
			(20*t, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(union(b(), c())), [
			(0*t, 0*t),
			(2*t + NANOSEC, 5*t - NANOSEC),
			(7*t, 7*t),
			(20*t, 40*t - NANOSEC),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(sym_diff(c(), a())), [
			(0*t + NANOSEC, 2*t - NANOSEC),
			(3*t + NANOSEC, 7*t - NANOSEC),
			(7*t + NANOSEC, 10*t - NANOSEC),
			(20*t, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, 50*t - 2*NANOSEC),
			(50*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(sym_diff(b(), c())), [
			(0*t, 0*t),
			(2*t + NANOSEC, 5*t - NANOSEC),
			(7*t, 7*t),
			(20*t, 40*t - NANOSEC),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, 50*t),
			(50*t + 2*NANOSEC, Time::MAX),
		]);
	}
}