//! Working with time.

use std::ops::{BitAnd, BitOr, BitXor, Not, RangeBounds};
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

/// Iterator types usable by [`TimeRanges`].
pub trait TimeRangeIter: Iterator<Item=TimeRange> + Send + Sync + Clone + 'static {}
impl<T: Iterator<Item=TimeRange> + Send + Sync + Clone + 'static> TimeRangeIter for T {}

/// Iterator types usable by [`TimeRangeBuilder`].
pub(crate) trait TimeIter: Iterator<Item=Time> + Send + Sync + Clone + 'static {}
impl<T: Iterator<Item=Time> + Send + Sync + Clone + 'static> TimeIter for T {}

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
///                                     :|0| |1| |2|:
/// TimeRange(Included(0), Included(2)) :|_________|:
/// TimeRange(Included(0), Excluded(2)) :|_______|  :
/// TimeRange(Excluded(0), Included(1)) :  |___|    :
/// TimeRange(Excluded(0), Excluded(1)) :  |_|      :
/// TimeRange(Excluded(0), Included(0)) :  |        :
/// TimeRange(Unbounded,   Included(1)) :______|    :
/// TimeRange(Unbounded,   Unbounded  ) :___________:
/// ```
#[derive(Copy, Clone, Debug)]
pub struct TimeRange(TimeBound, TimeBound);

impl TimeRange {
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
#[must_use]
#[derive(Clone)]
pub(crate) enum TimeRangeBuilder<I> {
	Inclusive(I, Option<Time>),
	Exclusive(I),
	Unbounded(Option<I>),
}

impl<I: TimeIter> Iterator for TimeRangeBuilder<I> {
	type Item = TimeRange;
	fn next(&mut self) -> Option<Self::Item> {
		match self {
			TimeRangeBuilder::Inclusive(iter, queued_time) => {
				if let Some(time) = queued_time.take() {
					let t = TimeBound::Included(time);
					Some(TimeRange(t, t))
				} else if let Some(time) = iter.next() {
					 // Ignore Repeated Times:
					for next_time in iter.by_ref() {
						if time != next_time {
							*queued_time = Some(next_time);
							break
						}
					}
					
					let t = TimeBound::Included(time);
					Some(TimeRange(t, t))
				} else {
					None
				}
			},
			TimeRangeBuilder::Exclusive(iter) => {
				if let Some(a) = iter.next().map(TimeBound::Excluded) {
					let b = iter.next().map(TimeBound::Excluded)
						.unwrap_or(TimeBound::Unbounded);
					
					if a == b {
						self.next()
					} else {
						Some(TimeRange(a, b))
					}
				} else {
					None
				}
			},
			TimeRangeBuilder::Unbounded(iter) => {
				let mut iter = iter.take().expect("this should always exist");
				let next = iter.next().map(TimeBound::Excluded)
					.unwrap_or(TimeBound::Unbounded);
				
				*self = TimeRangeBuilder::Exclusive(iter);
				
				Some(TimeRange(TimeBound::Unbounded, next))
			},
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		match self {
			TimeRangeBuilder::Inclusive(iter, queued_time) => {
				if queued_time.is_some() {
					(1, iter.size_hint().1.and_then(|x| x.checked_add(1)))
				} else {
					(0, iter.size_hint().1)
				}
			},
			TimeRangeBuilder::Exclusive(iter) |
			TimeRangeBuilder::Unbounded(Some(iter)) => {
				(0, iter.size_hint().1.map(|x| 1 + x/2))
			},
			TimeRangeBuilder::Unbounded(None) => unreachable!(),
		}
	}
}

/// Filters times for use with [`TimeRanges::into_filtered`].
#[derive(Clone)]
pub(crate) struct TimeFilter<I, F> {
	times: I,
	filter: F,
}

impl<I: TimeRangeIter, F> Iterator for TimeFilter<I, F>
where
	F: crate::kind::RootFilterMap + 'static
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

/// Iterator of [`TimeRange`]s.
#[must_use]
#[derive(Clone)]
pub struct TimeRanges<I> {
	times: I
}

impl TimeRanges<std::iter::Empty<TimeRange>> {
	pub fn empty() -> Self {
		Self {
			times: std::iter::empty()
		}
	}
}

impl TimeRanges<std::iter::Once<TimeRange>> {
	pub fn from_range(range: impl RangeBounds<Time>) -> Self {
		Self::try_from_range(range)
			.expect("must be true: lower bound <= upper bound")
	}
	
	pub fn try_from_range(range: impl RangeBounds<Time>) -> Option<Self> {
		let range = TimeRange(
			range.start_bound().cloned(),
			range.end_bound().cloned()
		);
		if range.cmp_lower_to_upper(&range).is_gt() {
			return None
		}
		Some(Self {
			times: std::iter::once(range)
		})
	}
}

impl<I: TimeIter> TimeRanges<TimeRangeBuilder<I>> {
	pub(crate) fn new(
		iter: impl IntoIterator<IntoIter=I>,
		initial_order: Ordering,
		order: Ordering,
	) -> TimeRanges<TimeRangeBuilder<I>>
	{
		let iter = iter.into_iter();
		TimeRanges {
			times: if order == initial_order {
				TimeRangeBuilder::Unbounded(Some(iter))
			} else if order.is_eq() {
				TimeRangeBuilder::Inclusive(iter, None)
			} else {
				TimeRangeBuilder::Exclusive(iter)
			}
		}
	}
}

impl<I: TimeRangeIter> TimeRanges<I> {
	pub(crate) fn into_filtered<F>(self, f: F) -> TimeRanges<TimeFilter<I, F>>
	where
		F: crate::kind::RootFilterMap + 'static
	{
		TimeRanges {
			times: TimeFilter {
				times: self.times,
				filter: f,
			}
		}
	}
	
	/// Decrements the lower bound of each range by 1 nanosecond.  
	pub fn pre(self) -> TimeRanges<impl TimeRangeIter> {
		self.into_filtered(|t: Time, is_end: bool| if is_end {
			Some(t)
		} else {
			t.checked_sub(NANOSEC)
		})
	}
}

impl From<Time> for TimeRanges<std::iter::Once<TimeRange>> {
	fn from(value: Time) -> Self {
		TimeRanges {
			times: std::iter::once(TimeRange(
				TimeBound::Included(value),
				TimeBound::Included(value)
			))
		}
	}
}

impl Default for TimeRanges<std::iter::Empty<TimeRange>> {
	fn default() -> Self {
		Self::empty()
	}
}

impl<I: TimeRangeIter> Iterator for TimeRanges<I> {
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

impl<A: TimeRangeIter, B: TimeRangeIter> BitAnd<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<InterTimeRanges<A, B>>;
	
	/// Intersection of ranges.
	fn bitand(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges {
			times: InterTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times)
			}
		}
	}
}

impl<A: TimeRangeIter, B: TimeRangeIter> BitOr<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<UnionTimeRanges<A, B>>;
	
	/// Union of ranges.
	fn bitor(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges {
			times: UnionTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times)
			}
		}
	}
}

impl<A: TimeRangeIter, B: TimeRangeIter> BitXor<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<DiffTimeRanges<A, B>>;
	
	/// Symmetric difference of ranges.
	fn bitxor(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges {
			times: DiffTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times),
				range: None,
			}
		}
	}
}

impl<I: TimeRangeIter> Not for TimeRanges<I> {
	type Output = TimeRanges<InvTimes<I>>;
	
	/// Inverse of ranges.
	fn not(self) -> Self::Output {
		TimeRanges {
			times: InvTimes {
				iter: self.times,
				prev: None,
			}
		}
	}
}

/// [`Not`] implementation for [`TimeRanges`].
#[derive(Clone)]
pub struct InvTimes<I> {
	iter: I,
	prev: Option<TimeBound>,
}

impl<I: TimeRangeIter> Iterator for InvTimes<I> {
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

/// [`BitAnd`] implementation for [`TimesRanges`].
#[derive(Clone)]
pub struct InterTimeRanges<A, B> {
	iter: OrdTimes<A, B>,
}

impl<A: TimeRangeIter, B: TimeRangeIter> Iterator for InterTimeRanges<A, B> {
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

/// [`BitOr`] implementation for [`TimesRanges`].
#[derive(Clone)]
pub struct UnionTimeRanges<A, B> {
	iter: OrdTimes<A, B>,
}

impl<A: TimeRangeIter, B: TimeRangeIter> Iterator for UnionTimeRanges<A, B> {
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

/// [`BitXor`] implementation for [`TimesRanges`].
#[derive(Clone)]
pub struct DiffTimeRanges<A, B> {
	iter: OrdTimes<A, B>,
	range: Option<TimeRange>,
}

impl<A: TimeRangeIter, B: TimeRangeIter> Iterator for DiffTimeRanges<A, B> {
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

/// Orders two [`TimeRanges`] iterators in parallel.
#[derive(Clone)]
struct OrdTimes<A, B> {
	a_iter: A,
	a_next: Option<TimeRange>,
	b_iter: B,
	b_next: Option<TimeRange>,
}

impl<A: TimeRangeIter, B: TimeRangeIter> OrdTimes<A, B> {
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

impl<A: TimeRangeIter, B: TimeRangeIter> Iterator for OrdTimes<A, B> {
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

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn range_logic() {
		assert!(TimeRanges::try_from_range(5*NANOSEC..4*NANOSEC).is_none());
		assert!(TimeRanges::try_from_range(5*NANOSEC..5*NANOSEC).is_some());
		let t = SEC;
		let a = TimeRanges::new([2*t, 3*t, 10*t, 20*t, 40*t, 40*t, 40*t, 40*t + NANOSEC], Ordering::Less, Ordering::Less);
		let b = TimeRanges::new([2*t, 5*t, 20*t, 40*t, 50*t, 50*t, 50*t], Ordering::Less, Ordering::Greater);
		let c = TimeRanges::new([0*t, 7*t, 20*t, 20*t, 50*t - NANOSEC, 50*t, 50*t + NANOSEC], Ordering::Greater, Ordering::Equal);
		assert_eq!(Vec::from_iter(a.clone()), [
			(Time::ZERO, 2*t - NANOSEC),
			(3*t + NANOSEC, 10*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(b.clone()), [
			(2*t + NANOSEC, 5*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(50*t + NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(c.clone()), [
			(0*t, 0*t),
			(7*t, 7*t),
			(20*t, 20*t),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, 50*t),
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(a.clone() & b.clone()), [
			(3*t + NANOSEC, 5*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(50*t + NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(a.clone() | b.clone()), [
			(0*t, 2*t - NANOSEC),
			(2*t + NANOSEC, 10*t - NANOSEC),
			(20*t + NANOSEC, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX)
		]);
		assert_eq!(Vec::from_iter(a.clone() ^ b.clone()), [
			(0*t, 2*t - NANOSEC),
			(2*t + NANOSEC, 3*t),
			(5*t, 10*t - NANOSEC),
			(40*t + 2*NANOSEC, 50*t)
		]);
		assert_eq!(Vec::from_iter(a.clone() & a.clone()), Vec::from_iter(a.clone()));
		assert_eq!(Vec::from_iter(b.clone() & b.clone()), Vec::from_iter(b.clone()));
		assert_eq!(Vec::from_iter(c.clone() & c.clone()), Vec::from_iter(c.clone()));
		assert_eq!(Vec::from_iter(a.clone() | a.clone()), Vec::from_iter(a.clone()));
		assert_eq!(Vec::from_iter(b.clone() | b.clone()), Vec::from_iter(b.clone()));
		assert_eq!(Vec::from_iter(c.clone() | c.clone()), Vec::from_iter(c.clone()));
		assert_eq!(Vec::from_iter(a.clone() ^ a.clone()), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(b.clone() ^ b.clone()), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(c.clone() ^ c.clone()), Vec::from_iter([]));
		assert_eq!(Vec::from_iter(!!!a.clone()), [
			(2*t, 3*t),
			(10*t, 20*t),
			(40*t, 40*t + NANOSEC)
		]);
		assert_eq!(Vec::from_iter(!b.clone()), [
			(Time::ZERO, 2*t),
			(5*t, 20*t),
			(40*t, 50*t)
		]);
		assert_eq!(Vec::from_iter(c.clone() & a.clone()), [
			(0*t, 0*t),
			(7*t, 7*t),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, 50*t),
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(b.clone() & c.clone()), [
			(50*t + NANOSEC, 50*t + NANOSEC),
		]);
		assert_eq!(Vec::from_iter(c.clone() | a.clone()), [
			(Time::ZERO, 2*t - NANOSEC),
			(3*t + NANOSEC, 10*t - NANOSEC),
			(20*t, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(b.clone() | c.clone()), [
			(0*t, 0*t),
			(2*t + NANOSEC, 5*t - NANOSEC),
			(7*t, 7*t),
			(20*t, 40*t - NANOSEC),
			(50*t - NANOSEC, 50*t - NANOSEC),
			(50*t, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(c.clone() ^ a.clone()), [
			(0*t + NANOSEC, 2*t - NANOSEC),
			(3*t + NANOSEC, 7*t - NANOSEC),
			(7*t + NANOSEC, 10*t - NANOSEC),
			(20*t, 40*t - NANOSEC),
			(40*t + 2*NANOSEC, 50*t - 2*NANOSEC),
			(50*t + 2*NANOSEC, Time::MAX),
		]);
		assert_eq!(Vec::from_iter(b.clone() ^ c.clone()), [
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