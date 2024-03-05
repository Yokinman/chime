//! Working with time.

use std::ops::{BitAnd, BitOr, BitXor, Bound, Not, RangeBounds};
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

/// Iterator types usable by [`Times`].
pub trait TimeIter: Iterator<Item=Time> + Send + Sync + Clone + 'static {}
impl<T: Iterator<Item=Time> + Send + Sync + Clone + 'static> TimeIter for T {}

/// Iterator of [`Time`] values.
#[must_use]
#[derive(Clone)]
pub(crate) struct Times<I> {
	iter: I,
	next_time: Option<Time>,
	
	/// Used to convert points of time into ranges.
	doubled_time: Option<Option<Time>>,
}

impl<I: TimeIter> Times<I> {
	pub fn new(iter: impl IntoIterator<IntoIter=I>) -> Self {
		Self::with_basis(iter)
	}
	
	pub fn with_basis(iter: impl IntoIterator<IntoIter=I>) -> Self {
		let mut times = Self {
			iter: iter.into_iter(),
			next_time: None,
			doubled_time: None,
		};
		times.pop();
		times
	}
	
	fn pop(&mut self) -> Option<Time> {
		if let Some((time, is_end)) = self.pop_doubled() {
			if is_end {
				time.checked_add(NANOSEC)
			} else {
				time.checked_sub(NANOSEC)
					.or_else(|| self.pop())
			}
		} else {
			std::mem::replace(&mut self.next_time, self.iter.next())
		}
	}
	
	fn pop_doubled(&mut self) -> Option<(Time, bool)> {
		if let Some(doubled_time) = self.doubled_time.take() {
			if let Some(doubled_time) = doubled_time {
				self.doubled_time = Some(None);
				return Some((doubled_time, true))
			}
			if let Some(time) = self.peek() {
				self.pop();
				
				 // Combine Repeated/Adjacent Times:
				let mut curr = time;
				while let Some(next) = self.peek() {
					if curr == next || curr.checked_add(NANOSEC) == Some(next) {
						curr = next;
						self.pop();
					} else {
						break
					}
				}
				self.doubled_time = Some(Some(curr));
				
				return Some((time, false))
			}
			self.doubled_time = Some(None);
		}
		None
	}
	
	pub(crate) fn peek(&self) -> Option<Time> {
		if let Some(doubled_time) = self.doubled_time {
			if let Some(doubled_time) = doubled_time {
				return doubled_time.checked_add(NANOSEC)
			} else if let Some(next_time) = self.next_time {
				return next_time.checked_sub(NANOSEC).or(Some(NANOSEC))
			}
		}
		self.next_time
	}
}

impl Default for Times<std::iter::Empty<Time>> {
	fn default() -> Self {
		Self::new(std::iter::empty())
	}
}

impl<I: TimeIter> Iterator for Times<I> {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(a) = self.pop() {
			while let Some(b) = self.peek() {
				if b > a {
					break
				}
				assert_eq!(a, b, "times must be ordered: {:?}", (a, b));
				self.pop();
			}
			return Some(a)
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		(0, self.iter.size_hint().1)
	}
}

/// Filters times for use with [`TimeRanges::into_filtered`].
#[derive(Clone)]
pub(crate) struct TimeFilter<I, F> {
	times: Times<I>,
	filter: F,
	is_end: bool,
}

impl<I: TimeIter, F> Iterator for TimeFilter<I, F>
where
	F: crate::kind::RootFilterMap<Output=Option<Time>> + 'static
{
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some((time, is_end)) = self.times.pop_doubled() {
			// !!! The whole concept of order for TimeRanges needs to be
			// reworked in some way to account for filters returning None.
			(self.filter)(time, is_end)
				.and_then(|t| if is_end {
					t.checked_add(NANOSEC)
				} else {
					t.checked_sub(NANOSEC)
						.or_else(|| self.next())
				})
		} else if let Some(time) = self.times.pop() {
			let is_end = self.is_end;
			self.is_end = !is_end;
			(self.filter)(time, is_end)
		} else {
			None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.times.size_hint()
	}
}

/// Iterator of [`Time`] ranges.
#[must_use]
#[derive(Clone)]
pub struct TimeRanges<I> {
	times: Times<I>,
	is_end: bool,
	basis: Option<Time>,
}

impl TimeRanges<std::iter::Empty<Time>> {
	pub fn empty() -> TimeRanges<std::iter::Empty<Time>> {
		Self::new(
			std::iter::empty(),
			Ordering::Greater,
			Ordering::Equal
		)
	}
}

impl TimeRanges<std::iter::Chain<std::option::IntoIter<Time>, std::option::IntoIter<Time>>> {
	pub fn from_range(range: impl RangeBounds<Time>) -> Self {
		Self::try_from_range(range)
			.expect("must be true: lower bound <= upper bound")
	}
	
	pub fn try_from_range(range: impl RangeBounds<Time>) -> Option<Self> {
		let a = match range.start_bound() {
			Bound::Included(&Time::ZERO) => Bound::Unbounded,
			Bound::Excluded(&Time::MAX) => return None,
			Bound::Included(&t) => Bound::Excluded(t - NANOSEC),
			Bound::Excluded(&t) => Bound::Excluded(t),
			Bound::Unbounded => Bound::Unbounded,
		};
		let b = match range.end_bound() {
			Bound::Included(&Time::MAX) => Bound::Unbounded,
			Bound::Excluded(&Time::ZERO) => return None,
			Bound::Included(&t) => Bound::Excluded(t + NANOSEC),
			Bound::Excluded(&t) => Bound::Excluded(t),
			Bound::Unbounded => Bound::Unbounded,
		};
		Some(match (a, b) {
			(Bound::Excluded(a), Bound::Excluded(b)) => {
				if a >= b || a == b - NANOSEC {
					return None
				}
				TimeRanges::new(
					Some(a).into_iter().chain(Some(b)),
					Ordering::Greater,
					Ordering::Less
				)
			},
			(Bound::Excluded(a), Bound::Unbounded) => TimeRanges::new(
				Some(a).into_iter().chain(None),
				Ordering::Greater,
				Ordering::Less
			),
			(Bound::Unbounded, Bound::Excluded(b)) => TimeRanges::new(
				Some(b).into_iter().chain(None),
				Ordering::Greater,
				Ordering::Greater
			),
			(Bound::Unbounded, Bound::Unbounded) => TimeRanges::new(
				None.into_iter().chain(None),
				Ordering::Greater,
				Ordering::Greater
			),
			_ => unreachable!()
		})
	}
}

impl<I: TimeIter> TimeRanges<I> {
	pub(crate) fn new(
		iter: impl IntoIterator<IntoIter=I>,
		initial_order: Ordering,
		border: Ordering,
	) -> TimeRanges<I>
	{
		let mut times = Times::with_basis(iter);
		let mut is_end = false;
		if border == initial_order {
			is_end = true;
		} else if border == Ordering::Equal {
			if times.next_time == Some(Time::ZERO) {
				is_end = true;
			}
			times.doubled_time = Some(None);
		}
		TimeRanges {
			times,
			is_end,
			basis: Some(Time::ZERO),
		}
	}
	
	pub(crate) fn into_filtered<F>(self, f: F) -> TimeRanges<TimeFilter<I, F>>
	where
		F: crate::kind::RootFilterMap<Output=Option<Time>> + 'static
	{
		TimeRanges {
			times: Times::new(TimeFilter {
				times: self.times,
				filter: f,
				is_end: self.is_end,
			}),
			is_end: self.is_end,
			basis: self.basis,
		}
	}
	
	/// Decrements the lower bound of each range by 1 nanosecond.  
	pub fn pre(self) -> TimeRanges<impl TimeIter> {
		self.into_filtered(|t, is_end| if is_end {
			Some(t)
		} else {
			t.checked_sub(NANOSEC)
		})
	}
}

impl From<Time> for TimeRanges<std::iter::Once<Time>> {
	fn from(value: Time) -> Self {
		TimeRanges::new(
			std::iter::once(value),
			Ordering::Greater,
			Ordering::Equal
		)
	}
}

impl Default for TimeRanges<std::iter::Empty<Time>> {
	fn default() -> Self {
		Self::empty()
	}
}

impl<I: TimeIter> Iterator for TimeRanges<I> {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		let basis = self.basis?;
		let range = if self.is_end {
			self.is_end = false;
			if let Some(b) = self.times.pop() {
				assert!(basis <= b, "times must be ordered: {:?}", (basis, b));
				if b == Time::ZERO {
					self.next()
				} else {
					Some((Time::ZERO, b - NANOSEC))
				}
			} else {
				Some((Time::ZERO, Time::MAX))
			}
		} else {
			let mut range = None;
			while let Some(a) = self.times.pop() {
				assert!(basis <= a, "times must be ordered: {:?}", (basis, a));
				range = if let Some(b) = self.times.pop() {
					assert!(a <= b, "times must be ordered: {:?}", (a, b));
					
					 // Ignore Zero-sized Ranges:
					if a == b || a == b - NANOSEC {
						continue
					}
					
					Some((a + NANOSEC, b - NANOSEC))
				} else if a == Time::MAX {
					None
				} else {
					Some((a + NANOSEC, Time::MAX))
				};
				break
			}
			range
		};
		
		if let Some((a, b)) = range {
			self.basis = b.checked_add(NANOSEC);
			Some((a, b))
		} else {
			None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, upper) = self.times.size_hint();
		(0, upper.map(|x| 1 + x/2))
	}
}

impl<A: TimeIter, B: TimeIter> BitAnd<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<JoinedTimeRanges<A, B>>;
	
	/// Intersection of ranges.
	fn bitand(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges::new(
			JoinedTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times),
				flag: (self.is_end as u8) | ((rhs.is_end as u8) << 1),
				has_extra: false,
			},
			Ordering::Greater,
			if self.is_end & rhs.is_end {
				Ordering::Greater
			} else {
				Ordering::Less
			},
		)
	}
}

impl<A: TimeIter, B: TimeIter> BitOr<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<JoinedTimeRanges<A, B>>;
	
	/// Union of ranges.
	fn bitor(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges::new(
			JoinedTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times),
				flag: !((self.is_end as u8) | ((rhs.is_end as u8) << 1)),
				has_extra: false,
			},
			Ordering::Greater,
			if self.is_end | rhs.is_end {
				Ordering::Greater
			} else {
				Ordering::Less
			},
		)
	}
}

impl<A: TimeIter, B: TimeIter> BitXor<TimeRanges<B>> for TimeRanges<A> {
	type Output = TimeRanges<SplitTimeRanges<A, B>>;
	
	/// Intersection of ranges.
	fn bitxor(self, rhs: TimeRanges<B>) -> Self::Output {
		TimeRanges::new(
			SplitTimeRanges {
				iter: OrdTimes::new(self.times, rhs.times),
				flag: (self.is_end as u8) | ((rhs.is_end as u8) << 1),
				extra: None,
			},
			Ordering::Greater,
			if self.is_end ^ rhs.is_end {
				Ordering::Greater
			} else {
				Ordering::Less
			},
		)
	}
}

impl<I: TimeIter> Not for TimeRanges<I> {
	type Output = TimeRanges<InvTimes<I>>;
	
	/// Inverse of ranges.
	fn not(self) -> Self::Output {
		TimeRanges::new(
			InvTimes {
				iter: self.times,
				flag: self.is_end,
				extra: None,
				prev: Time::ZERO,
			},
			Ordering::Greater,
			if !self.is_end {
				Ordering::Greater
			} else {
				Ordering::Less
			},
		)
	}
}

/// [`Not`] implementation for [`Times`] and [`TimeRanges`].
#[derive(Clone)]
pub struct InvTimes<I> {
	iter: Times<I>,
	flag: bool,
	extra: Option<Time>,
	prev: Time,
}

impl<I: TimeIter> Iterator for InvTimes<I> {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if self.extra.is_some() {
			return std::mem::take(&mut self.extra)
		}
		while let Some(mut time) = self.iter.pop() {
			self.flag = !self.flag;
			time = time.max(self.prev);
			
			 // Adjacent Roots:
			if let Some(next) = self.iter.peek() {
				if next == time || next == time + NANOSEC {
					if self.flag {
						self.flag = !self.flag;
						self.iter.pop();
					} else if self.extra.is_none() {
						self.extra = time.checked_sub(NANOSEC);
					}
					continue
				}
			}
			
			self.prev = time;
			
			let time = if self.flag {
				if self.extra.is_some() {
					std::mem::replace(&mut self.extra, time.checked_add(NANOSEC))
				} else {
					time.checked_add(NANOSEC)
				}
			} else if self.extra.is_some() {
				std::mem::take(&mut self.extra)
			} else {
				time.checked_sub(NANOSEC)
			};
			
			if let Some(time) = time {
				return Some(time)
			} else {
				continue
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, max) = self.iter.size_hint();
		let extra_num = self.extra.is_some() as usize;
		(
			extra_num,
			max.map(|x| x + extra_num)
		)
	}
}

/// [`BitAnd`] and [`BitOr`] implementation for [`TimesRanges`].
#[derive(Clone)]
pub struct JoinedTimeRanges<A, B> {
	iter: OrdTimes<A, B>,
	flag: u8,
	has_extra: bool,
}

impl<A, B> JoinedTimeRanges<A, B> {
	fn is_union(&self) -> bool {
		//! As opposed to intersection.
		self.flag >> 7 != 0
	}
}

impl<A: TimeIter, B: TimeIter> Iterator for JoinedTimeRanges<A, B> {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if self.has_extra {
			self.has_extra = false;
			return Some(self.iter.prev())
		}
		while let Some(next) = self.iter.next() {
			match next {
				OrdTime::Less(time) => {
					self.flag ^= 0b01;
					if self.flag & 0b10 != 0 {
						return Some(time)
					}
				},
				OrdTime::Greater(time) => {
					self.flag ^= 0b10;
					if self.flag & 0b01 != 0 {
						return Some(time)
					}
				},
				OrdTime::Equal(time) => {
					self.flag ^= 0b11;
					if self.flag & 0b11 == 0b00 || self.flag & 0b11 == 0b11 {
						return Some(time)
					} else if self.is_union() {
						self.has_extra = true;
						return Some(time)
					}
				},
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, max) = self.iter.size_hint();
		let extra_num = self.has_extra as usize;
		(
			extra_num,
			max.map(|x| x + extra_num),
		)
	}
}

/// [`BitXor`] implementation for [`TimesRanges`].
#[derive(Clone)]
pub struct SplitTimeRanges<A, B> {
	iter: OrdTimes<A, B>,
	flag: u8,
	extra: Option<Time>,
}

impl<A: TimeIter, B: TimeIter> SplitTimeRanges<A, B> {
	fn step(&mut self, flag: u8, a: Time) -> Option<Time> {
		self.flag ^= flag;
		
		 // Equal:
		if flag == 0b11 {
			return match self.flag {
				0b00 => None,
				0b11 => match self.iter.peek() {
					Some(OrdTime::Less(time)) if time == a => {
						self.flag ^= 0b01;
						self.iter.next();
						Some(time)
					},
					Some(OrdTime::Greater(time)) if time == a => {
						self.flag ^= 0b10;
						self.iter.next();
						Some(time)
					},
					Some(OrdTime::Equal(time)) if time == a => {
						self.flag ^= 0b11;
						self.iter.next();
						None
					},
					_ => None
				},
				_ => {
					self.extra = Some(a);
					Some(a)
				}
			}
		}
		
		 // Opposing Range isn't Active:
		if self.flag & !flag == 0 {
			return Some(a)
		}
		
		let b = if flag == 0b01 {
			self.iter.peek()
		} else {
			self.iter.peek().map(|t| t.reverse())
		};
		
		 // Adjacent Roots:
		if let Some(OrdTime::Less(next)) = b {
			if next == a || next == a + NANOSEC {
				if self.flag & flag != 0 {
					self.flag ^= flag;
					self.iter.next();
				} else if self.extra.is_none() {
					self.extra = a.checked_sub(NANOSEC);
				}
				return None
			}
		}
		
		if self.flag & flag != 0 {
			if a != Time::MAX && b == Some(OrdTime::Greater(a + NANOSEC)) {
				self.flag ^= !flag & 0b11;
				self.iter.next();
				std::mem::take(&mut self.extra)
			} else if self.extra.is_some() {
				std::mem::replace(&mut self.extra, a.checked_add(NANOSEC))
			} else {
				a.checked_add(NANOSEC)
			}
		} else if self.extra.is_some() {
			std::mem::take(&mut self.extra)
		} else {
			a.checked_sub(NANOSEC)
		}
	}
}

impl<A: TimeIter, B: TimeIter> Iterator for SplitTimeRanges<A, B> {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if self.extra.is_some() {
			return std::mem::take(&mut self.extra)
		}
		while let Some(time) = self.iter.next() {
			if let Some(time) = match time {
				OrdTime::Less(t)    => self.step(0b01, t),
				OrdTime::Greater(t) => self.step(0b10, t),
				OrdTime::Equal(t)   => self.step(0b11, t),
			} {
				return Some(time)
			}
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, max) = self.iter.size_hint();
		let extra_num = self.extra.is_some() as usize;
		(
			extra_num,
			max.map(|x| x + extra_num)
		)
	}
}

/// Ordering paired with a time.
#[derive(Copy, Clone, Debug, PartialEq)]
enum OrdTime {
	Less(Time),
	Greater(Time),
	Equal(Time),
}

impl OrdTime {
	fn reverse(self) -> Self {
		match self {
			OrdTime::Less(t) => OrdTime::Greater(t),
			OrdTime::Greater(t) => OrdTime::Less(t),
			t => t,
		}
	}
}

/// Orders two [`Times`] iterators in parallel.
#[derive(Clone)]
struct OrdTimes<A, B> {
	a_iter: Times<A>,
	b_iter: Times<B>,
	prev: Time,
}

impl<A: TimeIter, B: TimeIter> OrdTimes<A, B> {
	fn new(a_iter: Times<A>, b_iter: Times<B>) -> Self {
		Self {
			a_iter,
			b_iter,
			prev: Time::ZERO,
		}
	}
	
	fn peek(&mut self) -> Option<OrdTime> {
		match (self.a_iter.peek(), self.b_iter.peek()) {
			(Some(mut a), Some(mut b)) => {
				a = a.max(self.prev);
				b = b.max(self.prev);
				match a.cmp(&b) {
					Ordering::Less    => Some(OrdTime::Less(a)),
					Ordering::Greater => Some(OrdTime::Greater(b)),
					Ordering::Equal   => Some(OrdTime::Equal(a)),
				}
			},
			(Some(a), _) => Some(OrdTime::Less(a.max(self.prev))),
			(_, Some(b)) => Some(OrdTime::Greater(b.max(self.prev))),
			_ => None
		}
	}
	
	fn prev(&self) -> Time {
		self.prev
	}
}

impl<A: TimeIter, B: TimeIter> Iterator for OrdTimes<A, B> {
	type Item = OrdTime;
	fn next(&mut self) -> Option<Self::Item> {
		Some(match self.peek()? {
			a @ OrdTime::Less(time) => {
				self.prev = time;
				self.a_iter.pop();
				a
			},
			b @ OrdTime::Greater(time) => {
				self.prev = time;
				self.b_iter.pop();
				b
			},
			ab @ OrdTime::Equal(time) => {
				self.prev = time;
				self.a_iter.pop();
				self.b_iter.pop();
				ab
			}
		})
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (a_min, a_max) = self.a_iter.size_hint();
		let (b_min, b_max) = self.b_iter.size_hint();
		(
			a_min.max(b_min),
			if let (Some(a_max), Some(b_max)) = (a_max, b_max) {
				Some(a_max + b_max)
			} else {
				None
			}
		)
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	
	// #[test]
	// fn time_logic() {
	// 	let t = SEC;
	// 	let a = Times::new([1*t, 2*t, 3*t]);
	// 	let b = Times::new([2*t, 5*t]);
	// 	assert_eq!(Vec::from_iter(a.clone() & b.clone()), [2*t]);
	// 	assert_eq!(Vec::from_iter(a.clone() | b.clone()), [1*t, 2*t, 3*t, 5*t]);
	// 	assert_eq!(Vec::from_iter(a.clone() ^ b.clone()), [1*t, 3*t, 5*t]);
	// 	assert_eq!(Vec::from_iter(!!!b), [
	// 		(Time::ZERO, 2*t - NANOSEC),
	// 		(2*t + NANOSEC, 5*t - NANOSEC),
	// 		(5*t + NANOSEC, Time::MAX)
	// 	]);
	// 	assert_eq!(
	// 		Vec::from_iter(a.clone() & TimeRanges::new([], Time::ZERO, Ordering::Greater, Ordering::Greater)),
	// 		Vec::from_iter(a.clone())
	// 	);
	// }
	
	#[test]
	fn range_logic() {
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
			(50*t - NANOSEC, 50*t + NANOSEC),
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
			(50*t - NANOSEC, 50*t + NANOSEC),
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
			(50*t - NANOSEC, Time::MAX),
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
			(50*t - NANOSEC, 50*t),
			(50*t + 2*NANOSEC, Time::MAX),
		]);
	}
	
	#[ignore]
	#[test]
	fn order_fix() {
		// !!! This should be rewritten to just validate panicking.
		
		#[derive(Clone)]
		struct RootIter(std::vec::IntoIter<Time>);
		
		impl Iterator for RootIter {
			type Item = Time;
			fn next(&mut self) -> Option<Self::Item> {
				self.0.next()
			}
			fn size_hint(&self) -> (usize, Option<usize>) {
				// No upper bound, so the time iterator assumes it's endless.
				(0, None)
			}
		}
		
		let t = SEC;
		let roots = vec![10*t, 20*t, 30*t, 40*t, 50*t, 7*t];
		let b_roots = vec![10*t, 25*t, 30*t, 45*t, 5*t];
		
		let times = Times::new(RootIter(roots.clone().into_iter()));
		assert_eq!(times.collect::<Vec<_>>(), [10*t, 20*t, 30*t, 40*t, 50*t]);
		
		let ranges = TimeRanges::new(
			RootIter(roots.into_iter()),
			Ordering::Greater,
			Ordering::Less,
		);
		assert_eq!(ranges.clone().collect::<Vec<_>>(), [
			// SHOULD be (7..10, 20..30, 40..50), but the iterator doesn't know
			// all the roots. So instead, it corrects itself along the way:
			(10*t + NANOSEC, 20*t - NANOSEC),
			(20*t + NANOSEC, 30*t - NANOSEC),
			(40*t + NANOSEC, 50*t - NANOSEC),
		]);
		
		let b_ranges = TimeRanges::new(
			RootIter(b_roots.into_iter()),
			Ordering::Greater,
			Ordering::Less,
		);
		assert_eq!(b_ranges.clone().collect::<Vec<_>>(), [
			// This one is within range of all the roots by chance, so it works.
			(5*t + NANOSEC, 10*t - NANOSEC),
			(25*t + NANOSEC, 30*t - NANOSEC),
			(45*t + NANOSEC, Time::MAX),
		]);
		
		 // Operations:
		assert_eq!((ranges.clone() & b_ranges.clone()).collect::<Vec<_>>(), [
			(25*t + NANOSEC, 30*t - NANOSEC),
			(45*t + NANOSEC, 50*t - NANOSEC),
		]);
		assert_eq!((ranges.clone() ^ b_ranges.clone()).collect::<Vec<_>>(), [
			(10*t + NANOSEC, 20*t - NANOSEC),
			(20*t + NANOSEC, 25*t),
			(40*t + NANOSEC, 45*t),
			(50*t, Time::MAX),
		]);
		assert_eq!((!ranges.clone()).collect::<Vec<_>>(), [
			(0*t, 10*t),
			(20*t, 20*t),
			(30*t, 40*t),
			(50*t, Time::MAX),
		]);
	}
}