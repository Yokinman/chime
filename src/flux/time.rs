//! Working with time.

use std::ops::{BitAnd, BitOr, BitXor, Bound, Not, RangeBounds};
use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

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

/// Iterator trait object for [`Times`].
trait TimeIter: Iterator<Item=Time> + Send + Sync {
	fn clone_into_box(&self) -> Box<dyn TimeIter<Item=Time>>;
}
impl<T: Iterator<Item=Time> + Send + Sync + Clone + 'static> TimeIter for T {
	fn clone_into_box(&self) -> Box<dyn TimeIter<Item=Time>> {
		Box::new(self.clone())
	}
}

/// Iterator of [`Time`] values.
#[must_use]
pub(crate) struct Times {
	heap: BinaryHeap<Reverse<Time>>,
	iter: Box<dyn TimeIter<Item=Time>>,
	// !!! This could use generics, but it adds way too many tedious bounds
	// without the ability to use `impl Trait` syntax in traits. Also for my
	// purposes I would have to box the `Times<T>` regardless.
}

impl Times {
	pub fn new<T>(iter: T) -> Self
	where
		T: IntoIterator<Item=Time>,
		<T as IntoIterator>::IntoIter: Send + Sync + Clone + 'static,
	{
		Self {
			heap: BinaryHeap::new(),
			iter: Box::new(iter.into_iter()), // < This takes like 1 microsecond
		}
	}
	
	fn fill_heap(&mut self, basis: Time, order: &mut Ordering) {
		//! Consumes the iterator into a binary heap if it has a known upper
		//! limit. Otherwise, it just maintains a minimum heap size.
		
		let size = 256.min(match self.iter.size_hint() {
			(_, Some(size)) => size,
			(size, None) => (size + 4).saturating_sub(self.heap.len()),
		});
		
		if size != 0 {
			self.heap.reserve(size);
			for _ in 0..size {
				if let Some(t) = self.iter.next() {
					self.heap.push(Reverse(t));
					if t < basis {
						*order = order.reverse();
					}
				} else {
					break
				}
			}
		}
	}
	
	fn pop(&mut self) -> Option<Time> {
		if self.heap.is_empty() {
			self.fill_heap(Time::ZERO, &mut Ordering::Equal);
		}
		let time = self.heap.pop().map(|Reverse(t)| t);
		self.fill_heap(Time::ZERO, &mut Ordering::Equal);
		time
	}
	
	pub(crate) fn peek(&mut self) -> Option<Time> {
		if self.heap.is_empty() {
			self.fill_heap(Time::ZERO, &mut Ordering::Equal);
		}
		self.heap.peek().map(|&Reverse(t)| t)
	}
	
	fn into_ranges(self) -> TimeRanges {
		TimeRanges::new(self, Time::ZERO, Ordering::Greater, Ordering::Equal)
	}
}

impl Clone for Times {
	fn clone(&self) -> Self {
		Self {
			heap: self.heap.clone(),
			iter: self.iter.clone_into_box(),
		}
	}
}

impl Default for Times {
	fn default() -> Self {
		Self::new(std::iter::empty())
	}
}

impl Iterator for Times {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		if let Some(a) = self.pop() {
			while let Some(b) = self.peek() {
				if b > a {
					break
				}
				self.pop();
			}
			return Some(a)
		}
		None
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (_, upper) = self.iter.size_hint();
		(0, upper.map(|x| x + self.heap.len()))
	}
	fn count(self) -> usize {
		self.iter.count() + self.heap.len()
	}
}

impl BitAnd for Times {
	type Output = Self;
	
	/// Intersection of times.
	fn bitand(self, rhs: Self) -> Self::Output {
		Times::new(OrdTimes::new(self, rhs)
			.filter_map(|t| match t {
				OrdTime::Equal(t) => Some(t),
				_ => None
			}))
	}
}

impl BitOr for Times {
	type Output = Self;
	
	/// Union of times.
	fn bitor(self, rhs: Self) -> Self::Output {
		Times::new(OrdTimes::new(self, rhs)
			.map(|t| match t {
				OrdTime::Less(t)    |
				OrdTime::Greater(t) |
				OrdTime::Equal(t)   => t
			}))
	}
}

impl BitXor for Times {
	type Output = Self;
	
	/// Symmetric difference of times.
	fn bitxor(self, rhs: Self) -> Self::Output {
		Times::new(OrdTimes::new(self, rhs)
			.filter_map(|t| match t {
				OrdTime::Less(t) |
				OrdTime::Greater(t) => Some(t),
				_ => None
			}))
	}
}

impl Not for Times {
	type Output = TimeRanges;
	
	/// Inverse of times.
	fn not(self) -> Self::Output {
		!self.into_ranges()
	}
}

impl BitAnd<TimeRanges> for Times {
	type Output = Self;
	
	/// Intersection of times & ranges.
	fn bitand(self, rhs: TimeRanges) -> Self::Output {
		Times::new((self.into_ranges() & rhs)
			.map(|(a, b)| {
				debug_assert_eq!(a, b);
				a
			}))
	}
}

/// Iterator of [`Time`] ranges.
#[must_use]
#[derive(Clone)]
pub struct TimeRanges {
	times: Times,
	order: Ordering,
	basis: Option<Time>,
}

impl TimeRanges {
	pub fn empty() -> Self {
		Self::new(
			std::iter::empty(),
			Time::ZERO,
			Ordering::Greater,
			Ordering::Equal
		)
	}
	
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
				TimeRanges::new([a, b], Time::ZERO, Ordering::Greater, Ordering::Less)
			},
			(Bound::Excluded(a), Bound::Unbounded) => {
				TimeRanges::new([a], Time::ZERO, Ordering::Greater, Ordering::Less)
			},
			(Bound::Unbounded, Bound::Excluded(b)) => {
				TimeRanges::new([b], Time::ZERO, Ordering::Greater, Ordering::Greater)
			},
			(Bound::Unbounded, Bound::Unbounded) => {
				TimeRanges::new([], Time::ZERO, Ordering::Greater, Ordering::Greater)
			},
			_ => unreachable!()
		})
	}
	
	pub(crate) fn new<T>(iter: T, basis: Time, basis_order: Ordering, order: Ordering)
		-> Self
	where
		T: IntoIterator<Item=Time>,
		<T as IntoIterator>::IntoIter: Send + Sync + Clone + 'static,
	{
		let mut times = if basis_order == Ordering::Equal {
			Times::default()
		} else {
			Times::new(iter)
		};
		
		let mut order = if order == basis_order {
			Ordering::Greater
		} else if order == Ordering::Equal {
			Ordering::Equal
		} else {
			Ordering::Less
		};
		
		 // Find True Initial Order:
		if basis != Time::ZERO {
			debug_assert!(times.heap.is_empty()); 
			times.fill_heap(basis, &mut order);
		}
		
		Self {
			times,
			order,
			basis: Some(Time::ZERO),
		}
	}
	
	fn into_not_equal(mut self) -> TimeRanges {
		//! Converts points to ranges.
		if self.order == Ordering::Equal {
			if let Some(a) = self.times.peek() {
				let times = self
					.flat_map(|(a, b)| [
						a.checked_sub(NANOSEC),
						b.checked_add(NANOSEC)
					])
					.flatten();
				
				let order = if a == Time::ZERO {
					Ordering::Less
				} else {
					Ordering::Greater
				};
				
				return Self::new(times, Time::ZERO, Ordering::Less, order)
			}
			self.order = Ordering::Less;
		}
		self
	}
}

impl From<Time> for TimeRanges {
	fn from(value: Time) -> Self {
		TimeRanges::new([value], Time::ZERO, Ordering::Greater, Ordering::Equal)
	}
}

impl Default for TimeRanges {
	fn default() -> Self {
		Self::empty()
	}
}

impl Iterator for TimeRanges {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		let basis = self.basis?;
		let range = match self.order {
			Ordering::Equal => {
				if let Some(a) = self.times.next() {
					let mut b = a;
					
					 // Combine Adjacent Points:
					while b != Time::MAX && Some(b + NANOSEC) == self.times.peek() {
						b = self.times.next().unwrap();
					}
					
					Some((a, b))
				} else {
					None
				}
			},
			Ordering::Greater => {
				self.order = Ordering::Less;
				if let Some(b) = self.times.pop() {
					if b == Time::ZERO {
						self.next()
					} else {
						Some((Time::ZERO, b - NANOSEC))
					}
				} else {
					Some((Time::ZERO, Time::MAX))
				}
			},
			Ordering::Less => {
				let mut range = None;
				while let Some(mut a) = self.times.pop() {
					a = a.max(basis);
					range = if let Some(mut b) = self.times.pop() {
						 // Belated Order Correction:
						if b < a {
							let t = a;
							a = b.max(basis);
							b = t;
						}
						
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
			}
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

impl BitAnd for TimeRanges {
	type Output = Self;
	
	/// Intersection of ranges.
	fn bitand(self, rhs: Self) -> Self::Output {
		match (self.order, rhs.order) {
			(Ordering::Equal, Ordering::Equal) => TimeRanges::new(
				self.times & rhs.times,
				Time::ZERO,
				Ordering::Greater,
				Ordering::Equal,
			),
			
			(Ordering::Equal, _) => self.into_not_equal() & rhs,
			(_, Ordering::Equal) => self & rhs.into_not_equal(),
			
			(a_order, b_order) => Self::new(
				JoinedTimeRanges {
					iter: OrdTimes::new(self.times, rhs.times),
					flag: ((a_order == Ordering::Greater) as u8)
						| (((b_order == Ordering::Greater) as u8) << 1),
					has_extra: false,
				},
				Time::ZERO,
				Ordering::Greater,
				a_order.min(b_order),
			)
		}
	}
}

impl BitOr for TimeRanges {
	type Output = Self;
	
	/// Union of ranges.
	fn bitor(self, rhs: Self) -> Self::Output {
		match (self.order, rhs.order) {
			(Ordering::Equal, Ordering::Equal) => TimeRanges::new(
				self.times | rhs.times,
				Time::ZERO,
				Ordering::Greater,
				Ordering::Equal,
			),
			
			(Ordering::Equal, _) => self.into_not_equal() | rhs,
			(_, Ordering::Equal) => self | rhs.into_not_equal(),
			
			(a_order, b_order) => TimeRanges::new(
				JoinedTimeRanges {
					iter: OrdTimes::new(self.times, rhs.times),
					flag: !(((a_order == Ordering::Greater) as u8)
						| (((b_order == Ordering::Greater) as u8) << 1)),
					has_extra: false,
				},
				Time::ZERO,
				Ordering::Greater,
				a_order.max(b_order),
			)
		}
	}
}

impl BitXor for TimeRanges {
	type Output = Self;
	
	/// Intersection of ranges.
	fn bitxor(self, rhs: Self) -> Self::Output {
		match (self.order, rhs.order) {
			(Ordering::Equal, Ordering::Equal) => TimeRanges::new(
				self.times ^ rhs.times,
				Time::ZERO,
				Ordering::Greater,
				Ordering::Equal,
			),
			
			(Ordering::Equal, _) => self.into_not_equal() ^ rhs,
			(_, Ordering::Equal) => self ^ rhs.into_not_equal(),
			
			(a_order, b_order) => TimeRanges::new(
				SplitTimeRanges {
					iter: OrdTimes::new(self.times, rhs.times),
					flag: ((a_order == Ordering::Greater) as u8) | (((b_order == Ordering::Greater) as u8) << 1),
					extra: None,
				},
				Time::ZERO,
				Ordering::Greater,
				if a_order == b_order {
					Ordering::Less
				} else {
					Ordering::Greater
				}
			)
		}
	}
}

impl Not for TimeRanges {
	type Output = Self;
	
	/// Inverse of ranges.
	fn not(self) -> Self::Output {
		if self.order == Ordering::Equal {
			return !self.into_not_equal()
		}
		TimeRanges::new(
			InvTimes {
				iter: self.times,
				flag: self.order == Ordering::Greater,
				extra: None,
				prev: Time::ZERO,
			},
			Time::ZERO,
			Ordering::Greater,
			self.order.reverse(),
		)
	}
}

/// [`Not`] implementation for [`Times`] and [`TimeRanges`].
#[derive(Clone)]
struct InvTimes {
	iter: Times,
	flag: bool,
	extra: Option<Time>,
	prev: Time,
}

impl Iterator for InvTimes {
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
struct JoinedTimeRanges {
	iter: OrdTimes,
	flag: u8,
	has_extra: bool,
}

impl JoinedTimeRanges {
	fn is_union(&self) -> bool {
		//! As opposed to intersection.
		self.flag >> 7 != 0
	}
}

impl Iterator for JoinedTimeRanges {
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
struct SplitTimeRanges {
	iter: OrdTimes,
	flag: u8,
	extra: Option<Time>,
}

impl SplitTimeRanges {
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

impl Iterator for SplitTimeRanges {
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
struct OrdTimes {
	a_iter: Times,
	b_iter: Times,
	prev: Time,
}

impl OrdTimes {
	fn new(a_iter: Times, b_iter: Times) -> OrdTimes {
		OrdTimes {
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

impl Iterator for OrdTimes {
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
	
	#[test]
	fn time_logic() {
		let t = SEC;
		let a = Times::new([1*t, 2*t, 3*t]);
		let b = Times::new([2*t, 5*t]);
		assert_eq!(Vec::from_iter(a.clone() & b.clone()), [2*t]);
		assert_eq!(Vec::from_iter(a.clone() | b.clone()), [1*t, 2*t, 3*t, 5*t]);
		assert_eq!(Vec::from_iter(a.clone() ^ b.clone()), [1*t, 3*t, 5*t]);
		assert_eq!(Vec::from_iter(!!!b), [
			(Time::ZERO, 2*t - NANOSEC),
			(2*t + NANOSEC, 5*t - NANOSEC),
			(5*t + NANOSEC, Time::MAX)
		]);
		assert_eq!(
			Vec::from_iter(a.clone() & TimeRanges::new([], Time::ZERO, Ordering::Greater, Ordering::Greater)),
			Vec::from_iter(a.clone())
		);
	}
	
	#[test]
	fn range_logic() {
		let t = SEC;
		let a = TimeRanges::new([2*t, 3*t, 10*t, 20*t, 40*t, 40*t, 40*t, 40*t + NANOSEC], Time::ZERO, Ordering::Less, Ordering::Less);
		let b = TimeRanges::new([2*t, 5*t, 20*t, 40*t, 50*t, 50*t, 50*t], Time::ZERO, Ordering::Less, Ordering::Greater);
		let c = TimeRanges::new([0*t, 7*t, 20*t, 20*t, 50*t - NANOSEC, 50*t, 50*t + NANOSEC], Time::ZERO, Ordering::Greater, Ordering::Equal);
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
	
	#[test]
	fn order_fix() {
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
			Time::ZERO,
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
			Time::ZERO,
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