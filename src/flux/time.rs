//! Working with time.

use std::ops::{BitAnd, BitOr, BitXor, Not};
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
	pub const MIN: Time = Time::from_secs(60);
	/// 1 hour (h).
	pub const HOUR: Time = Time::from_secs(60 * 60);
}
pub use units::*;

/// Iterator of [`Time`] values.
#[derive(Clone, Default)]
#[must_use]
pub struct Times(std::vec::IntoIter<Time>);

impl Iterator for Times {
	type Item = Time;
	fn next(&mut self) -> Option<Self::Item> {
		self.0.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.0.size_hint()
	}
	fn count(self) -> usize {
		self.0.count()
	}
}

impl FromIterator<Time> for Times {
	fn from_iter<T: IntoIterator<Item=Time>>(iter: T) -> Self {
		let mut list = iter.into_iter().collect::<Vec<_>>();
		list.sort_unstable();
		list.dedup();
		Self(list.into_iter())
	}
}

impl BitAnd for Times {
	type Output = Self;
	
	/// Intersection of times.
	fn bitand(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, true)
			.filter_map(|t| match t {
				SomeOrd::Equal(t) => Some(t),
				_ => None
			})
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl BitOr for Times {
	type Output = Self;
	
	/// Union of times.
	fn bitor(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, false)
			.map(SomeOrd::into_inner)
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl BitXor for Times {
	type Output = Self;
	
	/// Symmetric difference of times.
	fn bitxor(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, false)
			.filter_map(|t| match t {
				SomeOrd::Less(t) |
				SomeOrd::Greater(t) => Some(t),
				_ => None
			})
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl Not for Times {
	type Output = TimeRanges;
	
	/// Inverse of times.
	fn not(self) -> Self::Output {
		!self.map(|t| (t, t)).collect::<TimeRanges>()
	}
}

impl BitAnd<TimeRanges> for Times {
	type Output = Self;
	
	/// Intersection of times & ranges.
	fn bitand(self, rhs: TimeRanges) -> Self::Output {
		(self.map(|t| (t, t)).collect::<TimeRanges>() & rhs)
			.map(|(a, b)| {
				debug_assert_eq!(a, b);
				a
			})
			.collect()
	}
}

/// Iterator of [`Time`] ranges.
#[derive(Clone, Default)]
#[must_use]
pub struct TimeRanges(std::vec::IntoIter<(Time, Time)>);

impl Iterator for TimeRanges {
	type Item = (Time, Time);
	fn next(&mut self) -> Option<Self::Item> {
		self.0.next()
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.0.size_hint()
	}
	fn count(self) -> usize {
		self.0.count()
	}
}

impl FromIterator<(Time, Time)> for TimeRanges {
	fn from_iter<T: IntoIterator<Item=(Time, Time)>>(iter: T) -> Self {
		let mut prev = None;
		Self(iter.into_iter()
			.inspect(|&(a, b)| {
				assert!(a <= b);
				if let Some(prev) = std::mem::replace(&mut prev, Some(b)) {
					assert!(prev < a);
				}
			})
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl BitAnd for TimeRanges {
	type Output = Self;
	
	/// Intersection of ranges.
	fn bitand(self, rhs: Self) -> Self::Output {
		let mut prev = None;
		Self(min_iter(self, rhs, false)
			.filter_map(|curr| if let Some(prev_range) = prev {
				match (prev_range, curr) {
					(SomeOrd::Less((_, prev_b)), SomeOrd::Greater((a, b))) |
					(SomeOrd::Greater((_, prev_b)), SomeOrd::Less((a, b))) => {
						let r = if a <= prev_b {
							Some((a, prev_b.min(b)))
						} else {
							None
						};
						if b >= prev_b {
							prev = Some(curr);
						}
						r
					},
					(_, SomeOrd::Equal(curr)) => {
						prev = None;
						Some(curr)
					},
					_ => {
						prev = Some(curr);
						None
					}
				}
			} else {
				prev = Some(curr);
				None
			})
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl BitOr for TimeRanges {
	type Output = Self;
	
	/// Union of ranges.
	fn bitor(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, false)
			.fold(Vec::new(), |mut ranges, range| {
				let (a, b) = range.into_inner();
				if let Some((_, prev_b)) = ranges.last_mut() {
					if *prev_b < a {
						ranges.push((a, b));
					} else {
						*prev_b = b;
					}
				} else {
					ranges.push((a, b));
				}
				ranges
			})
			.into_iter())
	}
}

impl BitXor for TimeRanges {
	type Output = Self;
	
	/// Intersection of ranges.
	fn bitxor(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, false)
			.filter_map(|t| match t {
				SomeOrd::Less(t) |
				SomeOrd::Greater(t) => Some(t),
				_ => None
			})
			.fold(Vec::new(), |mut ranges, (a, b)| {
				match ranges.last_mut() {
					Some((prev_a, prev_b))
					if a <= *prev_b => {
						if *prev_a != a {
							let prev = std::mem::replace(prev_b, a - NANOSEC);
							if prev != b {
								ranges.push((prev + NANOSEC, b));
							}
						} else if *prev_b != b {
							*prev_a = *prev_b + NANOSEC;
							*prev_b = b;
						}
					},
					_ => ranges.push((a, b))
				}
				ranges
			})
			.into_iter())
	}
}

impl Not for TimeRanges {
	type Output = Self;
	
	/// Inverse of ranges.
	fn not(self) -> Self::Output {
		let mut ranges = Vec::new();
		let mut prev_time = None;
		for (a, b) in self {
			if prev_time.unwrap_or_default() != a {
				let t = prev_time
					.map(|t| t + Time::new(0, 1))
					.unwrap_or_default();
				
				ranges.push((t, a - Time::new(0, 1)));
			}
			prev_time = Some(b);
		}
		if prev_time.unwrap_or_default() != Time::MAX {
			let t = prev_time
				.map(|t| t + Time::new(0, 1))
				.unwrap_or_default();
			
			ranges.push((t, Time::MAX));
		}
		Self(ranges.into_iter())
	}
}

/// ...
#[derive(Copy, Clone, Debug)]
enum SomeOrd<T> {
	Less(T),
	Greater(T),
	Equal(T),
}

impl<T> SomeOrd<T> {
	pub fn into_inner(self) -> T {
		match self {
			SomeOrd::Less(t)    |
			SomeOrd::Greater(t) |
			SomeOrd::Equal(t)   => t
		}
	}
}

/// An iterator of two iterators. Compares the next item of the first iterator
/// to the second iterator. If equal, both iterators are advanced at once.
/// Otherwise, the lesser iterator is advanced and returned.
struct OrdIter<A: Iterator, B: Iterator> {
	a_iter: A,
	a_next: Option<A::Item>,
	b_iter: B,
	b_next: Option<B::Item>,
	break_none: bool,
}

impl<A, B> Iterator for OrdIter<A, B>
where
	A: Iterator,
	B: Iterator<Item = A::Item>,
	A::Item: Ord
{
	type Item = SomeOrd<A::Item>;
	fn next(&mut self) -> Option<Self::Item> {
		match (&mut self.a_next, &mut self.b_next) {
			(a @ Some(_), b @ Some(_)) => match a.cmp(&b) {
				Ordering::Less => {
					let a = std::mem::replace(a, self.a_iter.next());
					Some(SomeOrd::Less(a?))
				},
				Ordering::Greater => {
					let b = std::mem::replace(b, self.b_iter.next());
					Some(SomeOrd::Greater(b?))
				},
				Ordering::Equal => {
					let a = std::mem::replace(a, self.a_iter.next());
					*b = self.b_iter.next();
					Some(SomeOrd::Equal(a?))
				},
			},
			
			_ if self.break_none => None,
			
			(a @ Some(_), _) => {
				let a = std::mem::replace(a, self.a_iter.next());
				Some(SomeOrd::Less(a?))
			},
			(_, b @ Some(_)) => {
				let b = std::mem::replace(b, self.b_iter.next());
				Some(SomeOrd::Greater(b?))
			},
			
			_ => None
		}
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		let (a_min, a_max) = self.a_iter.size_hint();
		let (b_min, b_max) = self.b_iter.size_hint();
		let peek_num = (self.a_next.is_some() as usize) + (self.b_next.is_some() as usize);
		(
			a_min + b_min + peek_num,
			if let (Some(a_max), Some(b_max)) = (a_max, b_max) {
				Some(a_max + b_max + peek_num)
			} else {
				None
			}
		)
	}
}

fn min_iter<A, B>(mut a_iter: A, mut b_iter: B, break_none: bool) -> OrdIter<A, B>
where
	A: Iterator,
	B: Iterator,
	OrdIter<A, B>: Iterator,
{
	let a_next = a_iter.next();
	let b_next = b_iter.next();
	OrdIter {
		a_iter,
		a_next,
		b_iter,
		b_next,
		break_none,
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	
	#[test]
	fn time_logic() {
		let t = SEC;
		let a = Times::from_iter([1*t, 2*t, 3*t]);
		let b = Times::from_iter([2*t, 5*t]);
		assert_eq!(Vec::from_iter(a.clone() & b.clone()), [2*t]);
		assert_eq!(Vec::from_iter(a.clone() | b.clone()), [1*t, 2*t, 3*t, 5*t]);
		assert_eq!(Vec::from_iter(a.clone() ^ b.clone()), [1*t, 3*t, 5*t]);
		assert_eq!(Vec::from_iter(!b), [
			(Time::ZERO, 2*t - NANOSEC),
			(2*t + NANOSEC, 5*t - NANOSEC),
			(5*t + NANOSEC, Time::MAX)
		]);
		assert_eq!(
			Vec::from_iter(a.clone() & TimeRanges::from_iter([(Time::ZERO, Time::MAX)])),
			Vec::from_iter(a.clone())
		);
	}
	
	#[test]
	fn range_logic() {
		let t = NANOSEC;
		let a = TimeRanges::from_iter(
			[(1*t, 2*t), (3*t, 10*t), (20*t, 40*t)]);
		let b = TimeRanges::from_iter(
			[(2*t, 5*t), (20*t, 40*t), (50*t, Time::MAX)]);
		assert_eq!(Vec::from_iter(a.clone() & b.clone()), [
			(2*t, 2*t),
			(3*t, 5*t),
			(20*t, 40*t)
		]);
		assert_eq!(Vec::from_iter(a.clone() | b.clone()), [
			(1*t, 10*t),
			(20*t, 40*t),
			(50*t, Time::MAX)
		]);
		assert_eq!(Vec::from_iter(a.clone() ^ b.clone()), [
			(1*t, 2*t - NANOSEC),
			(5*t + NANOSEC, 10*t),
			(50*t, Time::MAX)
		]);
		assert_eq!(Vec::from_iter(!b), [
			(Time::ZERO, 2*t - NANOSEC),
			(5*t + NANOSEC, 20*t - NANOSEC),
			(40*t + NANOSEC, 50*t - NANOSEC)
		]);
	}
}