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
#[derive(Clone)]
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
		Self(list.into_iter())
	}
}

impl BitAnd for Times {
	//! Intersection of times.
	
	type Output = Self;
	
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
	//! Union of times.
	
	type Output = Self;
	
	fn bitor(self, rhs: Self) -> Self::Output {
		Self(min_iter(self, rhs, false)
			.map(|t| match t {
				SomeOrd::Less(t) |
				SomeOrd::Greater(t) |
				SomeOrd::Equal(t) => t,
			})
			.collect::<Vec<_>>()
			.into_iter())
	}
}

impl BitXor for Times {
	//! Symmetric difference of times.
	
	type Output = Self;
	
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
	fn not(self) -> Self::Output {
		let mut ranges = Vec::new();
		let mut prev_time = None;
		for time in self {
			if prev_time.unwrap_or_default() != time {
				let t = prev_time
					.map(|t| t + Time::new(0, 1))
					.unwrap_or_default();
				
				ranges.push((t, time - Time::new(0, 1)));
			}
			prev_time = Some(time);
		}
		if prev_time.unwrap_or_default() != Time::MAX {
			let t = prev_time
				.map(|t| t + Time::new(0, 1))
				.unwrap_or_default();
			
			ranges.push((t, Time::MAX));
		}
		TimeRanges(ranges.into_iter())
	}
}

/// Iterator of [`Time`] ranges.
#[derive(Clone)]
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
		let mut list = iter.into_iter().collect::<Vec<_>>();
		list.sort_unstable_by_key(|(a, b)| {
			assert!(a <= b);
			*a
		});
		assert!(list.windows(2).all(|w| w[0].1 <= w[1].0));
		Self(list.into_iter())
	}
}

/// ...
enum SomeOrd<T> {
	Less(T),
	Greater(T),
	Equal(T),
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
		let a = [1*SEC, 2*SEC, 3*SEC].into_iter().collect::<Times>();
		let b = [2*SEC, 5*SEC].into_iter().collect::<Times>();
		assert_eq!((a.clone() & b.clone()).collect::<Vec<_>>(), [2*SEC]);
		assert_eq!((a.clone() | b.clone()).collect::<Vec<_>>(), [1*SEC, 2*SEC, 3*SEC, 5*SEC]);
		assert_eq!((a.clone() ^ b.clone()).collect::<Vec<_>>(), [1*SEC, 3*SEC, 5*SEC]);
		assert_eq!((!b).collect::<Vec<_>>(), [
			(Time::ZERO, 2*SEC - NANOSEC),
			(2*SEC + NANOSEC, 5*SEC - NANOSEC),
			(5*SEC + NANOSEC, Time::MAX)
		]);
	}
}