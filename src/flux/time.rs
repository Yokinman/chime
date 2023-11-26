//! ...

use impl_op::impl_op;

/// ...
pub type Time = std::time::Duration;

/// A nanosecond-based interface for creating & parsing instances of [Time].
#[derive(Debug, Copy, Clone)]
pub enum TimeUnit {
	/** Nanoseconds  (ns)          */ Nanosecs, 
	/** Microseconds (μs), 1,000ns */ Microsecs,
	/** Milliseconds (ms), 1,000μs */ Millisecs,
	/** Seconds      (s),  1,000ms */ Secs,
	/** Minutes      (m),  60s     */ Mins,
	/** Hours        (h),  60m     */ Hours,
}

impl TimeUnit {
	const fn length(self) -> u64 {
		match self {
			Self::Nanosecs  => 1,
			Self::Microsecs => 1000 * Self::Nanosecs.length(),
			Self::Millisecs => 1000 * Self::Microsecs.length(),
			Self::Secs      => 1000 * Self::Millisecs.length(),
			Self::Mins      => 60 * Self::Secs.length(),
			Self::Hours     => 60 * Self::Mins.length(),
		}
	}
	
	pub const fn time(self, value: u64) -> Time {
		Time::from_nanos(value.saturating_mul(self.length()))
	}
}

 // Conversion:
impl_op!{ a >> b -> u64 {
	(Time,     TimeUnit) => (a.as_nanos() / (b.length() as u128)) as u64,
	(TimeUnit, TimeUnit) => a.length() / b.length(),
}}

 // Scalar Multiplication:
impl_op!{ a * b -> Time {
	(u64, TimeUnit) => b.time(a),
}}

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