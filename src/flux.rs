//! Utilities for describing how a type changes over time.

use std::cmp::Ordering;
use std::ops::{Add, Sub};

use time::{Time, TimeUnit};
use crate::linear::*;
use crate::poly::*;

mod impls;
mod kind;
pub use self::impls::*;
pub use self::kind::*;

/// A value that can change over time.
pub trait FluxValue: Sized {
	/// The produced value.
	type Value;
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// The output accumulator of [`FluxValue::change`].
	type OutAccum<'a>: FluxAccum<'a, Self::Kind>
	where <Self::Kind as FluxKind>::Linear: 'a;
	
	/// Initial value.
	fn value(&self) -> <Self::Kind as FluxKind>::Linear;
	fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear);
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> Self::OutAccum<'a>;
	
	/// Apply change over time.
	fn time(&self) -> Time;
	fn set_time(&mut self, time: Time);
	
	fn calc_value(&self, time: Time) -> <Self::Kind as FluxKind>::Linear {
		let mut value = self.value();
		let value_accum = FluxAccumKind::Sum {
			sum: &mut value,
			depth: 0,
			time,
			offset: self.time(),
		};
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(value_accum));
		value
	}
	
	fn at(&self, time: Time) -> Self::Value
	where
		<Self::Kind as FluxKind>::Linear: LinearIso<Self::Value>,
	{
		self.calc_value(time).map()
	}
	
	fn poly(&self) -> Poly<Self::Kind> {
		let mut poly = Poly::default();
		poly.0 = self.value();
		let poly_accum = FluxAccumKind::Poly { poly: &mut poly, depth: 0 };
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(poly_accum));
		poly
	}
	
	/// ...
	fn when<B: FluxValue>(&self, cmp_order: Ordering, other: &B) -> TimeRanges
	where
		Self: Clone,
		B: Clone,
		When<Self::Kind, B::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let a_time = self.time();
		let b_time = other.time();
		let (a_poly, b_poly) = match a_time.cmp(&b_time) {
			Ordering::Less => {
				let mut a = self.clone();
				a.set_time(b_time);
				(a.poly(), other.poly())
			},
			Ordering::Greater => {
				let mut b = other.clone();
				b.set_time(a_time);
				(self.poly(), b.poly())
			},
			Ordering::Equal => (self.poly(), other.poly()),
		};
		When {
			a_poly,
			b_poly,
			cmp_order,
			time: a_time.max(b_time),
		}.into_iter()
	}
	
	/// ...
	fn when_eq<B: FluxValue>(&self, other: &B) -> Times
	where
		Self: Clone,
		B: Clone,
		WhenEq<Self::Kind, B::Kind>: IntoIterator<IntoIter=Times>
	{
		let a_time = self.time();
		let b_time = other.time();
		let (a_poly, b_poly) = match a_time.cmp(&b_time) {
			Ordering::Less => {
				let mut a = self.clone();
				a.set_time(b_time);
				(a.poly(), other.poly())
			},
			Ordering::Greater => {
				let mut b = other.clone();
				b.set_time(a_time);
				(self.poly(), b.poly())
			},
			Ordering::Equal => (self.poly(), other.poly()),
		};
		WhenEq {
			a_poly,
			b_poly,
			time: a_time.max(b_time),
		}.into_iter()
	}
}

/// Iterator of [`Time`] ranges.
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

/// Iterator of [`Time`] values.
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

#[derive(Copy, Clone, Debug)]
pub struct When<A: FluxKind, B: FluxKind> {
	a_poly: Poly<A>,
	b_poly: Poly<B>,
	cmp_order: Ordering,
	time: Time,
}

impl<A, B> IntoIterator for When<A, B>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Linear=A::Linear>,
	A::Linear: PartialOrd,
	<A as Add<B>>::Output: FluxKind<Linear=A::Linear> + Roots + PartialOrd,
	Poly<A>: Sub<Poly<B>, Output=Poly<<A as Add<B>>::Output>>,
{
	type Item = (Time, Time);
	type IntoIter = TimeRanges;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = self.a_poly - self.b_poly;
		
		 // Find Initial Order:
		let mut degree = <<A as Add<B>>::Output as FluxKind>::DEGREE;
		let initial_order = loop {
			if degree == 0 {
				break poly.constant().partial_cmp(&A::Linear::zero());
			} else {
				let coeff = poly.coeff(degree - 1).unwrap();
				let order = coeff.partial_cmp(&FluxKind::zero());
				if order != Some(Ordering::Equal) {
					break if degree % 2 == 0 {
						order
					} else {
						order.map(|o| o.reverse())
					}
				}
				degree -= 1;
			}
		};
		
		 // Convert Roots to Ranges:
		let range_list = match poly.real_roots() {
			Ok(roots) => {
				let mut list = Vec::with_capacity(
					if self.cmp_order == Ordering::Equal {
						roots.len()
					} else {
						1 + (roots.len() / 2)
					}
				);
				let mut prev_point = if initial_order == Some(self.cmp_order) {
					Some(f64::NEG_INFINITY)
				} else {
					None
				};
				for &point in roots.into_iter() {
					if let Some(prev) = prev_point {
						if point != prev {
							list.push((prev, point));
						}
						prev_point = None;
					} else if self.cmp_order == Ordering::Equal {
						list.push((point, point));
					} else {
						prev_point = Some(point);
					}
				}
				if let Some(point) = prev_point {
					list.push((point, f64::INFINITY));
				}
				list.into_iter()
			},
			Err(_) => Default::default(),
		};
		
		 // Convert Ranges to Times:
		let time = self.time.as_secs_f64();
		let max_time = Time::MAX.as_secs_f64();
		let list: Vec<(Time, Time)> = range_list
			.filter_map(|(a, b)| {
				debug_assert!(a <= b);
				let a = a + time;
				let b = b + time;
				if b < 0.0 || a >= max_time {
					None
				} else {
					Some((
						Time::try_from_secs_f64(a).unwrap_or(Time::ZERO),
						Time::try_from_secs_f64(b).unwrap_or(Time::MAX)
					))
				}
			})
			.collect();
		
		TimeRanges(list.into_iter())
	}
}

/// [`FluxValue`] predictive equality comparison.
#[derive(Copy, Clone, Debug)]
pub struct WhenEq<A: FluxKind, B: FluxKind> {
	a_poly: Poly<A>,
	b_poly: Poly<B>,
	time: Time,
}

impl<A, B> IntoIterator for WhenEq<A, B>
where
	A: FluxKind + Add<B>,
	B: FluxKind<Linear=A::Linear>,
	A::Linear: PartialEq,
	<A as Add<B>>::Output: FluxKind<Linear=A::Linear> + Roots + PartialEq,
	Poly<A>: Sub<Poly<B>, Output=Poly<<A as Add<B>>::Output>>,
{
	type Item = Time;
	type IntoIter = Times;
	
	fn into_iter(self) -> Self::IntoIter {
		let poly = self.a_poly - self.b_poly;
		let mut real_roots = poly.real_roots().unwrap_or(Box::default())
			.into_vec();
		
		 // Constant Equality:
		if
			real_roots.is_empty()
			&& poly.constant().is_zero()
			&& poly.coeff_iter().all(FluxKind::is_zero)
		{
			real_roots.push(0.0);
		}
		
		 // Convert Roots to Times:
		let time = self.time.as_secs_f64();
		let max_time = Time::MAX.as_secs_f64();
		let vec: Vec<Time> = real_roots.into_iter()
			.filter_map(|t| {
				let t = t + time;
				if t < 0.0 || t >= max_time {
					None
				} else {
					Some(Time::from_secs_f64(t))
				}
			})
			.collect();
		
		Times(vec.into_iter())
	}
}

#[allow(type_alias_bounds)]
pub type Changes<'a, T: FluxValue> = <T::Kind as FluxKind>::Accum<'a>;

/// Change accumulator.
/// 
/// Converts a discrete pattern of change into a desired form.
pub trait FluxAccum<'a, K: FluxKind> {
	fn from_kind(kind: FluxAccumKind<'a, K>) -> Self;
}

/// General accumulator arguments.
#[non_exhaustive]
pub enum FluxAccumKind<'a, K: FluxKind> {
	Sum {
		sum: &'a mut K::Linear,
		depth: usize,
		time: Time,
		offset: Time,
	},
	Poly {
		poly: &'a mut Poly<K>,
		depth: usize,
	},
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	
	#[derive(Debug, Default, Clone)] struct Pos   { time: Time, value: f64, spd: Spd, misc: Vec<Spd> }
	#[derive(Debug, Default, Clone)] struct Spd   { time: Time, value: f64, fric: Fric, accel: Accel }
	#[derive(Debug, Default, Clone)] struct Fric  { time: Time, value: f64 }
	#[derive(Debug, Default, Clone)] struct Accel { time: Time, value: f64, jerk: Jerk }
	#[derive(Debug, Default, Clone)] struct Jerk  { time: Time, value: f64, snap: Snap }
	#[derive(Debug, Default, Clone)] struct Snap  { time: Time, value: f64 }
	
	impl FluxValue for Pos {
		type Value = i64;
		type Kind = Sum<f64, 4>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.add(&self.spd, TimeUnit::Secs)
				.add(&self.misc, TimeUnit::Secs)
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
			self.spd.set_time(time);
			self.misc.set_time(time);
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	impl FluxValue for Spd {
		type Value = i64;
		type Kind = Sum<f64, 3>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.sub(&self.fric, TimeUnit::Secs)
				.add(&self.accel, TimeUnit::Secs)
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
			self.fric.set_time(time);
			self.accel.set_time(time);
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	impl FluxValue for Fric {
		type Value = i64;
		type Kind = Sum<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	impl FluxValue for Accel {
		type Value = i64;
		type Kind = Sum<f64, 2>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.jerk, TimeUnit::Secs)
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
			self.jerk.set_time(time);
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	impl FluxValue for Jerk {
		type Value = i64;
		type Kind = Sum<f64, 1>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.snap, TimeUnit::Secs)
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
			self.snap.set_time(time);
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	impl FluxValue for Snap {
		type Value = i64;
		type Kind = Sum<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			self.value
		}
		fn set_value(&mut self, value: <Self::Kind as FluxKind>::Linear) {
			self.value = value;
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn set_time(&mut self, time: Time) {
			self.value = self.calc_value(time);
			self.time = time;
		}
		fn time(&self) -> Time {
			self.time
		}
	}
	
	fn position() -> Pos {
		Pos {
			time: Time::default(),
			value: 32.0, 
			spd: Spd {
				time: Time::default(),
				value: 0.0,
				fric: Fric {
					time: Time::default(),
					value: 3.5,
				},
				accel: Accel {
					time: Time::default(),
					value: 0.3,
					jerk: Jerk {
						time: Time::default(),
						value: 0.4,
						snap: Snap {
							time: Time::default(),
							value: -0.01
						},
					},
				},
			},
			misc: Vec::new(),
		}
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		
		 // Values:
		assert_eq!(pos.at(0*Secs), 32);
		assert_eq!(pos.at(10*Secs), -63);
		assert_eq!(pos.at(20*Secs), -113);
		assert_eq!(pos.at(100*Secs), 8339);
		assert_eq!(pos.at(200*Secs), -209779);
		
		 // Update:
		pos.set_time(20*Secs);
		assert_eq!(pos.at(20*Secs), -113);
		assert_eq!(pos.at(100*Secs), 8339);
		assert_eq!(pos.at(200*Secs), -209779);
		pos.set_time(55*Secs);
		assert_eq!(pos.at(100*Secs), 8339);
		assert_eq!(pos.at(200*Secs), -209779);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_eq!(
			pos.poly(),
			Poly(32.0, [
				Sum(-1.4691666666666667),
				Sum(-1.4045833333333335),
				Sum(0.06416666666666666),
				Sum(-0.00041666666666666664),
			])
		);
		for _ in 0..2 {
			pos.set_time(20*Secs);
			assert_eq!(
				pos.poly(),
				Poly(-112.55000000000007, [
					Sum(6.014166666666661),
					Sum(1.44541666666666666),
					Sum(0.030833333333333334),
					Sum(-0.00041666666666666664),
				])
			);
		}
		pos.set_time(0*Secs);
		assert_eq!(
			pos.poly(),
			Poly(32.00000000000006, [
				Sum(-1.4691666666666667),
				Sum(-1.4045833333333335),
				Sum(0.06416666666666666),
				Sum(-0.00041666666666666664),
			])
		);
	}
	
	#[test]
	fn when() {
		let pos = position();
		let acc = position().spd.accel;
		
		let pos_when: Vec<(Time, Time)> = pos.when(Ordering::Greater, &acc).collect();
		let pos_when_eq: Vec<Time> = pos.when_eq(&acc).collect();
		
		assert_eq!(pos_when, [
			(0*Nanosecs, 4560099744*Nanosecs),
			(26912076291*Nanosecs, 127394131312*Nanosecs)
		]);
		assert_eq!(pos_when_eq, [
			4560099744*Nanosecs,
			26912076291*Nanosecs,
			127394131312*Nanosecs
		]);
		
		// pos.set_time(20*Secs);
		// pos.borrow_mut().spd.value = -20.0;
		// pos.borrow_mut().spd.accel.jerk.value = 0.3;
		// 
		// assert_eq!(pos_when.into_iter().collect::<Vec<(Time, Time)>>(), [
		// 	(0*Nanosecs, 4560099744*Nanosecs),
		// 	(33544693273*Nanosecs, 157824014330*Nanosecs)
		// ]);
		// assert_eq!(pos.when(Ordering::Greater, &Snap { value: -116.0 }).into_iter().collect::<Vec<(Time, Time)>>(), [
		// 	(0*Nanosecs, 16114479737*Nanosecs),
		// 	(19315363058*Nanosecs, 20188850711*Nanosecs),
		// 	(29523688931*Nanosecs, 157716291466*Nanosecs)
		// ]);
	}
	
	#[test]
	fn when_borrow() {
		let mut a_pos = position();
		let mut b_pos = position();
		
		 // Check Before:
		let vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		assert_eq!(vec, []);
		assert_eq!(vec_eq, [0*Nanosecs]);
		
		 // Apply Changes:
		a_pos.set_time(20*Secs);
		b_pos.misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.misc.push(Spd {
			value: 12.25,
			fric: Fric {
				time: Time::default(),
				value: 0.5,
			},
			..Default::default()
		});
		
		 // Check After:
		let vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		assert_eq!(vec, [(0*Secs, 0*Secs), (58*Secs, Time::MAX)]);
		assert_eq!(vec_eq, [0*Secs, 58*Secs]);
	}
}