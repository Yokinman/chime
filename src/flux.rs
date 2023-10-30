//! Utilities for describing how a type changes over time.

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, Deref, DerefMut, Sub};

use time::{Time, TimeUnit};
use crate::linear::*;
use crate::poly::*;

mod impls;
mod kind;
pub use self::impls::*;
pub use self::kind::*;

/// A value that can change over time.
pub trait FluxValue: Sized {
	type Moment: FluxMoment<Self>;
	
	/// The kind of change over time.
	type Kind: FluxKind;
	
	/// The output accumulator of [`FluxValue::change`].
	type OutAccum<'a>: FluxAccum<'a, Self::Kind>
	where <Self::Kind as FluxKind>::Linear: 'a;
	
	/// Initial value.
	fn value(&self) -> <Self::Kind as FluxKind>::Linear;
	
	/// ...
	fn time(&self) -> Time;
	
	/// Accumulates change over time.
	fn change<'a>(&self, accum: <Self::Kind as FluxKind>::Accum<'a>) -> Self::OutAccum<'a>;
	
	/// The moment of this value at the given time.
	fn at(&self, time: Time) -> Self::Moment;
	
	/// Sets the moment of this value at the given time (affects all moments).
	fn set_at(&mut self, time: Time, moment: Self::Moment) {
		*self = moment.to_value(time);
	}
	
	/// Returns a mutable moment in time.
	/// 
	/// ```text
	/// let mut moment = self.at_mut(time);
	/// // modifications
	/// ```
	/// 
	/// equivalent to:
	/// 
	/// ```text
	/// let mut moment = self.at(time);
	/// // modifications
	/// self.set_moment(time, moment);
	/// ```
	fn at_mut(&mut self, time: Time) -> FluxMut<Self> {
		let moment = Some(self.at(time));
		FluxMut {
			inner: self,
			time,
			moment,
		}
	}
	
	/// The evaluation of this value at the given time.
	fn value_at(&self, time: Time) -> <Self::Kind as FluxKind>::Linear {
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
	
	/// A polynomial description of this flux value relative to `self.time()`.
	fn poly(&self) -> Poly<Self::Kind> {
		let mut poly = Poly::default();
		poly.0 = self.value();
		let poly_accum = FluxAccumKind::Poly { poly: &mut poly, depth: 0 };
		self.change(<Self::Kind as FluxKind>::Accum::from_kind(poly_accum));
		poly
	}
	
	/// Ranges of time when this value will be above/below/equal to the given value.
	fn when<O: FluxValue>(&self, cmp_order: Ordering, other: &O) -> TimeRanges
	where
		When<Self::Kind, O::Kind>: IntoIterator<IntoIter=TimeRanges>
	{
		let a_time = self.time();
		let b_time = other.time();
		let (a_poly, b_poly) = match a_time.cmp(&b_time) {
			Ordering::Less => (
				self.at(b_time).to_value(b_time).poly(), // !!! This may be lossy
				other.poly()
			),
			Ordering::Greater => (
				self.poly(),
				other.at(a_time).to_value(a_time).poly()
			),
			Ordering::Equal => (self.poly(), other.poly()),
		};
		When {
			a_poly,
			b_poly,
			cmp_order,
			time: a_time.max(b_time),
		}.into_iter()
	}
	
	/// Times when this value will be equal to the given value.
	fn when_eq<O: FluxValue>(&self, other: &O) -> Times
	where
		WhenEq<Self::Kind, O::Kind>: IntoIterator<IntoIter=Times>
	{
		let a_time = self.time();
		let b_time = other.time();
		let (a_poly, b_poly) = match a_time.cmp(&b_time) {
			Ordering::Less => (
				self.at(b_time).to_value(b_time).poly(),
				other.poly()
			),
			Ordering::Greater => (
				self.poly(),
				other.at(a_time).to_value(a_time).poly()
			),
			Ordering::Equal => (self.poly(), other.poly()),
		};
		WhenEq {
			a_poly,
			b_poly,
			time: a_time.max(b_time),
		}.into_iter()
	}
}

/// Moment-in-time interface for [`FluxValue::at`] / [`FluxValue::at_mut`].
pub trait FluxMoment<V: FluxValue<Moment=Self>> {
	/// Constructs the entirety of a [`FluxValue`] from a single moment.
	fn to_value(self, time: Time) -> V;
}

/// Mutable moment-in-time interface for [`FluxValue::at_mut`].
pub struct FluxMut<'v, V: FluxValue> {
	inner: &'v mut V,
	time: Time,
	moment: Option<V::Moment>,
}

impl<V: FluxValue> Drop for FluxMut<'_, V> {
	fn drop(&mut self) {
		if let Some(moment) = std::mem::take(&mut self.moment) {
			self.inner.set_at(self.time, moment);
		}
	}
}

impl<V: FluxValue> Deref for FluxMut<'_, V> {
	type Target = V::Moment;
	fn deref(&self) -> &Self::Target {
		if let Some(moment) = self.moment.as_ref() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<V: FluxValue> DerefMut for FluxMut<'_, V> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		if let Some(moment) = self.moment.as_mut() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<V: FluxValue> Debug for FluxMut<'_, V>
where
	V::Moment: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<V::Moment as Debug>::fmt(self, f)
	}
}

impl<V: FluxValue> Display for FluxMut<'_, V>
where
	V::Moment: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<V::Moment as Display>::fmt(self, f)
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

/// [`FluxValue::when`] predictive comparison.
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

/// [`FluxValue::when_eq`] predictive equality comparison.
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

/// Convenience for defining a generic self-isomorphic type.
pub trait Iso<A: LinearIso<B>, B: InvLinearIso<A> = A> {
	type Value;
}

impl<A: LinearIso<B>, B: LinearValue + InvLinearIso<A>> Iso<A, B> for Flux<A> {
	type Value = Self;
}

impl<A: LinearIso<B>, B: LinearValue + InvLinearIso<A>> Iso<A, B> for B {
	type Value = Self;
}

/// Convenience for encapsulating the unmapped values of a [`FluxValue`] type
/// with their time.
#[derive(Copy, Clone, Debug, Default)]
pub struct Flux<T: LinearValue>(Time, T);

impl<T: LinearValue> Flux<T> {
	pub fn new(time: Time, value: T) -> Self {
		Self(time, value)
	}
	
	pub fn time(&self) -> Time {
		self.0
	}
}

impl<T: LinearValue> Deref for Flux<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.1
	}
}

impl<T: LinearValue> DerefMut for Flux<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.1
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use TimeUnit::*;
	
	macro_rules! make {
		($name:ident<$t:ident> { $($sub_name:ident: $sub_type:ty),* }) => {
			#[derive(Debug, Default, Clone)]
			struct $name<$t: Iso<f64>> {
				value: $t::Value,
				$($sub_name: $sub_type),*
			}
			
			impl Deref for $name<f64> {
				type Target = f64;
				fn deref(&self) -> &Self::Target {
					&self.value
				}
			}
			
			impl DerefMut for $name<f64> {
				fn deref_mut(&mut self) -> &mut Self::Target {
					&mut self.value
				}
			}
		}
	}
	
	make!(Pos   <I> { spd: Spd<I>, misc: Vec<Spd<I>> });
	make!(Spd   <I> { fric: Fric<I>, accel: Accel<I> });
	make!(Fric  <I> { });
	make!(Accel <I> { jerk: Jerk<I> });
	make!(Jerk  <I> { snap: Snap<I> });
	make!(Snap  <I> { });
	
	impl FluxValue for Pos<Flux<f64>> {
		type Moment = Pos<f64>;
		type Kind = Sum<f64, 4>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.add(&self.spd, TimeUnit::Secs)
				.add(&self.misc, TimeUnit::Secs)
		}
		fn at(&self, time: Time) -> Self::Moment {
			Pos {
				value: self.value_at(time).inv_map(),
				spd: self.spd.at(time),
				misc: self.misc.at(time),
			}
		}
	}
	impl FluxMoment<Pos<Flux<f64>>> for Pos<f64> {
		fn to_value(self, time: Time) -> Pos<Flux<f64>> {
			Pos {
				value: Flux::new(time, self.value.inv_map()),
				spd: self.spd.to_value(time),
				misc: self.misc.to_value(time),
			}
		}
	}
	
	impl FluxValue for Spd<Flux<f64>> {
		type Moment = Spd<f64>;
		type Kind = Sum<f64, 3>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
				.sub(&self.fric, TimeUnit::Secs)
				.add(&self.accel, TimeUnit::Secs)
		}
		fn at(&self, time: Time) -> Self::Moment {
			Spd {
				value: self.value_at(time).inv_map(),
				fric: self.fric.at(time),
				accel: self.accel.at(time),
			}
		}
	}
	impl FluxMoment<Spd<Flux<f64>>> for Spd<f64> {
		fn to_value(self, time: Time) -> Spd<Flux<f64>> {
			Spd {
				value: Flux::new(time, self.value.inv_map()),
				fric: self.fric.to_value(time),
				accel: self.accel.to_value(time),
			}
		}
	}
	
	impl FluxValue for Fric<Flux<f64>> {
		type Moment = Fric<f64>;
		type Kind = Sum<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn at(&self, time: Time) -> Self::Moment {
			Fric {
				value: self.value_at(time).inv_map(),
			}
		}
	}
	impl FluxMoment<Fric<Flux<f64>>> for Fric<f64> {
		fn to_value(self, time: Time) -> Fric<Flux<f64>> {
			Fric {
				value: Flux::new(time, self.value.inv_map()),
			}
		}
	}
	
	impl FluxValue for Accel<Flux<f64>> {
		type Moment = Accel<f64>;
		type Kind = Sum<f64, 2>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.jerk, TimeUnit::Secs)
		}
		fn at(&self, time: Time) -> Self::Moment {
			Accel {
				value: self.value_at(time).inv_map(),
				jerk: self.jerk.at(time)
			}
		}
	}
	impl FluxMoment<Accel<Flux<f64>>> for Accel<f64> {
		fn to_value(self, time: Time) -> Accel<Flux<f64>> {
			Accel {
				value: Flux::new(time, self.value.inv_map()),
				jerk: self.jerk.to_value(time),
			}
		}
	}
	
	impl FluxValue for Jerk<Flux<f64>> {
		type Moment = Jerk<f64>;
		type Kind = Sum<f64, 1>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes.add(&self.snap, TimeUnit::Secs)
		}
		fn at(&self, time: Time) -> Self::Moment {
			Jerk {
				value: self.value_at(time).inv_map(),
				snap: self.snap.at(time),
			}
		}
	}
	impl FluxMoment<Jerk<Flux<f64>>> for Jerk<f64> {
		fn to_value(self, time: Time) -> Jerk<Flux<f64>> {
			Jerk {
				value: Flux::new(time, self.value.inv_map()),
				snap: self.snap.to_value(time),
			}
		}
	}
	
	impl FluxValue for Snap<Flux<f64>> {
		type Moment = Snap<f64>;
		type Kind = Sum<f64, 0>;
		type OutAccum<'a> = SumAccum<'a, Self::Kind>;
		fn value(&self) -> <Self::Kind as FluxKind>::Linear {
			*self.value
		}
		fn time(&self) -> Time {
			self.value.time()
		}
		fn change<'a>(&self, changes: Changes<'a, Self>) -> Self::OutAccum<'a> {
			changes
		}
		fn at(&self, time: Time) -> Self::Moment {
			Snap {
				value: self.value_at(time).inv_map(),
			}
		}
	}
	impl FluxMoment<Snap<Flux<f64>>> for Snap<f64> {
		fn to_value(self, time: Time) -> Snap<Flux<f64>> {
			Snap {
				value: Flux::new(time, self.value.inv_map()),
			}
		}
	}
	
	fn position() -> Pos<Flux<f64>> {
		let pos = Pos {
			value: 32.0, 
			spd: Spd {
				value: 0.0,
				fric: Fric { value: 3.5 },
				accel: Accel {
					value: 0.3,
					jerk: Jerk {
						value: 0.4,
						snap: Snap { value: -0.01 },
					},
				},
			},
			misc: Vec::new(),
		};
		pos.to_value(Time::default())
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		
		 // Values:
		assert_eq!(pos.at(0*Secs).round(), 32.);
		assert_eq!(pos.at(10*Secs).round(), -63.);
		assert_eq!(pos.at(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
		
		 // Update:
		assert_eq!(pos.at_mut(20*Secs).round(), -113.);
		assert_eq!(pos.at(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
		assert_eq!(pos.at_mut(100*Secs).round(), 8339.);
		assert_eq!(pos.at(200*Secs).round(), -209779.);
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
			pos.at_mut(20*Secs);
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
		pos.at_mut(0*Secs);
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
		
		// pos.at_mut(20*Secs);
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
	fn when_at_mut() {
		let mut a_pos = position();
		let mut b_pos = position();
		
		 // Check Before:
		let vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		assert_eq!(vec, []);
		assert_eq!(vec_eq, [0*Nanosecs]);
		a_pos.at_mut(20*Secs);
		
		 // Apply Changes:
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.at_mut(0*Secs).misc.push(Spd {
			value: 12.25,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.at_mut(10*Secs).value -= 100.0;
		
		 // Check After:
		let mut vec: Vec<(Time, Time)> = a_pos.when(Ordering::Greater, &b_pos).collect();
		let mut vec_eq: Vec<Time> = a_pos.when_eq(&b_pos).collect();
		for (a, b) in &mut vec {
			*a = Time::from_secs(a.as_secs_f64().round() as u64);
			*b = Time::from_secs(b.as_secs_f64().round() as u64);
		}
		for t in &mut vec_eq {
			*t = Time::from_secs(t.as_secs_f64().round() as u64);
		}
		assert_eq!(vec, [(0*Secs, 8*Secs), (50*Secs, vec[1].1)]);
		assert_eq!(vec_eq[0..2], [8*Secs, 50*Secs]);
		assert!(vec[1].1 > 2_u64.pow(53)*Secs);
		// https://www.desmos.com/calculator
	}
}