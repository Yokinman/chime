//! Defining a kind of change over time.

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, Deref, DerefMut, Mul, Sub};

use crate::linear::{Scalar, Vector, LinearIso};
use crate::time::Time;
use crate::{Constant, Flux, Moment, MomentMut, ToMoment, ToMomentMut};
use crate::kind::{FluxIntegral, FluxKind, KindLinear, ops as kind_ops, Roots};
use crate::pred::{When, WhenDis, WhenDisEq, WhenEq};

/// A [`FluxKind`] paired with a basis time.
/// e.g. `Temporal<1 + 2x>` => `1 + 2(x-time)`.
#[derive(Copy, Clone, Debug, Default)]
#[cfg_attr(feature = "bevy", derive(
	bevy_ecs::component::Component,
	bevy_ecs::system::Resource,
))]
pub struct Temporal<T> {
	pub(crate) inner: T,
	pub(crate) time: Time,
}

/// ... [`<Temporal as IntoIterator>::IntoIter`]
pub struct TemporalIter<T> {
	iter: T,
	time: Time,
}

impl<K> PartialEq for Temporal<K>
where
	K: PartialEq
{
	fn eq(&self, other: &Self) -> bool {
		// ??? Deriving PartialEq, Eq could count `f(t) = 1 + 2t` and
 		// `g(t) = 3 + 2(t-basis_time)` as the same if `basis_time = 1`.
		self.inner == other.inner && self.time == other.time
	}
}
	
impl<T> Deref for Temporal<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.inner
	}
}

impl<T> DerefMut for Temporal<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.inner
	}
}

impl<T> Temporal<T> {
	pub fn new(inner: T, time: Time) -> Self {
		Self { inner, time }
	}
	
	pub fn map<U>(&self, f: impl Fn(&T) -> &U) -> Temporal<&'_ U> {
		Temporal {
			inner: f(&self.inner),
			time: self.time,
		}
	}
	
	fn secs(&self, time: Time) -> Scalar {
		Scalar::from(if time > self.time {
			(time - self.time).as_secs_f64()
		} else {
			-(self.time - time).as_secs_f64()
		})
	}
}

impl<T: Flux> Temporal<T> {
	/// An evaluation of this flux at some point in time.
	pub fn basis(&self) -> <T::Kind as FluxKind>::Basis {
		self.inner.basis()
	}
	
	/// The time of [`Flux::basis`].
	pub fn basis_time(&self) -> Time {
		self.time
	}
	
	/// Conversion into a standard representation.
	pub fn to_kind(&self) -> Temporal<T::Kind> {
		Temporal {
			inner: self.inner.to_kind(),
			time: self.time,
		}
	}
	
	/// A polynomial description of this flux at the given time.
	pub fn poly(&self, time: Time) -> Temporal<T::Kind> {
		let mut inner = self.inner.to_kind();
		let _ = inner.to_moment_mut(self.secs(time));
		Temporal { inner, time }
	}
	
	/// Ranges when this is above/below/equal to another flux.
	pub fn when<U>(&self, cmp: Ordering, other: &Temporal<U>)
		-> <Temporal<T::Kind> as When<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<T::Kind>: When<U::Kind>
	{
		let time = self.basis_time();
		When::when(self.poly(time), cmp, other.poly(time))
	}
	
	/// Times when this is equal to another flux.
	pub fn when_eq<U>(&self, other: &Temporal<U>)
		-> <Temporal<T::Kind> as WhenEq<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<T::Kind>: WhenEq<U::Kind>
	{
		let time = self.basis_time();
		WhenEq::when_eq(self.poly(time), other.poly(time))
	}
	
	/// Ranges when this is above/below/equal to a constant.
	pub fn when_constant<U>(&self, cmp: Ordering, other: U)
		-> <Temporal<T::Kind> as When<Constant<KindLinear<T::Kind>>>>::Pred
	where
		U: LinearIso<KindLinear<T::Kind>>,
		Temporal<T::Kind>: When<Constant<KindLinear<T::Kind>>>
	{
		self.when(cmp, &Temporal::new(Constant::from(U::into_linear(other)), Time::ZERO))
	}
	
	/// Times when this is equal to a constant.
	pub fn when_eq_constant<U>(&self, other: U)
		-> <Temporal<T::Kind> as WhenEq<Constant<KindLinear<T::Kind>>>>::Pred
	where
		U: LinearIso<KindLinear<T::Kind>>,
		Temporal<T::Kind>: WhenEq<Constant<KindLinear<T::Kind>>>
	{
		self.when_eq(&Temporal::new(Constant::from(U::into_linear(other)), Time::ZERO))
	}
}

impl<K: FluxKind> Temporal<K> {
	pub fn eval(&self, time: Time) -> K::Basis {
		self.inner.eval(self.secs(time))
	}
	
	pub fn initial_order(&self, time: Time) -> Option<Ordering>
	where
		KindLinear<K>: PartialOrd
	{
		self.inner.initial_order(self.secs(time))
	}
	
	pub fn deriv(mut self) -> Self {
		self.inner = self.inner.deriv();
		self
	}
	
	pub fn integ(self) -> Temporal<K::Integ>
	where
		K: FluxIntegral,
	{
		Temporal::new(self.inner.integ(), self.time)
	}
}

impl<T: ToMoment> Temporal<T> {
	/// See [`ToMoment::to_moment`].
	pub fn to_moment(&self, time: Time) -> T::Moment<'_> {
		self.inner.to_moment(self.secs(time))
	}
	
	pub fn moment(&self, time: Time) -> Moment<T> {
		Moment {
			moment: self.to_moment(time),
			borrow: std::marker::PhantomData,
		}
	}
}

impl<T: ToMomentMut> Temporal<T> {
	/// See [`ToMomentMut::to_moment_mut`].
	pub fn to_moment_mut(&mut self, time: Time) -> T::MomentMut<'_> {
		let secs = self.secs(time);
		self.time = time;
		self.inner.to_moment_mut(secs)
	}
	
	pub fn moment_mut(&mut self, time: Time) -> MomentMut<T> {
		MomentMut {
			moment: self.to_moment_mut(time),
			borrow: std::marker::PhantomData,
		}
	}
	
	pub fn at_time(mut self, time: Time) -> Self {
		if self.time != time {
			let _ = self.inner.to_moment_mut(self.secs(time));
			self.time = time;
		}
		self
	}
}

impl<K: FluxKind> Temporal<K> {
	pub fn sqr(self) -> Temporal<<K as kind_ops::Sqr>::Output>
	where
		K: kind_ops::Sqr
	{
		Temporal::new(self.inner.sqr(), self.time)
	}
	
	/// Ranges when the sign is greater than, less than, or equal to zero.
	pub(crate) fn when_sign<F>(self, order: Ordering, filter: F) -> crate::pred::PredFilter<crate::pred::Pred<K>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialEq,
		KindLinear<K>: PartialOrd,
	{
		let pred = crate::pred::Pred {
			poly: self,
			order
		};
		crate::pred::PredFilter { pred, filter }
	}
	
	/// Times when the value is equal to zero.
	pub(crate) fn when_zero<F>(self, filter: F) -> crate::pred::PredFilter<crate::pred::PredEq<K>, F>
	where
		F: crate::pred::TimeFilterMap,
		K: Roots + PartialEq,
		KindLinear<K>: PartialEq,
	{
		crate::pred::PredFilter {
			pred: crate::pred::PredEq { poly: self },
			filter,
		}
	}
}

impl<K: FluxKind> From<K> for Temporal<K> {
	fn from(value: K) -> Self {
		Self::new(value, Time::ZERO)
	}
}

impl<A: FluxKind, B: FluxKind> Add<Temporal<B>> for Temporal<A>
where
	A: kind_ops::Add<B>
{
	type Output = Temporal<<A as kind_ops::Add<B>>::Output>;
	fn add(self, rhs: Temporal<B>) -> Self::Output {
		Temporal {
			inner: self.inner.add(rhs.at_time(self.time).inner),
			time: self.time,
		}
	}
}

impl<A: FluxKind, B: FluxKind> Sub<Temporal<B>> for Temporal<A>
where
	A: kind_ops::Sub<B>
{
	type Output = Temporal<<A as kind_ops::Sub<B>>::Output>;
	fn sub(self, rhs: Temporal<B>) -> Self::Output {
		Temporal {
			inner: self.inner.sub(rhs.at_time(self.time).inner),
			time: self.time,
		}
	}
}

impl<K> Mul<Scalar> for Temporal<K>
where
	K: Mul<Scalar, Output=K>
{
	type Output = Self;
	fn mul(mut self, rhs: Scalar) -> Self::Output {
		self.inner = self.inner * rhs;
		self
	}
}

impl<K, const SIZE: usize> Vector<SIZE> for Temporal<K>
where
	K: Vector<SIZE>,
{
	type Output = Temporal<K::Output>;
	fn index(&self, index: usize) -> Self::Output {
		Temporal {
			inner: self.inner.index(index),
			time: self.time,
		}
	}
}

impl<K> IntoIterator for Temporal<K>
where
	K: IntoIterator<Item: FluxKind>,
{
	type Item = Temporal<K::Item>;
	type IntoIter = TemporalIter<K::IntoIter>;
	fn into_iter(self) -> Self::IntoIter {
		TemporalIter {
			iter: self.inner.into_iter(),
			time: self.time,
		}
	}
}

impl<T> Iterator for TemporalIter<T>
where
	T: Iterator<Item: FluxKind>,
{
	type Item = Temporal<T::Item>;
	fn next(&mut self) -> Option<Self::Item> {
		self.iter.next()
			.map(|x| Temporal::new(x, self.time))
	}
	fn size_hint(&self) -> (usize, Option<usize>) {
		self.iter.size_hint()
	}
}

/// Multidimensional change over time.
pub trait TemporalVector<const SIZE: usize> {
	type Kind: Vector<SIZE, Output: FluxKind> + 'static;
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<U, D>(&self, other: &Temporal<U>, cmp: Ordering, dis: &Temporal<D>)
		-> <Temporal<Self::Kind> as WhenDis<SIZE, U::Kind, D::Kind>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Temporal<Self::Kind>: WhenDis<SIZE, U::Kind, D::Kind>,
	;
	
	/// Ranges when the distance to another vector is equal to X.
	fn when_dis_eq<U, D>(&self, other: &Temporal<U>, dis: &Temporal<D>)
		-> <Temporal<Self::Kind> as WhenDisEq<SIZE, U::Kind, D::Kind>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Temporal<Self::Kind>: WhenDisEq<SIZE, U::Kind, D::Kind>,
	;
	
	/// Ranges when a component is above/below/equal to another flux.
	fn when_index<U>(&self, index: usize, cmp: Ordering, other: &Temporal<U>)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as When<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: When<U::Kind>
	;
	
	/// Times when a component is equal to another flux.
	fn when_index_eq<U>(&self, index: usize, other: &Temporal<U>)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<U::Kind>
	;
	
	/// Ranges when the distance to another vector is above/below/equal to a constant.
	fn when_dis_constant<U, C>(&self, other: &Temporal<U>, cmp: Ordering, dis: C)
		-> <Temporal<Self::Kind> as WhenDis<SIZE, U::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		C: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Temporal<Self::Kind>: WhenDis<SIZE, U::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	{
		self.when_dis(other, cmp, &Temporal::new(Constant::from(C::into_linear(dis)), Time::ZERO))
	}
	
	/// Ranges when the distance to another vector is equal to a constant.
	fn when_dis_eq_constant<U, C>(&self, other: &Temporal<U>, dis: C)
		-> <Temporal<Self::Kind> as WhenDisEq<SIZE, U::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		C: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Temporal<Self::Kind>: WhenDisEq<SIZE, U::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	{
		self.when_dis_eq(other, &Temporal::new(Constant::from(C::into_linear(dis)), Time::ZERO))
	}
	
	/// Ranges when a component is above/below/equal to a constant.
	fn when_index_constant<C>(&self, index: usize, cmp: Ordering, other: C)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		C: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	{
		self.when_index(index, cmp, &Temporal::new(Constant::from(C::into_linear(other)), Time::ZERO))
	}
	
	/// Times when a component is equal to a constant.
	fn when_index_eq_constant<C>(&self, index: usize, other: C)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		C: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	{
		self.when_index_eq(index, &Temporal::new(Constant::from(C::into_linear(other)), Time::ZERO))
	}
	
	// !!!
	// - to rotate line by a fixed angle, multiply axial polynomials by Scalar?
	// - to fit a line segment, filter predicted times.
	// - for non-fixed angle line segments, use a double distance check?
	// - rotating point-line may be handleable iteratively, find the bounds in
	//   which the roots may be and iterate through it.
}

impl<T, const SIZE: usize> TemporalVector<SIZE> for Temporal<T>
where
	T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
{
	type Kind = T::Kind;
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<U, D>(&self, other: &Temporal<U>, cmp: Ordering, dis: &Temporal<D>)
		-> <Temporal<Self::Kind> as WhenDis<SIZE, U::Kind, D::Kind>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Temporal<Self::Kind>: WhenDis<SIZE, U::Kind, D::Kind>,
	{
		self.poly(self.time)
			.when_dis(other.poly(self.time), cmp, dis.poly(self.time))
	}
	
	/// Ranges when the distance to another vector is equal to X.
	fn when_dis_eq<U, D>(&self, other: &Temporal<U>, dis: &Temporal<D>)
		-> <Temporal<Self::Kind> as WhenDisEq<SIZE, U::Kind, D::Kind>>::Pred
	where
		U: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Temporal<Self::Kind>: WhenDisEq<SIZE, U::Kind, D::Kind>,
	{
		self.poly(self.time)
			.when_dis_eq(other.poly(self.time), dis.poly(self.time))
	}
	
	/// Ranges when a component is above/below/equal to another flux.
	fn when_index<U>(&self, index: usize, cmp: Ordering, other: &Temporal<U>)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as When<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: When<U::Kind>
	{
		self.poly(self.time).index(index)
			.when(cmp, other.poly(self.time))
	}
	
	/// Times when a component is equal to another flux.
	fn when_index_eq<U>(&self, index: usize, other: &Temporal<U>)
		-> <Temporal<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<U::Kind>>::Pred
	where
		U: Flux,
		Temporal<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<U::Kind>
	{
		self.poly(self.time).index(index)
			.when_eq(other.poly(self.time))
	}
}