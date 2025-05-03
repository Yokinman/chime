//! Utilities for describing how a type changes over time.

#![recursion_limit = "2048"]

pub mod change;
pub mod temporal;
pub mod time;
pub mod linear;
pub mod poly;
pub mod pred;
pub mod exp;

mod impls;

use linear::*;
use change::*;

use time::Time;

pub use chime_flux_proc_macro::{Flux, ToMoment, ToMomentMut};

#[cfg(feature = "bevy")]
pub use chime_flux_proc_macro::{TemporalComponent, TemporalResource};

pub use temporal::Temporal;

/// Context for a `bevy_time::Time`.
#[cfg_attr(not(feature = "bevy"), doc(hidden))]
#[derive(Default)]
pub struct Chime;

/// Immutable moment-in-time interface for [`Temporal::moment`].
pub struct Moment<'a, T: ToMoment> {
	pub(crate) moment: T::Moment<'a>,
	pub(crate) borrow: std::marker::PhantomData<&'a T>,
}

mod _moment_ref_impls {
	use std::fmt::{Debug, Display, Formatter};
	use std::ops::Deref;
	use super::{ToMoment, Moment};
	
	impl<'a, T: ToMoment> Deref for Moment<'a, T> {
		type Target = T::Moment<'a>;
		fn deref(&self) -> &Self::Target {
			&self.moment
		}
	}
	
	impl<'a, T: ToMoment> Debug for Moment<'a, T>
	where
		T::Moment<'a>: Debug
	{
		fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
			<T::Moment<'a> as Debug>::fmt(self, f)
		}
	}
	
	impl<'a, T: ToMoment> Display for Moment<'a, T>
	where
		T::Moment<'a>: Display
	{
		fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
			<T::Moment<'a> as Display>::fmt(self, f)
		}
	}
}

/// Mutable moment-in-time interface for [`Temporal::moment_mut`].
pub struct MomentMut<'a, T: ToMomentMut> {
	pub(crate) moment: T::MomentMut<'a>,
	pub(crate) borrow: std::marker::PhantomData<&'a mut T>,
}

mod _moment_mut_impls {
	use std::fmt::{Debug, Display, Formatter};
	use std::ops::{Deref, DerefMut};
	use super::{ToMomentMut, MomentMut};
	
	impl<'a, T: ToMomentMut> Deref for MomentMut<'a, T> {
		type Target = T::MomentMut<'a>;
		fn deref(&self) -> &Self::Target {
			&self.moment
		}
	}
	
	impl<T: ToMomentMut> DerefMut for MomentMut<'_, T> {
		fn deref_mut(&mut self) -> &mut Self::Target {
			&mut self.moment
		}
	}
	
	impl<'a, T: ToMomentMut> Debug for MomentMut<'a, T>
	where
		T::MomentMut<'a>: Debug
	{
		fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
			<T::MomentMut<'a> as Debug>::fmt(self, f)
		}
	}
	
	impl<'a, T: ToMomentMut> Display for MomentMut<'a, T>
	where
		T::MomentMut<'a>: Display
	{
		fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
			<T::MomentMut<'a> as Display>::fmt(self, f)
		}
	}
}

#[cfg(feature = "bevy")]
mod bevy_moment {
	use std::ops::{Deref, DerefMut};
	use bevy_ecs::archetype::Archetype;
	use bevy_ecs::component::{Component, ComponentId, Components, Tick};
	use bevy_ecs::entity::Entity;
	use bevy_ecs::query::{FilteredAccess, QueryData, ReadOnlyQueryData, WorldQuery};
	use bevy_ecs::storage::{Table, TableRow};
	use bevy_ecs::system::{Res, ResMut, Resource, SystemMeta, SystemParam};
	use bevy_ecs::world::{World, unsafe_world_cell::UnsafeWorldCell};
	use super::*;
	
	type ChimeTime = bevy_time::Time<Chime>;
	
	type Ref<'b, M> = &'b Temporal<M>;
	type Mut<'b, M> = &'b mut Temporal<M>;
	
	/// SAFETY: `Self` is the same as `Self::ReadOnly`.
	unsafe impl<M> QueryData for Moment<'_, M>
	where
		M: ToMoment,
		Temporal<M>: Component,
	{
		type ReadOnly = Self;
	}
	
	/// SAFETY: access is read only.
	unsafe impl<M> ReadOnlyQueryData for Moment<'_, M>
	where
		M: ToMoment,
		Temporal<M>: Component,
	{}
	
	/// SAFETY: access of `Moment<T>` is a subset of `MomentMut<T>`.
	unsafe impl<'b, M> QueryData for MomentMut<'b, M>
	where
		M: ToMomentMut,
		Temporal<M>: Component,
	{
		type ReadOnly = Moment<'b, M>;
	}
	
	/// SAFETY:
	/// Most safety guarantees are inherited from `impl QueryData for &T`.
	/// This additionally requests a `Time` resource and takes care of adding
	/// read access for that resource in `update_component_access`.
	unsafe impl<'b, M> WorldQuery for Moment<'b, M>
	where
		M: ToMoment,
		Temporal<M>: Component,
	{
		type Item<'a> = Moment<'a, M>;
		type Fetch<'a> = (Time, <Ref<'b, M> as WorldQuery>::Fetch<'a>);
		type State = (ComponentId, <Ref<'b, M> as WorldQuery>::State);
		fn shrink<'wlong: 'wshort, 'wshort>(item: Self::Item<'wlong>) -> Self::Item<'wshort> {
			unsafe {
				// SAFETY: This is not safe, but I'm refactoring and this is the
				// only problem so I'll fix it later if it causes issues. IDGAF!!!!
				std::mem::transmute(item)
			}
		}
		fn shrink_fetch<'wlong: 'wshort, 'wshort>(fetch: Self::Fetch<'wlong>) -> Self::Fetch<'wshort> {
			fetch
		}
		unsafe fn init_fetch<'w>(world: UnsafeWorldCell<'w>, (_, state): &Self::State, last_run: Tick, this_run: Tick) -> Self::Fetch<'w> {
			// !!! For parallelization, `last_run` could be used as a nanosecond
			// offset. It's a hack but it's also the most usable option.
			let time = world.get_resource::<ChimeTime>()
				.expect("bevy_time::Time resource should exist in the world")
				.elapsed();
			(time, <Ref<'b, M> as WorldQuery>::init_fetch(world, state, last_run, this_run))
		}
		const IS_DENSE: bool = <Ref<'b, M> as WorldQuery>::IS_DENSE;
		unsafe fn set_archetype<'w>((_, fetch): &mut Self::Fetch<'w>, (_, state): &Self::State, archetype: &'w Archetype, table: &'w Table) {
			<Ref<'b, M> as WorldQuery>::set_archetype(fetch, state, archetype, table)
		}
		unsafe fn set_table<'w>((_, fetch): &mut Self::Fetch<'w>, (_, state): &Self::State, table: &'w Table) {
			<Ref<'b, M> as WorldQuery>::set_table(fetch, state, table)
		}
		unsafe fn fetch<'w>((time, fetch): &mut Self::Fetch<'w>, entity: Entity, table_row: TableRow) -> Self::Item<'w> {
			<Ref<'b, M> as WorldQuery>::fetch(fetch, entity, table_row)
				.moment(*time)
		}
		fn update_component_access((time_id, state): &Self::State, access: &mut FilteredAccess<ComponentId>) {
	        assert!(
	            !access.access().has_resource_write(*time_id),
	            "&{} conflicts with a previous access in this query. Shared access cannot coincide with exclusive access.",
	                std::any::type_name::<ChimeTime>(),
	        );
	        access.access_mut().add_resource_read(*time_id);
			<Ref<'b, M> as WorldQuery>::update_component_access(state, access)
		}
		fn init_state(world: &mut World) -> Self::State {
			(world.init_resource::<ChimeTime>(),
				<Ref<'b, M> as WorldQuery>::init_state(world))
		}
		fn get_state(components: &Components) -> Option<Self::State> {
			components.resource_id::<ChimeTime>()
				.zip(<Ref<'b, M> as WorldQuery>::get_state(components))
		}
		fn matches_component_set((_, state): &Self::State, set_contains_id: &impl Fn(ComponentId) -> bool) -> bool {
			<Ref<'b, M> as WorldQuery>::matches_component_set(state, set_contains_id)
		}
	}
	
	/// SAFETY:
	/// Most safety guarantees are inherited from `impl QueryData for &mut T`.
	/// This additionally requests a `Time` resource and takes care of adding
	/// read access for that resource in `update_component_access`.
	unsafe impl<'b, M> WorldQuery for MomentMut<'b, M>
	where
		M: ToMomentMut,
		Temporal<M>: Component,
	{
		type Item<'a> = MomentMut<'a, M>;
		type Fetch<'a> = (Time, <Mut<'b, M> as WorldQuery>::Fetch<'a>);
		type State = (ComponentId, <Mut<'b, M> as WorldQuery>::State);
		fn shrink<'wlong: 'wshort, 'wshort>(item: Self::Item<'wlong>) -> Self::Item<'wshort> {
			unsafe {
				// SAFETY: This is not safe, but I'm refactoring and this is the
				// only problem so I'll fix it later if it causes issues. IDGAF!!!!
				std::mem::transmute(item)
			}
		}
		fn shrink_fetch<'wlong: 'wshort, 'wshort>(fetch: Self::Fetch<'wlong>) -> Self::Fetch<'wshort> {
			fetch
		}
		unsafe fn init_fetch<'w>(world: UnsafeWorldCell<'w>, (_, state): &Self::State, last_run: Tick, this_run: Tick) -> Self::Fetch<'w> {
			let time = world.get_resource::<ChimeTime>()
				.expect("bevy_time::Time resource should exist in the world")
				.elapsed();
			(time, <Mut<'b, M> as WorldQuery>::init_fetch(world, state, last_run, this_run))
		}
		const IS_DENSE: bool = <Mut<'b, M> as WorldQuery>::IS_DENSE;
		unsafe fn set_archetype<'w>((_, fetch): &mut Self::Fetch<'w>, (_, state): &Self::State, archetype: &'w Archetype, table: &'w Table) {
			<Mut<'b, M> as WorldQuery>::set_archetype(fetch, state, archetype, table)
		}
		unsafe fn set_table<'w>((_, fetch): &mut Self::Fetch<'w>, (_, state): &Self::State, table: &'w Table) {
			<Mut<'b, M> as WorldQuery>::set_table(fetch, state, table)
		}
		unsafe fn fetch<'w>((time, fetch): &mut Self::Fetch<'w>, entity: Entity, table_row: TableRow) -> Self::Item<'w> {
			<Mut<'b, M> as WorldQuery>::fetch(fetch, entity, table_row)
				.into_inner()
				.moment_mut(*time)
		}
		fn update_component_access((time_id, state): &Self::State, access: &mut FilteredAccess<ComponentId>) {
	        assert!(
	            !access.access().has_resource_write(*time_id),
	            "&{} conflicts with a previous access in this query. Shared access cannot coincide with exclusive access.",
	                std::any::type_name::<ChimeTime>(),
	        );
	        access.access_mut().add_resource_read(*time_id);
			<Mut<'b, M> as WorldQuery>::update_component_access(state, access);
		}
		fn init_state(world: &mut World) -> Self::State {
			(world.init_resource::<ChimeTime>(),
				<Mut<'b, M> as WorldQuery>::init_state(world))
		}
		fn get_state(components: &Components) -> Option<Self::State> {
			components.resource_id::<ChimeTime>()
				.zip(<Mut<'b, M> as WorldQuery>::get_state(components))
		}
		fn matches_component_set((_, state): &Self::State, set_contains_id: &impl Fn(ComponentId) -> bool) -> bool {
			<Mut<'b, M> as WorldQuery>::matches_component_set(state, set_contains_id)
		}
	}
	
	/// `SystemParam` for fetching a resource at the current `Time<Chime>`.
	pub struct ResMoment<'w, M: ToMoment> {
		inner: Moment<'w, M>,
	}
	
	impl<'w, M: ToMoment> Deref for ResMoment<'w, M> {
		type Target = Moment<'w, M>;
		fn deref(&self) -> &Self::Target {
			&self.inner
		}
	}
	
	/// `SystemParam` for fetching a mutable resource at the current `Time<Chime>`.
	pub struct ResMomentMut<'w, M: ToMomentMut> {
		inner: MomentMut<'w, M>,
	}
	
	impl<'w, M: ToMomentMut> Deref for ResMomentMut<'w, M> {
		type Target = MomentMut<'w, M>;
		fn deref(&self) -> &Self::Target {
			&self.inner
		}
	}
	
	impl<M: ToMomentMut> DerefMut for ResMomentMut<'_, M> {
		fn deref_mut(&mut self) -> &mut Self::Target {
			&mut self.inner
		}
	}
	
	/// SAFETY:
	/// Safety guarantees are inherited from `impl SystemParam for Res<T>`.
	unsafe impl<'w, M> SystemParam for ResMoment<'w, M>
	where
		M: ToMoment,
		Temporal<M>: Resource,
	{
		type State = (<Res<'w, ChimeTime> as SystemParam>::State, <Res<'w, Temporal<M>> as SystemParam>::State);
		type Item<'world, 'state> = ResMoment<'world, M>;
		fn init_state(world: &mut World, system_meta: &mut SystemMeta) -> Self::State {
			let time_state = <Res<'w, ChimeTime> as SystemParam>::init_state(world, system_meta);
			let res_state = <Res<'w, Temporal<M>> as SystemParam>::init_state(world, system_meta);
			(time_state, res_state)
		}
		unsafe fn get_param<'world, 'state>((time_state, res_state): &'state mut Self::State, system_meta: &SystemMeta, world: UnsafeWorldCell<'world>, change_tick: Tick) -> Self::Item<'world, 'state> {
			let time = <Res<'w, ChimeTime> as SystemParam>::get_param(time_state, system_meta, world, change_tick)
				.elapsed();
			let inner = <Res<'w, Temporal<M>> as SystemParam>::get_param(res_state, system_meta, world, change_tick)
				.into_inner()
				.moment(time);
			ResMoment { inner }
		}
	}
	
	/// SAFETY:
	/// Safety guarantees are inherited from `impl SystemParam for Res<T>` and
	/// `impl SystemParam for ResMut<T>`.
	unsafe impl<'w, M> SystemParam for ResMomentMut<'w, M>
	where
		M: ToMomentMut,
		Temporal<M>: Resource,
	{
		type State = (<Res<'w, ChimeTime> as SystemParam>::State, <ResMut<'w, Temporal<M>> as SystemParam>::State);
		type Item<'world, 'state> = ResMomentMut<'world, M>;
		fn init_state(world: &mut World, system_meta: &mut SystemMeta) -> Self::State {
			let time_state = <Res<'w, ChimeTime> as SystemParam>::init_state(world, system_meta);
			let res_state = <ResMut<'w, Temporal<M>> as SystemParam>::init_state(world, system_meta);
			(time_state, res_state)
		}
		unsafe fn get_param<'world, 'state>((time_state, res_state): &'state mut Self::State, system_meta: &SystemMeta, world: UnsafeWorldCell<'world>, change_tick: Tick) -> Self::Item<'world, 'state> {
			let time = <Res<'w, ChimeTime> as SystemParam>::get_param(time_state, system_meta, world, change_tick)
				.elapsed();
			let inner = <ResMut<'w, Temporal<M>> as SystemParam>::get_param(res_state, system_meta, world, change_tick)
				.into_inner()
				.moment_mut(time);
			ResMomentMut { inner }
		}
	}
}

#[cfg(feature = "bevy")]
pub use bevy_moment::{ResMoment, ResMomentMut};

/// Change applied per unit of time.
/// 
/// Generally constructed using [`Flux::per`] and applied to an
/// [accumulator](ChangeAccum) in the [`Flux::change`] method.
#[derive(Copy, Clone, Debug, Default, Ord, PartialOrd, Eq, PartialEq)]
pub struct Rate<T> {
	pub amount: T,
	pub unit: Time,
}

mod _rate_impls {
	use crate::{Flux, ToMoment, ToMomentMut};
	use crate::change::{ApplyChange, Change, ChangeUp};
	use super::Rate;
	
	impl<T> Rate<T> {
		pub fn as_ref(&self) -> Rate<&T> {
			Rate {
				amount: &self.amount,
				unit: self.unit,
			}
		}
	}
	
	impl<T, F, const OP: char> ApplyChange<OP, Rate<F>> for T
	where
		T: ApplyChange<OP, <F::Change as ChangeUp<OP>>::Up>,
		F: Flux<Change: ChangeUp<OP>>,
	{
		type Output = T::Output;
		fn apply_change(self, rhs: Rate<F>) -> Self::Output {
			let change = rhs.amount.change()
				.up(rhs.amount.basis())
				.scale(crate::linear::Linear::from_f64(rhs.unit.as_secs_f64().recip()));
			self.apply_change(change)
		}
	}
	
	impl<'a, T, F, const OP: char> ApplyChange<OP, &'a Rate<F>> for T
	where
		T: ApplyChange<OP, Rate<&'a F>>,
	{
		type Output = T::Output;
		fn apply_change(self, rhs: &'a Rate<F>) -> Self::Output {
			self.apply_change(rhs.as_ref())
		}
	}
	
	impl<T: ToMoment> ToMoment for Rate<T> {
		type Moment<'a> = Rate<T::Moment<'a>> where Self: 'a;
		fn to_moment(&self, time: f64) -> Self::Moment<'_> {
			Rate {
				amount: self.amount.to_moment(time),
				unit: self.unit,
			}
		}
	}
	
	impl<T: ToMomentMut> ToMomentMut for Rate<T> {
		type MomentMut<'a> = Rate<T::MomentMut<'a>> where Self: 'a;
		fn to_moment_mut(&mut self, time: f64) -> Self::MomentMut<'_> {
			Rate {
				amount: self.amount.to_moment_mut(time),
				unit: self.unit,
			}
		}
	}
}

/// Defining values that change over time.
/// 
/// This is similar to [`ToMoment`] in that both describe change over time, but
/// specific to abstracting the timeline. User types will often implement all of
/// `Flux`, [`ToMoment`], and [`ToMomentMut`].
/// 
/// # Describing a Change
/// 
/// For simple cases, derive `Flux`:
/// 
/// ```
/// use chime::Flux;
/// use chime::time::SEC;
/// 
/// // pos = pos + spd
/// #[derive(Flux)]
/// struct Position {
///     #[basis]
///     pos: f64,
///     #[change(add_per(SEC))]
///     spd: Speed,
/// }
/// 
/// // spd = spd + acc - fric
/// #[derive(Flux)]
/// struct Speed {
///     #[basis]
///     spd: f64,
///     #[change(add_per(SEC))]
///     acc: f64,
///     #[change(sub_per(SEC))]
///     fric: f64,
/// }
/// ```
/// 
/// # Evaluating a State
/// 
/// ```
/// use chime::Flux;
/// use chime::poly::Poly;
/// use chime::time::SEC;
/// 
/// #[derive(Flux)]
/// struct Num {
///     #[basis]
///     num: f64,
///     #[change(add_per(SEC))]
///     add: f64,
/// }
/// 
/// let poly = Num { num: 3., add: 1.5 }.to_poly();
/// assert_eq!(poly.eval(0.), 3.);
/// assert_eq!(poly.eval(1.), 4.5);
/// assert_eq!(poly.eval(2.), 6.);
/// ```
pub trait Flux {
	/// The type of value (e.g. `f64`, `f32`, `[f64; N]`, etc.).
	type Basis: Basis;
	
	/// The type of change over time (e.g. [`Nil<T>`](constant::Nil),
	/// [`Sum<T, D>`](sum::Sum), etc.).
	type Change: Change<Basis = Self::Basis>;
	
	/// The initial value.
	fn basis(&self) -> Self::Basis;
	
	/// The applied change over time.
	fn change(&self) -> Self::Change;
	
	/// Conversion into a standard representation.
	fn to_poly(&self) -> <Self::Change as Change>::Poly {
		self.change().into_poly(self.basis())
	}
	
	/// Used with [`Self::per`] to manually implement [`Self::change`].
	fn accum(&self) -> ChangeAccum<constant::Nil<Self::Basis>> {
		ChangeAccum::default()
	}
	
	/// Used to construct a [`Rate`] for convenient change-over-time operations.
	/// 
	/// `1 + 2.per(time_unit::SEC)` 
	fn per(&self, unit: Time) -> Rate<&Self> where Self: Sized {
		Rate {
			amount: &self,
			unit,
		}
	}
}

/// Types that represent a timeline of moments.
/// 
/// This is similar to [`Flux`] in that both describe change over time, but
/// specific to interfacing with the timeline. User types will often implement
/// all of [`Flux`], `ToMoment`, and [`ToMomentMut`].
/// 
/// Data structures such as `HashMap<T> where T: Flux + ToMoment` only
/// implement `ToMoment`, since they represent multiple values that change
/// over time.
pub trait ToMoment {
	/// The interface for a moment in the timeline.
	type Moment<'a> where Self: 'a;
	
	/// Produces a moment in the timeline.
	/// 
	/// This generally returns an owned value. However, the value should be
	/// treated as a reference, as modifying it has no effect on the timeline.
	/// This is enforced through [`Temporal::moment`] ([`Moment`]).
	fn to_moment(&self, time: f64) -> Self::Moment<'_>;
}

/// Types that represent a mutable timeline of moments.
/// 
/// This is similar to [`Flux`] in that both describe change over time, but
/// specific to interfacing with the timeline. User types will often implement
/// all of [`Flux`], [`ToMoment`], and `ToMomentMut`.
pub trait ToMomentMut: ToMoment {
	/// The mutable interface for a moment in the timeline.
	type MomentMut<'a> where Self: 'a;
	
	/// Produces a mutable moment in the timeline.
	/// 
	/// In general, this should permanently shift the basis of the value and
	/// return the moment by reference, unlike [`ToMoment::to_moment`]. This
	/// is enforced through [`Temporal::moment_mut`] ([`MomentMut`]).
	fn to_moment_mut(&mut self, time: f64) -> Self::MomentMut<'_>;
}

#[cfg(test)]
mod tests {
	use std::ops::Deref;
	use super::*;
	use crate::time;
	use crate::time::{SEC, TimeRanges};
	use crate::change::{Nil, Sum2};
	use crate::pred::Prediction;
	use std::cmp::Ordering;
	
	use crate as chime;
	
	#[derive(Clone, Debug, Default, ToMoment, ToMomentMut)]
	struct Pos {
		#[basis]
		value: f64,
		spd: Spd,
		misc: Vec<Spd>,
	}
	
	impl Flux for Pos {
		type Basis = f64;
		type Change = Sum2<Nil<f64>, Sum2<Nil<f64>, Sum2<Nil<f64>, Sum2<Nil<f64>, Nil<f64>>>>>;
		fn basis(&self) -> Self::Basis {
			self.value
		}
		fn change(&self) -> Self::Change {
			let mut accum = self.accum().into_change() + self.spd.per(SEC);
			for spd in &self.misc {
				accum = accum + spd.per(SEC);
			}
			accum
		}
	}
	
	#[derive(Clone, Debug, Default, Flux, ToMoment, ToMomentMut)]
	struct Spd {
		#[basis]
		value: f64,
		#[change(sub_per(SEC))]
		fric: Fric,
		#[change(add_per(SEC))]
		accel: Accel,
	}
	
	#[derive(Clone, Debug, Default, Flux, ToMoment, ToMomentMut)]
	struct Fric {
		#[basis]
		value: f64,
	}
	
	#[derive(Clone, Debug, Default, Flux, ToMoment, ToMomentMut)]
	struct Accel {
		#[basis]
		value: f64,
		#[change(add_per(SEC))]
		jerk: Jerk,
	}
	
	#[derive(Clone, Debug, Default, Flux, ToMoment, ToMomentMut)]
	struct Jerk {
		#[basis]
		value: f64,
		#[change(add_per(SEC))]
		snap: Snap,
	}
	
	#[derive(Clone, Debug, Default, Flux, ToMoment, ToMomentMut)]
	struct Snap {
		#[basis]
		value: f64,
	}
	
	impl Deref for Pos {
		type Target = f64;
		fn deref(&self) -> &Self::Target {
			&self.value
		}
	}
	
	fn position() -> Temporal<Pos> {
		let mut pos = Temporal::from(Pos {
			value: 32.0, 
			spd: Spd {
				value: -4.4075 / 3.0,
				fric: Fric { value: 3.5 },
				accel: Accel {
					value: 2.0725 / 3.0,
					jerk: Jerk {
						value: 0.385,
						snap: Snap { value: -0.01 },
					},
				},
			},
			misc: Vec::new(),
		});
		*pos = pos.to_moment(10*SEC);
		pos.time = 10*SEC;
		pos
	}
	
	macro_rules! assert_times {
		($times:expr, []) => {
			assert_time_ranges!($times, [(Time::ZERO, Time::MAX)])
		};
		($times:expr, $cmp_times:expr) => {
			assert_time_ranges!(
				$times,
				$cmp_times.into_iter()
					.map(|x| (x, x))
					.collect::<Vec<_>>()
			)
			// let times: Times = $times;
			// let cmp_times = Times::new($cmp_times);
			// assert_eq!(
			// 	times.clone().count(),
			// 	cmp_times.clone().count(),
			// 	"a: {:?}, b: {:?}",
			// 	times.collect::<Box<[_]>>(),
			// 	cmp_times.collect::<Box<[_]>>(),
			// );
			// for (a, b) in times.zip(cmp_times) {
			// 	let time = ((a.max(b) - b.min(a)).as_secs_f64() * 60.).floor();
			// 	assert_eq!(time, 0., "a: {:?}, b: {:?}", a, b);
			// }
		};
	}
	
	macro_rules! assert_time_ranges {
		($ranges:expr, $cmp_ranges:expr) => {{
			let ranges = Vec::<(Time, Time)>::from_iter($ranges
				.into_ranges(Time::ZERO)
				.inclusive());
			let cmp_ranges = Vec::<(Time, Time)>::from_iter($cmp_ranges);
			assert_eq!(
				ranges.len(),
				cmp_ranges.len(),
				"a: {:?}, b: {:?}",
				ranges, cmp_ranges,
			);
			for ((a, b), (x, y)) in ranges.into_iter().zip(cmp_ranges) {
				let a_time = ((a.max(x) - x.min(a)).as_secs_f64() * 60.).floor();
				let b_time = ((b.max(y) - y.min(b)).as_secs_f64() * 60.).floor();
				assert_eq!(a_time, 0., "a: {:?}, x: {:?}", a, x);
				assert_eq!(b_time, 0., "b: {:?}, y: {:?}", b, y);
			}
		}};
	}
	
	macro_rules! assert_poly {
		($poly:expr, $cmp_poly:expr) => {{
			let poly = $poly;
			let cmp_poly = $cmp_poly;
			let dif_poly = poly.sub_poly(cmp_poly);
			for i in 0..100 {
				assert_eq!((dif_poly.eval(SEC * i) * 2e8).round(), 0.,
					"\npoly: {poly:?}\ncmp_poly: {cmp_poly:?}");
			}
		}};
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		assert_eq!(pos.basis_time(), 10*SEC);
		
		 // Values:
		assert_eq!(pos.moment(0*SEC).round(), 32.);
		assert_eq!(pos.moment(10*SEC).round(), -63.);
		assert_eq!(pos.moment(20*SEC).round(), -113.);
		assert_eq!(pos.moment(100*SEC).round(), 8339.);
		assert_eq!(pos.moment(200*SEC).round(), -209778.);
		
		 // Update:
		assert_eq!(pos.moment_mut(20*SEC).round(), -113.);
		assert_eq!(pos.moment(100*SEC).round(), 8339.);
		assert_eq!(pos.moment(200*SEC).round(), -209778.);
		assert_eq!(pos.moment_mut(100*SEC).round(), 8339.);
		assert_eq!(pos.moment(200*SEC).round(), -209778.);
	}
	
	#[test]
	fn poly() {
		use num_traits::Pow;
		let mut pos = position();
		assert_poly!(
			pos.to_poly(),
			Temporal::new(symb_poly::Invar(crate::constant::Constant2(-63.15))
				+ symb_poly::Invar(crate::constant::Constant2(-11.9775)) * <symb_poly::Var>::default()
				+ symb_poly::Invar(crate::constant::Constant2(0.270416666666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P2>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(0.0475)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P3>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(-0.00041666666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P4>::default()))
			, 10*SEC)
		);
		for _ in 0..2 {
			pos.moment_mut(20*SEC);
			assert_poly!(
				pos.to_poly(),
				Temporal::new(symb_poly::Invar(crate::constant::Constant2(-112.55))
					+ symb_poly::Invar(crate::constant::Constant2(6.0141666666666666666)) * <symb_poly::Var>::default()
					+ symb_poly::Invar(crate::constant::Constant2(1.4454166666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P2>::default()))
					+ symb_poly::Invar(crate::constant::Constant2(0.0308333333333333333)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P3>::default()))
					+ symb_poly::Invar(crate::constant::Constant2(-0.0004166666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P4>::default()))
				, 20*SEC)
			);
		}
		pos.moment_mut(0*SEC);
		assert_poly!(
			pos.to_poly(),
			Temporal::new(symb_poly::Invar(crate::constant::Constant2(32.))
					+ symb_poly::Invar(crate::constant::Constant2(-1.4691666666666666666)) * <symb_poly::Var>::default()
					+ symb_poly::Invar(crate::constant::Constant2(-1.4045833333333333333)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P2>::default()))
					+ symb_poly::Invar(crate::constant::Constant2(0.0641666666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P3>::default()))
					+ symb_poly::Invar(crate::constant::Constant2(-0.0004166666666666666)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P4>::default()))
			, 0*SEC)
		);
		
		assert_poly!(
			Temporal::new(symb_poly::Invar(crate::constant::Constant2(-0.8_f64))
				+ symb_poly::Invar(crate::constant::Constant2(-2.7)) * <symb_poly::Var>::default()
				+ symb_poly::Invar(crate::constant::Constant2(3.4)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P2>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(2.8)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P3>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(0.3)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P4>::default()))
			, 4000*time::MILLISEC).at_time(320*time::MILLISEC),
			Temporal::new(symb_poly::Invar(crate::constant::Constant2(-29.341750272_f64))
				+ symb_poly::Invar(crate::constant::Constant2(26.2289216)) * <symb_poly::Var>::default()
				+ symb_poly::Invar(crate::constant::Constant2(-3.13568)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P2>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(-1.616)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P3>::default()))
				+ symb_poly::Invar(crate::constant::Constant2(0.3)) * <symb_poly::Var>::default().pow(symb_poly::Invar(symb_poly::Num::<typenum::P4>::default()))
			, 320*time::MILLISEC)
		);
	}
	
	#[test]
	fn when() {
		let pos = position();
		let acc = Temporal::from(Accel {
			value: 0.3,
			jerk: Jerk {
				value: 0.4,
				snap: Snap { value: -0.01 },
			},
		});
		
		assert_time_ranges!(pos.to_poly().when(Ordering::Greater, acc.to_poly()), [
			(Time::ZERO, Time::from_secs_f64(4.56)),
			(Time::from_secs_f64(26.912), Time::from_secs_f64(127.394))
		]);
		assert_times!(pos.to_poly().when_eq(acc.to_poly()), [
			Time::from_secs_f64(4.56),
			Time::from_secs_f64(26.912),
			Time::from_secs_f64(127.394)
		]);
	}
	
	#[test]
	fn when_at_mut() {
		let mut a_pos = position();
		let mut b_pos = position();
		
		 // Check Before:
		assert_time_ranges!(a_pos.to_poly().when(Ordering::Greater, b_pos.to_poly()), []);
		assert_time_ranges!(a_pos.to_poly().when(Ordering::Equal, b_pos.to_poly()), [(Time::ZERO, Time::MAX)]);
		assert_times!(a_pos.to_poly().when_eq(b_pos.to_poly()), []);
		a_pos.moment_mut(20*SEC);
		
		 // Apply Changes:
		b_pos.moment_mut(0*SEC).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.moment_mut(0*SEC).misc.push(Spd {
			value: 12.,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.moment_mut(10*SEC).value -= 100.0;
		
		 // Check After:
		assert_eq!(
			Vec::from_iter(a_pos.to_poly().when(Ordering::Greater, b_pos.to_poly())
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(0*SEC, 8*SEC - time::NANOSEC),
				(50*SEC /*+ time::NANOSEC*/, Time::MAX)
			]
		);
		assert_times!(a_pos.to_poly().when_eq(b_pos.to_poly()), [8*SEC, 50*SEC]);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq, Copy, Clone, Debug, Flux, ToMoment, ToMomentMut)]
		struct Pos {
			#[basis]
			value: Iso<f64, i64>,
			#[change(add_per(time::MINUTE))]
			spd: Spd,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone, Debug, Flux, ToMoment, ToMomentMut)]
		struct Spd {
			#[basis]
			value: Iso<f64, i64>,
			#[change(add_per(SEC))]
			acc: Acc,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone, Debug, Flux, ToMoment, ToMomentMut)]
		struct Acc {
			#[basis]
			value: Iso<f64, i64>,
		}

		todo!("re-add Roots for Iso impls");
		// let a_pos = Temporal::from([
		// 	Pos { value: Iso::new(3), spd: Spd { value: Iso::new(300), acc: Acc { value: Iso::new(-240) } } },
		// 	Pos { value: Iso::new(-4), spd: Spd { value: Iso::new(120), acc: Acc { value: Iso::new(1080) } } }
		// ]);
		// let b_pos = Temporal::from([
		// 	Pos { value: Iso::new(8), spd: Spd { value: Iso::new(330), acc: Acc { value: Iso::new(-300) } } },
		// 	Pos { value: Iso::new(4), spd: Spd { value: Iso::new(600), acc: Acc { value: Iso::new(720) } } }
		// ]);
		//
		// let dis = Temporal::from(Spd { value: Iso::new(10), acc: Acc { value: Iso::new(0) } });
		//
		// assert_time_ranges!(
		// 	a_pos.to_poly().when_dis(b_pos.to_poly(), Ordering::Less, dis.to_poly()),
		// 	// https://www.desmos.com/calculator/spxyoloyx9
		// 	[
		// 		(Time::ZERO, Time::from_secs_f64(0.0823337)),
		// 		(Time::from_secs_f64(2.246506069), Time::from_secs_f64(4.116193987)),
		// 	]
		// );
		//
		// let b_pos = Temporal::new(b_pos.to_moment(Time::ZERO), SEC);
		// assert_time_ranges!(
		// 	a_pos.to_poly().when_dis_eq_constant(b_pos.to_poly(), 2.into()),
		// 	[
		// 		(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
		// 		(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191)),
		// 	]
		// );
		// assert_time_ranges!(
		// 	a_pos.to_poly().when_dis_constant(b_pos.to_poly(), Ordering::Equal, 2.into()),
		// 	[
		// 		(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
		// 		(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191))
		// 	]
		// );
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}