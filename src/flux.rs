//! Utilities for describing how a type changes over time.

pub mod time;
pub mod kind;
pub mod pred;
mod impls;

use std::cmp::Ordering;
use std::fmt::{Debug, Display, Formatter};
use std::marker::PhantomData;
use std::ops::{Deref, DerefMut};

use crate::{
	linear::*,
	kind::*,
	pred::*,
};

use self::time::Time;

pub use chime_flux_proc_macro::flux;

/// Context for a `bevy_time::Time`.
#[derive(Default)]
pub struct Chime;

/// A discrete interface for a value that changes over time.
pub trait Moment {
	type Flux: Flux<Moment = Self>;
	
	/// Constructs the entirety of a [`Flux`] from a single moment.
	fn to_flux(self, time: Time) -> Self::Flux;
	
	fn to_flux_value(self, time: Time) -> FluxValue<Self::Flux>
	where
		Self: Sized
	{
		FluxValue::new(self.to_flux(time), time)
	}
}

#[allow(type_alias_bounds)]
pub type FluxOf<T: Moment> = FluxValue<<T as Moment>::Flux>;

impl<A: Flux> FluxValue<A> {
	// !!! Deriving PartialEq, Eq should count `f(t) = 1 + 2t` and
	// `g(t) = 3 + 2(t-basis_time)` as the same Flux if `basis_time = 1`.
	
	pub fn into_inner_flux(self) -> A {
		self.inner
	}
	
	/// An evaluation of this flux at some point in time.
	pub fn basis(&self) -> <<A::Kind as FluxKind>::Basis as Basis>::Inner {
		self.inner.basis()
	}
	
	/// The time of [`Flux::basis`].
	pub fn basis_time(&self) -> Time {
		self.time
	}
	
	/// A moment in the timeline.
	pub fn to_moment(self, time: Time) -> A::Moment {
		self.inner.to_moment(self.time, time)
	}
	
	/// Sets a moment in the timeline (affects all moments).
	pub fn set_moment(&mut self, time: Time, moment: A::Moment)
	where
		Self: Sized
	{
		*self = moment.to_flux_value(time);
	}
	
	/// A reference to a moment in the timeline.
	pub fn at(&self, time: Time) -> MomentRef<A::Moment>
	where
		Self: Clone
	{
		let moment = self.clone().to_moment(time);
		MomentRef {
			moment,
			borrow: PhantomData,
		}
	}
	
	/// A unique, mutable reference to a moment in the timeline.
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
	pub fn at_mut(&mut self, time: Time) -> MomentMut<A::Moment>
	where
		Self: Clone
	{
		let moment = Some(self.clone().to_moment(time));
		MomentMut {
			moment,
			time,
			borrow: self,
		}
	}
	
	/// A point in the timeline.
	/// 
	/// `self.eval(self.basis_time()) == self.basis()`
	pub fn eval(&self, time: Time) -> <<A::Kind as FluxKind>::Basis as Basis>::Inner {
		let basis_time = self.basis_time();
		if time == basis_time {
			return self.basis()
		}
		self.poly(basis_time).eval(time)
	}
	
	/// A polynomial description of this flux at the given time.
	pub fn poly(&self, time: Time) -> Poly<A::Kind> {
		self.inner
			.change(FluxAccum {
				poly: Poly::new(Constant::from_value(self.eval(time)), time),
				time: self.time,
			})
			.poly
	}
	
	/// Ranges when this is above/below/equal to another flux.
	pub fn when<T>(&self, order: Ordering, other: &FluxValue<T>)
		-> <Poly<A::Kind> as When<T::Kind>>::Pred
	where
		T: Flux,
		Poly<A::Kind>: When<T::Kind>
	{
		let time = self.basis_time();
		self.poly(time).when(order, other.poly(time))
	}
	
	/// Times when this is equal to another flux.
	pub fn when_eq<T>(&self, other: &FluxValue<T>)
		-> <Poly<A::Kind> as WhenEq<T::Kind>>::Pred
	where
		T: Flux,
		Poly<A::Kind>: WhenEq<T::Kind>
	{
		let time = self.basis_time();
		self.poly(time).when_eq(other.poly(time))
	}
	
	/// Ranges when this is above/below/equal to a constant.
	pub fn when_constant<T>(&self, order: Ordering, other: T)
		-> <Poly<A::Kind> as When<Constant<KindLinear<A::Kind>>>>::Pred
	where
		T: LinearIso<KindLinear<A::Kind>>,
		Poly<A::Kind>: When<Constant<KindLinear<A::Kind>>>
	{
		self.when(order, &FluxValue::new(Constant::from(T::into_linear(other)), Time::ZERO))
	}
	
	/// Times when this is equal to a constant.
	pub fn when_eq_constant<T>(&self, other: T)
		-> <Poly<A::Kind> as WhenEq<Constant<KindLinear<A::Kind>>>>::Pred
	where
		T: LinearIso<KindLinear<A::Kind>>,
		Poly<A::Kind>: WhenEq<Constant<KindLinear<A::Kind>>>
	{
		self.when_eq(&FluxValue::new(Constant::from(T::into_linear(other)), Time::ZERO))
	}
}

/// Immutable moment-in-time interface for [`Flux::at`].
pub struct MomentRef<'b, M: Moment> {
	moment: M,
	borrow: PhantomData<&'b M::Flux>,
}

impl<M: Moment> Deref for MomentRef<'_, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		&self.moment
	}
}

impl<M: Moment> Debug for MomentRef<'_, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<M: Moment> Display for MomentRef<'_, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

/// Mutable moment-in-time interface for [`Flux::at_mut`].
pub struct MomentMut<'b, M: Moment> {
	moment: Option<M>,
	time: Time,
	borrow: &'b mut FluxValue<M::Flux>,
}

impl<M: Moment> Drop for MomentMut<'_, M> {
	fn drop(&mut self) {
		if let Some(moment) = std::mem::take(&mut self.moment) {
			self.borrow.set_moment(self.time, moment);
		}
	}
}

impl<M: Moment> Deref for MomentMut<'_, M> {
	type Target = M;
	fn deref(&self) -> &Self::Target {
		if let Some(moment) = self.moment.as_ref() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<M: Moment> DerefMut for MomentMut<'_, M> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		if let Some(moment) = self.moment.as_mut() {
			moment
		} else {
			unreachable!()
		}
	}
}

impl<M: Moment> Debug for MomentMut<'_, M>
where
	M: Debug
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Debug>::fmt(self, f)
	}
}

impl<M: Moment> Display for MomentMut<'_, M>
where
	M: Display
{
	fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
		<M as Display>::fmt(self, f)
	}
}

#[cfg(feature = "bevy")]
mod bevy_moment {
	use bevy_ecs::archetype::Archetype;
	use bevy_ecs::component::{Component, ComponentId, Components, Tick};
	use bevy_ecs::entity::Entity;
	use bevy_ecs::query::{FilteredAccess, QueryData, ReadOnlyQueryData, WorldQuery};
	use bevy_ecs::storage::{Table, TableRow};
	use bevy_ecs::system::{Res, ResMut, Resource, SystemMeta, SystemParam};
	use bevy_ecs::world::{World, unsafe_world_cell::UnsafeWorldCell};
	use super::*;
	
	type ChimeTime = bevy_time::Time<Chime>;
	
	type Ref<'b, M> = &'b FluxValue<<M as Moment>::Flux>;
	type Mut<'b, M> = &'b mut FluxValue<<M as Moment>::Flux>;
	
	/// SAFETY: `Self` is the same as `Self::ReadOnly`.
	unsafe impl<'b, M> QueryData for MomentRef<'b, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Component + Clone,
	{
		type ReadOnly = Self;
	}
	
	/// SAFETY: access is read only.
	unsafe impl<'b, M> ReadOnlyQueryData for MomentRef<'b, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Component + Clone,
	{}
	
	/// SAFETY: access of `MomentRef<T>` is a subset of `MomentMut<T>`.
	unsafe impl<'b, M> QueryData for MomentMut<'b, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Component + Clone,
	{
		type ReadOnly = MomentRef<'b, M>;
	}
	
	/// SAFETY:
	/// Most safety guarantees are inherited from `impl QueryData for &T`.
	/// This additionally requests a `Time` resource and takes care of adding
	/// read access for that resource in `update_component_access`.
	unsafe impl<'b, M> WorldQuery for MomentRef<'b, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Component + Clone,
	{
		type Item<'a> = MomentRef<'a, M>;
		type Fetch<'a> = (Time, <Ref<'b, M> as WorldQuery>::Fetch<'a>);
		type State = (ComponentId, <Ref<'b, M> as WorldQuery>::State);
		fn shrink<'wlong: 'wshort, 'wshort>(item: Self::Item<'wlong>) -> Self::Item<'wshort> {
			item
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
				.at(*time)
		}
		fn update_component_access((time_id, state): &Self::State, access: &mut FilteredAccess<ComponentId>) {
	        assert!(
	            !access.access().has_write(*time_id),
	            "&{} conflicts with a previous access in this query. Shared access cannot coincide with exclusive access.",
	                std::any::type_name::<ChimeTime>(),
	        );
	        access.access_mut().add_read(*time_id);
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
		M: Moment,
		FluxValue<M::Flux>: Component + Clone,
	{
		type Item<'a> = MomentMut<'a, M>;
		type Fetch<'a> = (Time, <Mut<'b, M> as WorldQuery>::Fetch<'a>);
		type State = (ComponentId, <Mut<'b, M> as WorldQuery>::State);
		fn shrink<'wlong: 'wshort, 'wshort>(item: Self::Item<'wlong>) -> Self::Item<'wshort> {
			item
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
				.at_mut(*time)
		}
		fn update_component_access((time_id, state): &Self::State, access: &mut FilteredAccess<ComponentId>) {
	        assert!(
	            !access.access().has_write(*time_id),
	            "&{} conflicts with a previous access in this query. Shared access cannot coincide with exclusive access.",
	                std::any::type_name::<ChimeTime>(),
	        );
	        access.access_mut().add_read(*time_id);
			<Mut<'b, M> as WorldQuery>::update_component_access(state, access)
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
	pub struct ResMoment<'w, M: Moment> {
		inner: MomentRef<'w, M>,
	}
	
	impl<'w, M: Moment> Deref for ResMoment<'w, M> {
		type Target = MomentRef<'w, M>;
		fn deref(&self) -> &Self::Target {
			&self.inner
		}
	}
	
	/// `SystemParam` for fetching a mutable resource at the current `Time<Chime>`.
	pub struct ResMomentMut<'w, M: Moment> {
		inner: MomentMut<'w, M>,
	}
	
	impl<'w, M: Moment> Deref for ResMomentMut<'w, M> {
		type Target = MomentMut<'w, M>;
		fn deref(&self) -> &Self::Target {
			&self.inner
		}
	}
	
	impl<'w, M: Moment> DerefMut for ResMomentMut<'w, M> {
		fn deref_mut(&mut self) -> &mut Self::Target {
			&mut self.inner
		}
	}
	
	/// SAFETY:
	/// Safety guarantees are inherited from `impl SystemParam for Res<T>`.
	unsafe impl<'w, M> SystemParam for ResMoment<'w, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Resource + Clone,
	{
		type State = (<Res<'w, ChimeTime> as SystemParam>::State, <Res<'w, FluxValue<M::Flux>> as SystemParam>::State);
		type Item<'world, 'state> = ResMoment<'world, M>;
		fn init_state(world: &mut World, system_meta: &mut SystemMeta) -> Self::State {
			let time_state = <Res<'w, ChimeTime> as SystemParam>::init_state(world, system_meta);
			let res_state = <Res<'w, FluxValue<M::Flux>> as SystemParam>::init_state(world, system_meta);
			(time_state, res_state)
		}
		unsafe fn get_param<'world, 'state>((time_state, res_state): &'state mut Self::State, system_meta: &SystemMeta, world: UnsafeWorldCell<'world>, change_tick: Tick) -> Self::Item<'world, 'state> {
			let time = <Res<'w, ChimeTime> as SystemParam>::get_param(time_state, system_meta, world, change_tick)
				.elapsed();
			let inner = <Res<'w, FluxValue<M::Flux>> as SystemParam>::get_param(res_state, system_meta, world, change_tick)
				.into_inner()
				.at(time);
			ResMoment { inner }
		}
	}
	
	/// SAFETY:
	/// Safety guarantees are inherited from `impl SystemParam for Res<T>` and
	/// `impl SystemParam for ResMut<T>`.
	unsafe impl<'w, M> SystemParam for ResMomentMut<'w, M>
	where
		M: Moment,
		FluxValue<M::Flux>: Resource + Clone,
	{
		type State = (<Res<'w, ChimeTime> as SystemParam>::State, <ResMut<'w, FluxValue<M::Flux>> as SystemParam>::State);
		type Item<'world, 'state> = ResMomentMut<'world, M>;
		fn init_state(world: &mut World, system_meta: &mut SystemMeta) -> Self::State {
			let time_state = <Res<'w, ChimeTime> as SystemParam>::init_state(world, system_meta);
			let res_state = <ResMut<'w, FluxValue<M::Flux>> as SystemParam>::init_state(world, system_meta);
			(time_state, res_state)
		}
		unsafe fn get_param<'world, 'state>((time_state, res_state): &'state mut Self::State, system_meta: &SystemMeta, world: UnsafeWorldCell<'world>, change_tick: Tick) -> Self::Item<'world, 'state> {
			let time = <Res<'w, ChimeTime> as SystemParam>::get_param(time_state, system_meta, world, change_tick)
				.elapsed();
			let inner = <ResMut<'w, FluxValue<M::Flux>> as SystemParam>::get_param(res_state, system_meta, world, change_tick)
				.into_inner()
				.at_mut(time);
			ResMomentMut { inner }
		}
	}
}

#[cfg(feature = "bevy")]
pub use bevy_moment::{ResMoment, ResMomentMut};

/// Multidimensional change over time.
pub trait FluxVector<const SIZE: usize> {
	type Kind: Vector<SIZE, Output: FluxKind> + 'static;
	
	fn index_poly(&self, index: usize, time: Time)
		-> Poly<<Self::Kind as Vector<SIZE>>::Output>;
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<T, D>(&self, other: &FluxValue<T>, order: Ordering, dis: &FluxValue<D>)
		-> <Poly<Self::Kind> as WhenDis<SIZE, T::Kind, D::Kind>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Poly<Self::Kind>: WhenDis<SIZE, T::Kind, D::Kind>,
	;
	
	/// Ranges when the distance to another vector is equal to X.
	fn when_dis_eq<T, D>(&self, other: &FluxValue<T>, dis: &FluxValue<D>)
		-> <Poly<Self::Kind> as WhenDisEq<SIZE, T::Kind, D::Kind>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Poly<Self::Kind>: WhenDisEq<SIZE, T::Kind, D::Kind>,
	;
	
	/// Ranges when the distance to another vector is above/below/equal to a constant.
	fn when_dis_constant<T, D>(&self, other: &FluxValue<T>, order: Ordering, dis: D)
		-> <Poly<Self::Kind> as WhenDis<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<Self::Kind>: WhenDis<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	;
	
	/// Ranges when the distance to another vector is equal to a constant.
	fn when_dis_eq_constant<T, D>(&self, other: &FluxValue<T>, dis: D)
		-> <Poly<Self::Kind> as WhenDisEq<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<Self::Kind>: WhenDisEq<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	;
	
	/// Ranges when a component is above/below/equal to another flux.
	fn when_index<T>(&self, index: usize, order: Ordering, other: &FluxValue<T>)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as When<T::Kind>>::Pred
	where
		T: Flux,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: When<T::Kind>
	;
	
	/// Times when a component is equal to another flux.
	fn when_index_eq<T>(&self, index: usize, other: &FluxValue<T>)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<T::Kind>>::Pred
	where
		T: Flux,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<T::Kind>
	;
	
	/// Ranges when a component is above/below/equal to a constant.
	fn when_index_constant<T>(&self, index: usize, order: Ordering, other: T)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	;
	
	/// Times when a component is equal to a constant.
	fn when_index_eq_constant<T>(&self, index: usize, other: T)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	;
}

impl<A, const SIZE: usize> FluxVector<SIZE> for FluxValue<A>
where
	A: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
{
	type Kind = A::Kind;
	
	fn index_poly(&self, index: usize, time: Time) -> Poly<<Self::Kind as Vector<SIZE>>::Output> {
		self.poly(time).index(index)
	}
	
	/// Ranges when the distance to another vector is above/below/equal to X.
	fn when_dis<T, D>(&self, other: &FluxValue<T>, order: Ordering, dis: &FluxValue<D>)
		-> <Poly<Self::Kind> as WhenDis<SIZE, T::Kind, D::Kind>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Poly<Self::Kind>: WhenDis<SIZE, T::Kind, D::Kind>,
	{
		let time = self.basis_time();
		self.poly(time)
			.when_dis(other.poly(time), order, dis.poly(time))
	}
	
	/// Ranges when the distance to another vector is equal to X.
	fn when_dis_eq<T, D>(&self, other: &FluxValue<T>, dis: &FluxValue<D>)
		-> <Poly<Self::Kind> as WhenDisEq<SIZE, T::Kind, D::Kind>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: Flux,
		Poly<Self::Kind>: WhenDisEq<SIZE, T::Kind, D::Kind>,
	{
		let time = self.basis_time();
		self.poly(time)
			.when_dis_eq(other.poly(time), dis.poly(time))
	}
	
	/// Ranges when the distance to another vector is above/below/equal to a constant.
	fn when_dis_constant<T, D>(&self, other: &FluxValue<T>, order: Ordering, dis: D)
		-> <Poly<Self::Kind> as WhenDis<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<Self::Kind>: WhenDis<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	{
		self.when_dis(other, order, &FluxValue::new(Constant::from(D::into_linear(dis)), Time::ZERO))
	}
	
	/// Ranges when the distance to another vector is equal to a constant.
	fn when_dis_eq_constant<T, D>(&self, other: &FluxValue<T>, dis: D)
		-> <Poly<Self::Kind> as WhenDisEq<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: Flux<Kind: Vector<SIZE, Output: FluxKind>>,
		D: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<Self::Kind>: WhenDisEq<SIZE, T::Kind, Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>,
	{
		self.when_dis_eq(other, &FluxValue::new(Constant::from(D::into_linear(dis)), Time::ZERO))
	}
	
	/// Ranges when a component is above/below/equal to another flux.
	fn when_index<T>(&self, index: usize, order: Ordering, other: &FluxValue<T>)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as When<T::Kind>>::Pred
	where
		T: Flux,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: When<T::Kind>
	{
		let time = self.basis_time();
		self.index_poly(index, time)
			.when(order, other.poly(time))
	}
	
	/// Times when a component is equal to another flux.
	fn when_index_eq<T>(&self, index: usize, other: &FluxValue<T>)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<T::Kind>>::Pred
	where
		T: Flux,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<T::Kind>
	{
		let time = self.basis_time();
		self.index_poly(index, time)
			.when_eq(other.poly(time))
	}
	
	/// Ranges when a component is above/below/equal to a constant.
	fn when_index_constant<T>(&self, index: usize, order: Ordering, other: T)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: When<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	{
		self.when_index(index, order, &FluxValue::new(Constant::from(T::into_linear(other)), Time::ZERO))
	}
	
	/// Times when a component is equal to a constant.
	fn when_index_eq_constant<T>(&self, index: usize, other: T)
		-> <Poly<<Self::Kind as Vector<SIZE>>::Output> as WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>>::Pred
	where
		T: LinearIso<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>,
		Poly<<Self::Kind as Vector<SIZE>>::Output>: WhenEq<Constant<KindLinear<<Self::Kind as Vector<SIZE>>::Output>>>
	{
		self.when_index_eq(index, &FluxValue::new(Constant::from(T::into_linear(other)), Time::ZERO))
	}
	
	// !!!
	// - to rotate line by a fixed angle, multiply axial polynomials by Scalar?
	// - to fit a line segment, filter predicted times.
	// - for non-fixed angle line segments, use a double distance check?
	// - rotating point-line may be handleable iteratively, find the bounds in
	//   which the roots may be and iterate through it.
}

/// Used to construct a [`Change`] for convenient change-over-time operations.
/// 
/// `1 + 2.per(time_unit::SEC)` 
pub trait Per: Sized {
	fn per(self, unit: Time) -> Change<Self> {
		Change {
			rate: self,
			unit
		}
	}
}

impl<T> Per for T {}

/// A description of a change over time for use with arithmetic operators.
#[derive(Copy, Clone, Debug, Default, Ord, PartialOrd, Eq, PartialEq)]
pub struct Change<T> {
	pub rate: T,
	pub unit: Time,
}

impl<T> Change<T> {
	pub fn as_ref(&self) -> Change<&T> {
		Change {
			rate: &self.rate,
			unit: self.unit,
		}
	}
}

impl<T: Moment> Moment for Change<T> {
	type Flux = Change<T::Flux>;
	fn to_flux(self, time: Time) -> Self::Flux {
		Change {
			rate: self.rate.to_flux(time),
			unit: self.unit,
		}
	}
}

impl<T: Flux> Flux for Change<T> {
	type Moment = Change<T::Moment>;
	type Kind = T::Kind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
		self.rate.basis()
	}
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
		self.rate.change(accum)
	}
	fn to_moment(self, basis_time: Time, time: Time) -> Self::Moment {
		Change {
			rate: self.rate.to_moment(basis_time, time),
			unit: self.unit,
		}
	}
}

/// Wrapper for partial [`Flux`] types created by the [`flux`] macro.
#[derive(Copy, Clone, Debug, Default)]
#[cfg_attr(feature = "bevy", derive(
	bevy_ecs::component::Component,
	bevy_ecs::system::Resource,
))]
pub struct FluxValue<T> {
	inner: T,
	time: Time,
}

impl<T> FluxValue<T> {
	pub fn new(inner: T, time: Time) -> Self {
		Self { inner, time }
	}
	
	pub fn time(&self) -> Time {
		self.time
	}
	
	pub fn map<U>(&self, f: impl Fn(&T) -> &U) ->  FluxValue<&'_ U>
	where
		U: Flux,
	{
		FluxValue {
			inner: f(&self.inner),
			time: self.time,
		}
	}
}

impl<T> Deref for FluxValue<T> {
	type Target = T;
	fn deref(&self) -> &Self::Target {
		&self.inner
	}
}

impl<T> DerefMut for FluxValue<T> {
	fn deref_mut(&mut self) -> &mut Self::Target {
		&mut self.inner
	}
}

/// The continuous interface for a value that changes over time.
pub trait Flux {
	type Moment: Moment<Flux = Self>;
	type Kind: FluxKind;
	fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner;
	fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind>;
	fn to_moment(self, basis_time: Time, time: Time) -> Self::Moment;
}

/// No change over time.
/// 
/// Equivalent "constant" flux kinds should implement both `Into<Constant<T>>`
/// and `From<Constant<T>>` (e.g. `Sum<T, 0>`).
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct Constant<T>(pub T);

/// ... [`<Constant as IntoIterator>::IntoIter`]
pub struct ConstantIter<T>(T);

mod _constant_impls {
	use std::ops::{Deref, DerefMut, Mul};
	use crate::{Flux, Moment};
	use crate::kind::{EmptyFluxAccum, FluxAccum, FluxKind};
	use crate::linear::{Basis, Linear, Scalar, Vector};
	use crate::time::Time;
	use super::{Constant, ConstantIter};
	
	impl<T> Deref for Constant<T> {
		type Target = T;
		fn deref(&self) -> &Self::Target {
			&self.0
		}
	}
	
	impl<T> DerefMut for Constant<T> {
		fn deref_mut(&mut self) -> &mut Self::Target {
			&mut self.0
		}
	}
	
	impl<T: Basis> From<T> for Constant<T> {
		fn from(value: T) -> Self {
			Constant(value)
		}
	}
	
	impl<T: Basis> Moment for Constant<T> {
		type Flux = Self;
		fn to_flux(self, _time: Time) -> Self::Flux {
			self
		}
	}
	
	impl<T: Basis> Flux for Constant<T> {
		type Moment = Self;
		type Kind = Self;
		fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
			self.0.clone().into_inner()
		}
		fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
			accum
		}
		fn to_moment(self, _base_time: Time, _time: Time) -> Self::Moment {
			self
		}
	}
	
	impl<T: Basis> FluxKind for Constant<T> {
		type Basis = T;
		const DEGREE: usize = 0;
		fn from_value(value: <Self::Basis as Basis>::Inner) -> Self {
			Constant(T::from_inner(value))
		}
		fn add_value(mut self, value: <Self::Basis as Basis>::Inner) -> Self {
			self.0 = T::from_inner(self.0.into_inner().add(value));
			self
		}
		fn deriv(self) -> Self {
			Self::zero()
		}
		fn eval(&self, _time: Scalar) -> <Self::Basis as Basis>::Inner {
			self.0.clone().into_inner()
		}
	}
	
	impl<T, const SIZE: usize> Vector<SIZE> for Constant<T>
	where
		T: Vector<SIZE>
	{
		type Output = Constant<T::Output>;
		fn index(&self, index: usize) -> Self::Output {
			Constant(self.0.index(index))
		}
	}
	
	impl<T: Basis> Mul<Scalar> for Constant<T> {
		type Output = Self;
		fn mul(mut self, rhs: Scalar) -> Self::Output {
			self.0 = T::from_inner(self.0.into_inner().mul_scalar(rhs));
			self
		}
	}
	
	impl<T> Mul<Constant<T>> for Constant<T>
	where
		T: Mul<Output = T>
	{
		type Output = Self;
		fn mul(mut self, rhs: Self) -> Self {
			self.0 = self.0 * rhs.0;
			self
		}
	}
	
	impl<A, B> std::ops::Add<B> for Constant<A>
	where
		A: Basis,
		B: FluxKind<Basis = A>,
	{
		type Output = B;
		fn add(self, rhs: B) -> Self::Output {
			rhs.add_value(self.0.into_inner())
		}
	}
	
	impl<T> IntoIterator for Constant<T>
	where
		T: IntoIterator,
	{
		type Item = Constant<T::Item>;
		type IntoIter = ConstantIter<T::IntoIter>;
		fn into_iter(self) -> Self::IntoIter {
			ConstantIter(self.0.into_iter())
		}
	}
	
	impl<T> Iterator for ConstantIter<T>
	where
		T: Iterator,
	{
		type Item = Constant<T::Item>;
		fn next(&mut self) -> Option<Self::Item> {
			self.0.next().map(Constant)
		}
		fn size_hint(&self) -> (usize, Option<usize>) {
			self.0.size_hint()
		}
	}
}

#[cfg(test)]
mod tests {
	use super::*;
	use super::time::{SEC, TimeRanges};
	use crate::sum::Sum;
	
	// #[flux(
	// 	kind = Sum<f64, 4>,
	// 	value = value,
	// 	change = |c| c + spd.per(SEC) + misc.per(SEC),
	// 	crate = crate,
	// )]
	#[derive(Clone, Debug, Default)]
	struct Pos {
		value: f64,
		spd: Spd,
		misc: Vec<Spd>,
	}
	
	impl Moment for Pos {
		type Flux = Self;
		fn to_flux(self, _time: Time) -> Self::Flux {
			self
		}
	}
	
	impl Flux for Pos {
		type Moment = Self;
		type Kind = Sum<f64, 4>;
		fn basis(&self) -> <<Self::Kind as FluxKind>::Basis as Basis>::Inner {
			self.value
		}
		fn change(&self, accum: EmptyFluxAccum<Self::Kind>) -> FluxAccum<Self::Kind> {
			let mut accum = accum + (&self.spd).per(SEC);
			for spd in &self.misc {
				accum = accum + spd.per(SEC);
			}
			accum
		}
		fn to_moment(self, basis_time: Time, time: Time) -> Self::Moment {
			Self {
				value: FluxValue::new(&self, basis_time).eval(time),
				spd: self.spd.to_moment(basis_time, time),
				misc: self.misc.to_moment(basis_time, time),
			}
		}
	}
	
	#[flux(
		kind = Sum<f64, 3>,
		value = value,
		change = |c| c - fric.per(SEC) + accel.per(SEC),
		crate = crate,
	)]
	#[derive(Clone, Debug, Default)]
	struct Spd {
		value: f64,
		fric: Fric,
		accel: Accel,
	}
	
	#[flux(
		kind = Sum<f64, 0>,
		value = value,
		crate = crate,
	)]
	#[derive(Clone, Debug, Default)]
	struct Fric {
		value: f64,
	}
	
	#[flux(
		kind = Sum<f64, 2>,
		value = value,
		change = |c| c + jerk.per(SEC),
		crate = crate,
	)]
	#[derive(Clone, Debug, Default)]
	struct Accel {
		value: f64,
		jerk: Jerk,
	}
	
	#[flux(
		kind = Sum<f64, 1>,
		value = value,
		change = |c| c + snap.per(SEC),
		crate = crate,
	)]
	#[derive(Clone, Debug, Default)]
	struct Jerk {
		value: f64,
		snap: Snap,
	}
	
	#[flux(
		kind = Constant<f64>,
		value = value,
		crate = crate,
	)]
	#[derive(Clone, Debug, Default)]
	struct Snap {
		value: f64,
	}
	
	impl Deref for Pos {
		type Target = f64;
		fn deref(&self) -> &Self::Target {
			&self.value
		}
	}
	
	fn position() -> FluxOf<<Pos as Moment>::Flux> {
		let pos = Pos {
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
		};
		let mut pos = pos.to_flux_value(Time::ZERO);
		*pos = pos.clone().to_moment(10*SEC);
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
			let dif_poly = poly - cmp_poly;
			for (degree, coeff) in dif_poly.kind.into_iter().enumerate() {
				assert_eq!((coeff * 2_f64.powi((1 + degree) as i32)).round(), 0.,
					"poly: {:?}, cmp_poly: {:?}", poly, cmp_poly);
			}
		}};
	}
	
	#[test]
	fn value() {
		let mut pos = position();
		assert_eq!(pos.basis_time(), 10*SEC);
		
		 // Values:
		assert_eq!(pos.at(0*SEC).round(), 32.);
		assert_eq!(pos.at(10*SEC).round(), -63.);
		assert_eq!(pos.at(20*SEC).round(), -113.);
		assert_eq!(pos.at(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
		
		 // Update:
		assert_eq!(pos.at_mut(20*SEC).round(), -113.);
		assert_eq!(pos.at(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
		assert_eq!(pos.at_mut(100*SEC).round(), 8339.);
		assert_eq!(pos.at(200*SEC).round(), -209778.);
	}
	
	#[test]
	fn poly() {
		let mut pos = position();
		assert_poly!(
			pos.poly(10*SEC),
			Poly::new(Sum::new(-63.15, [
				-11.9775,
				0.270416666666666,
				0.0475,
				-0.000416666666666,
			]), 10*SEC)
		);
		for _ in 0..2 {
			pos.at_mut(20*SEC);
			assert_poly!(
				pos.poly(pos.basis_time()),
				Poly::new(Sum::new(-112.55, [
					 6.0141666666666666666,
					 1.4454166666666666666,
					 0.0308333333333333333,
					-0.0004166666666666666,
				]), 20*SEC)
			);
		}
		pos.at_mut(0*SEC);
		assert_poly!(
			pos.poly(pos.basis_time()),
			Poly::new(Sum::new(32., [
				-1.4691666666666666666,
				-1.4045833333333333333,
				 0.0641666666666666666,
				-0.0004166666666666666,
			]), 0*SEC)
		);
		
		assert_poly!(
			Poly::new(Sum::new(-0.8_f64, [-2.7, 3.4, 2.8, 0.3]), 4000*time::MILLISEC).to_time(320*time::MILLISEC),
			Poly::new(Sum::new(-29.341750272_f64, [26.2289216, -3.13568, -1.616, 0.3]), 320*time::MILLISEC)
		);
	}
	
	#[test]
	fn when() {
		let pos = position();
		let acc = Accel {
			value: 0.3,
			jerk: Jerk {
				value: 0.4,
				snap: Snap { value: -0.01 },
			},
		}.to_flux_value(Time::ZERO);
		
		assert_time_ranges!(pos.when(Ordering::Greater, &acc), [
			(Time::ZERO, Time::from_secs_f64(4.56)),
			(Time::from_secs_f64(26.912), Time::from_secs_f64(127.394))
		]);
		assert_times!(pos.when_eq(&acc), [
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
		assert_time_ranges!(a_pos.when(Ordering::Greater, &b_pos), []);
		assert_time_ranges!(a_pos.when(Ordering::Equal, &b_pos), [(Time::ZERO, Time::MAX)]);
		assert_times!(a_pos.when_eq(&b_pos), []);
		a_pos.at_mut(20*SEC);
		
		 // Apply Changes:
		b_pos.at_mut(0*SEC).misc.push(Spd {
			value: 2.5,
			..Default::default()
		});
		b_pos.at_mut(0*SEC).misc.push(Spd {
			value: 12.,
			fric: Fric { value: 0.5 },
			..Default::default()
		});
		b_pos.at_mut(10*SEC).value -= 100.0;
		
		 // Check After:
		assert_eq!(
			Vec::from_iter(a_pos.when(Ordering::Greater, &b_pos)
				.into_ranges(Time::ZERO)
				.inclusive()),
			[
				(0*SEC, 8*SEC - time::NANOSEC),
				(50*SEC + time::NANOSEC, Time::MAX)
			]
		);
		assert_times!(a_pos.when_eq(&b_pos), [8*SEC, 50*SEC]);
	}
	
	#[test]
	fn distance() {
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = Sum<Iso<f64, i64>, 2>,
			value = value,
			change = |c| c + spd.per(time::MINUTE),
			crate = crate,
		)]
		#[derive(Debug)]
		struct Pos {
			value: Iso<f64, i64>,
			spd: Spd,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = Sum<Iso<f64, i64>, 1>,
			value = value,
			change = |c| c + acc.per(SEC),
			crate = crate,
		)]
		#[derive(Debug)]
		struct Spd {
			value: Iso<f64, i64>,
			acc: Acc,
		}
		
		#[derive(PartialOrd, PartialEq, Copy, Clone)]
		#[flux(
			kind = Constant<Iso<f64, i64>>,
			value = value,
			crate = crate,
		)]
		#[derive(Debug)]
		struct Acc {
			value: Iso<f64, i64>,
		}
		
		let a_pos = [
			Pos { value: Iso::new(3), spd: Spd { value: Iso::new(300), acc: Acc { value: Iso::new(-240) } } },
			Pos { value: Iso::new(-4), spd: Spd { value: Iso::new(120), acc: Acc { value: Iso::new(1080) } } }
		].to_flux_value(Time::ZERO);
		let b_pos = [
			Pos { value: Iso::new(8), spd: Spd { value: Iso::new(330), acc: Acc { value: Iso::new(-300) } } },
			Pos { value: Iso::new(4), spd: Spd { value: Iso::new(600), acc: Acc { value: Iso::new(720) } } }
		].to_flux_value(Time::ZERO);
		
		let dis = Spd { value: Iso::new(10), acc: Acc { value: Iso::new(0) } }
			.to_flux_value(Time::ZERO);
		
		assert_time_ranges!(
			a_pos.when_dis(&b_pos, Ordering::Less, &dis),
			// https://www.desmos.com/calculator/spxyoloyx9
			[
				(Time::ZERO, Time::from_secs_f64(0.0823337)),
				(Time::from_secs_f64(2.246506069), Time::from_secs_f64(4.116193987)),
			]
		);
		
		let b_pos = b_pos.to_moment(Time::ZERO).to_flux_value(SEC);
		assert_time_ranges!(
			a_pos.when_dis_eq_constant(&b_pos, 2),
			[
				(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
				(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191)),
			]
		);
		assert_time_ranges!(
			a_pos.when_dis_constant(&b_pos, Ordering::Equal, 2),
			[
				(Time::from_secs_f64(0.229597034), Time::from_secs_f64(0.414068993)),
				(Time::from_secs_f64(0.689701729), Time::from_secs_f64(0.84545191))
			]
		);
		
		// https://www.desmos.com/calculator/23ic1ikyt3
	}
}