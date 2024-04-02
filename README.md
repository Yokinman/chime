## Goals

- [x] Can statically define how a type changes over time.
- [x] Values have a continuous representation based on time; change occurs implicitly, not explicitly.
- [x] Can predict when a value will match a given state with high accuracy.
- [ ] (WIP) Can define systems that schedule events based on a given predictive predicate.
- [x] Events can use their scheduled time to fetch the momentary state of values.
- [x] Events are ignorant of the rate of time. Time can be advanced at any speed without creating inconsistent behavior.
- [x] Events run in order of scheduled time, avoiding ambiguous discrete jumps that may "miss" events.
- [x] Modifying a value automatically reschedules all events that used a predicate involving the value.
- [ ] Rescheduling events is fairly speedy. The vast majority of logic is dedicated to running scheduled events and revising previous predictions in a consistent and highly parallelizable way.
- [ ] ? Values can share change over time without requiring a timer-based update event.
- [ ] ? Can move a limited distance back in time, undoing events and returning values to their past state.



## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.