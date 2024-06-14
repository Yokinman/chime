## Goals

- [x] Can define how any value changes over time.[^a]
- [x] All change over time is inherently continuous.[^b] 
- [x] Can accurately predict when a value will match a given state.[^c]
- [x] Predictions can be combined and used to schedule events.[^d]
- [x] Predictions and scheduled events remain in sync with the values used to create them[^e].
- [x] Events can fetch and modify values, acting as a "moment" in time.[^f]
- [x] Scheduled events run in a consistent order that avoids discrete jump-like behavior.[^g]
- [x] Making predictions and scheduling events is fairly fast.[^h]
- [ ] ? Values can share change over time without requiring a timer-based update event.
- [ ] ? Can move a limited distance back in time, undoing events and returning values to their past state.

[^a]: Rust's trait system is used to describe how types change over time. All values of a type are assumed to be capable of the same change over time. This is fast and catches errors instantly, but requires a little boilerplate for each type that changes over time.

[^b]: Every value is associated with the structure of a polynomial, and this can be used to calculate or replace each value's state at any point in time. As a consequence, every value has a continuous representation in which change occurs implicitly as time advances. It's worth noting that this generally provides less control and can make it more difficult (or impossible) to describe highly complex change.

[^c]: Polynomial root-finding is used to predict the times at which a given value will be in a given state. This is generally fast and allows for complex predictions, such as when multiple values described by different polynomials are equal to, above, below, or within a given distance of each other.

[^d]: Special systems are assigned to the Bevy world, each with a scheduler function that runs over the permutations of its parameters. For each case, predictions are collected and an event is scheduled. All predictions involve ranges of time with two endpoints, and therefore each event involves two functions - "begin" and "end".

[^e]: The scheduler function[^d] only runs over "updated" permutations, in which at least one of the items is newly added or changed. Each permuted case is linked to its scheduled event by a unique ID, which is used to replace that event when making a new prediction.

[^f]: The "begin" and "end" functions of each event are essentially Bevy systems; interfacing with the world through arbitrary parameters. Events have knowledge of their scheduled time(s), and can use that time to fetch parameters at the current "moment".

[^g]: At its simplest, events are executed in order of scheduled time and are ignorant of the rate of time. Thus, time can be advanced at any speed without inconsistent behavior. For performance reasons, sequential events that don't functionally overlap may be executed out of order.

[^h]: Most logic is dedicated to executing events in a linear order and revising predictions in a consistent and highly parallelizable way.

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