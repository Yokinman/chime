//! Change over time.

pub mod kind;
pub mod temporal;
pub mod time;
pub mod linear;
pub mod pred;
pub mod exp;

mod flux;
mod impls;

pub use flux::*;