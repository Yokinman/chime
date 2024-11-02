//! Change over time.

pub mod kind;
pub mod temporal;
pub mod linear;
pub mod exp;

mod flux;
mod impls;

pub use flux::*;