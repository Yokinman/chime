//! Change over time.

pub mod constant;
pub mod kind;
pub mod temporal;
pub mod linear;
pub mod sum;
pub mod exp;

mod flux;

pub use flux::*;