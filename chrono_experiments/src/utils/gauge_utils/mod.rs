//! Gauge utilities for the E02 chrono_gauge experiment.
//!
//! This module provides utilities for:
//! - Lattice construction from satellite data
//! - Observable computation from ChronoGauge fields
//! - J2 oblateness estimation from latitude sweep

mod lattice_utils;
mod observable_utils;

pub use lattice_utils::*;
pub use observable_utils::*;
