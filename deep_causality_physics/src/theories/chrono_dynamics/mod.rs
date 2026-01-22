//! Chrono-Dynamics module implementing Chrono-Gauge Lattice Theory (CGLT).
//!
//! This module provides:
//! - `ChronoGaugeWitness`: HKT4 witness for RiemannMap operations
//! - `ChronoVector`: Concrete vector type for curvature calculations
//! - `ChronoGaugeOps`: Trait for chrono-gauge specific operations
//!
//! # Architecture
//!
//! The chrono-gauge theory recasts gravitational time dilation as a lattice gauge
//! field with U(1) × SU(2) gauge group. Key mappings:
//!
//! | Physical Observable | Lattice Gauge Equivalent |
//! |---------------------|--------------------------|
//! | Gravitational potential | U(1) temporal link phase |
//! | Frame-dragging | SU(2) spatial plaquettes |
//! | Mass density (Laplacian) | Wilson action |
//! | Vorticity | Polyakov loop winding |
//!
//! # Safety Note
//!
//! The HKT implementations use unsafe pointer casting to work around Rust's
//! current GAT trait solver limitations. See module documentation for details.

mod chrono_hkt;
mod chrono_ops;
mod chrono_ops_impl;
mod chrono_ops_monte_carlo;
mod chrono_ops_monte_carlo_impl;
mod chrono_params;
mod chrono_vector;

pub use chrono_hkt::*;
pub use chrono_ops::*;
pub use chrono_ops_monte_carlo::*;
pub use chrono_params::*;
pub use chrono_vector::*;
