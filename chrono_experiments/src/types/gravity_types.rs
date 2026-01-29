//! Gravity-related types for GQCD experiments.
//!
//! These types are generic over `R: RealField` to support both `f64` and
//! `DoubleFloat` precision.

use deep_causality_num::RealField;

/// Result from solving for GM using time dilation measurements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GravitySolverResult<R: RealField> {
    /// The geocentric gravitational constant GM, derived purely from clock measurements.
    /// This is the TRUE measurement. G and M are only recovered by assuming the other.
    pub gm: R,
    pub term_potential: R, // The geometric contribution (1/r - 1/r)
    pub term_kinetic: R,   // The velocity contribution (v^2 - v^2)
    pub term_time: R,      // The clock contribution (rate - rate)
}

impl<R: RealField> GravitySolverResult<R> {
    /// Creates a new result with all components.
    pub fn new(gm: R, term_potential: R, term_kinetic: R, term_time: R) -> Self {
        Self {
            gm,
            term_potential,
            term_kinetic,
            term_time,
        }
    }
}

/// Structure to hold GM result and associated altitude for correlation analysis
#[derive(Debug, Clone, Copy)]
pub struct GmDataPoint<R: RealField> {
    pub altitude: R,
    pub gm: R,
    pub radial_velocity: R,
    pub latitude: R,
    pub longitude: R,
}

impl<R: RealField> GmDataPoint<R> {
    pub fn new(altitude: R, gm: R, radial_velocity: R, latitude: R, longitude: R) -> Self {
        Self {
            altitude,
            gm,
            radial_velocity,
            latitude,
            longitude,
        }
    }
}
