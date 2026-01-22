//! Gravity-related types for GQCD experiments.
//!
//! These types are generic over `R: RealField` to support both `f64` and
//! `DoubleFloat` precision.

use deep_causality_num::RealField;
use deep_causality_physics::SPEED_OF_LIGHT;

/// Represents a point in 4D Space-Time with associated kinematic and clock data.
///
/// This is the fundamental unit of data for GQCD. It maps a measurable "Clock Event"
/// (timestamp + bias) to a "Geometric Event" (position + velocity) ensuring that
/// $M \leftrightarrow T$ (Mass-Time Equivalence) can be calculated.
///
/// # Type Parameter
///
/// - `R`: Real field type (e.g., `f64`, `DoubleFloat`)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpaceTimeCoordinate<R: RealField> {
    /// The UTC timestamp in seconds (Unix Epoch).
    pub timestamp: u64,
    /// Satellite ID (e.g., E14)
    pub sat_id: u32,
    /// Radius from Earth's center of mass (meters) $|r|$
    pub r_m: R,
    /// Velocity magnitude relative to Earth (m/s) $|v|$
    pub v_ms: R,
    /// Raw clock bias (seconds), uncorrected for relativistic effects.
    pub clock_bias_s: R,
    /// Full 3D position vector [x, y, z] (ITRF frame)
    pub position: [R; 3],
    /// Full 3D velocity vector [vx, vy, vz] (ITRF frame)
    pub velocity: [R; 3],
    /// Rate of change of clock bias (seconds/second), i.e., frequency offset.
    /// This corresponds to $\mathcal{T} - 1$.
    pub clock_drift_rate: R,
}

impl<R: RealField + From<f64>> SpaceTimeCoordinate<R> {
    /// Helper to restore relativistic effects removed by IGS.
    /// Calculates $\Delta t_{periodic} = -2(\vec{r} \cdot \vec{v}) / c^2$
    pub fn get_total_bias(&self) -> R {
        let dot_rv = self.position[0] * self.velocity[0]
            + self.position[1] * self.velocity[1]
            + self.position[2] * self.velocity[2];
        let c_sq = R::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);
        let two = R::from(2.0);
        let rel_correction = -two * dot_rv / c_sq;
        self.clock_bias_s + rel_correction
    }
}

/// Implements the `SpaceTimeCoord` trait from CGLT for direct use with
/// `ChronoGaugeOps::source()` and other physics operations.
impl<R: RealField + From<f64>> deep_causality_physics::SpaceTimeCoord<R>
    for SpaceTimeCoordinate<R>
{
    fn clock_drift_rate(&self) -> R {
        self.clock_drift_rate
    }

    fn radius_m(&self) -> R {
        self.r_m
    }

    fn inertial_velocity_magnitude(&self) -> R {
        self.v_ms
    }

    fn z_m(&self) -> R {
        self.position[2]
    }
}

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
