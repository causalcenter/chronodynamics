//! Type aliases for Chrono-Gauge Lattice Theory.
//!
//! This module provides convenient type aliases for the chrono-gauge lattice fields
//! used in chrono-dynamics experiments.
//!
//! # Gauge Group: U(1) × SU(2) (Electroweak)
//!
//! The chrono-gauge theory uses the electroweak gauge group:
//! - **U(1) sector**: Scalar gravitational potential (Newtonian time dilation)
//! - **SU(2) sector**: Gravitomagnetic effects (frame-dragging, vorticity)
//!
//! # Precision
//!
//! For high-precision clocks (optical lattice clocks with $10^{-18}$ to $10^{-21}$
//! fractional stability), use `DoubleFloat` for quad-precision.
//!
use crate::SpaceTimeCoordinate;
use deep_causality_num::Complex;
use deep_causality_topology::{LatticeGaugeField, SE3, SU2_U1};

const DIM_4D_SPACE_TIME: usize = 4; // 4D spacetime (3 Space + 1 Time)
const DIM_2D_SPACE_TIME: usize = 2; // 2D spacetime (1 Space + 1 Time)

/// Chrono-Gauge lattice field: U(1) × SU(2) gauge theory on 4D spacetime.
///
/// This is the primary type for chrono-dynamics computations on GNSS satellite
/// constellations and other gravitational timing experiments.
///
/// # Type Parameter
///
/// * `FloatType` - The floating-point type for calculations. Must implement
///   `RealField` from `deep_causality_num`. For high-precision clocks, use
///   `DoubleFloat`; for prototyping, `f64` is acceptable.
///
/// # Physics
///
/// Link variables $U_\mu(x)$ encode the time dilation factor between adjacent
/// lattice sites. The U(1) phase represents the scalar gravitational potential,
/// while the SU(2) matrix encodes gravitomagnetic (frame-dragging) effects.
///
///
pub type ChronoGauge<FloatType, FieldSource = Vec<SpaceTimeCoordinate<FloatType>>> =
    LatticeGaugeField<SU2_U1, DIM_4D_SPACE_TIME, Complex<FloatType>, FloatType, FieldSource>;

pub type KineticGauge<FloatType> = LatticeGaugeField<SE3, DIM_2D_SPACE_TIME, FloatType, FloatType>;
