//! Physical parameters and constants for Chrono-Gauge Lattice Theory.
//!
//! This module provides:
//! - Lattice parameters for chrono-gauge field construction
//! - β coupling interpretation (chrono-gravitational stiffness)
use crate::{
    EARTH_GM, EARTH_RADIUS, NEWTONIAN_CONSTANT_OF_GRAVITATION as GRAVITATIONAL_CONSTANT,
    SPEED_OF_LIGHT,
};
use deep_causality_num::RealField;

/// Speed of light in m/s.
/// Earth's geocentric gravitational constant (GM) in m³/s².
/// Earth's mean radius in meters.
/// Parameters for constructing a chrono-gauge lattice.
///
/// These parameters control the physical interpretation and numerical
/// properties of the latrice gauge field.
#[derive(Debug, Clone)]
pub struct ChronoLatticeParams<R> {
    /// Lattice spacing in meters.
    ///
    /// Typical values:
    /// - GNSS constellation: 1000m (1 km) to 10000m (10 km)
    /// - Ground-based clocks: 1m to 100m
    pub spacing_m: R,

    /// Inverse coupling parameter β = c² / (G · a²).
    ///
    /// # Physical Interpretation
    ///
    /// - **Large β → weak coupling → near-flat spacetime**: Plaquettes near identity
    /// - **Small β → strong coupling → highly curved**: Large deviations from flatness
    ///
    /// At GNSS scales, β is very large (~10²¹), indicating spacetime is extremely stiff.
    ///
    /// # Practical Usage
    ///
    /// For numerical stability, set β = 1.0 and interpret observables as fractional
    /// deviations from flat spacetime. The absolute GM value is recovered by rescaling.
    pub beta: R,

    /// Reference time dilation (typically 1.0 for geocentric reference).
    pub tau_ref: R,
}

impl ChronoLatticeParams<f64> {
    /// Default parameters for GNSS satellite-scale lattice.
    ///
    /// Uses 1 km spacing with dimensionless β = 1.0.
    #[inline]
    pub fn gnss_default() -> Self {
        Self {
            spacing_m: 1_000.0, // 1 km spacing
            beta: 1.0,          // Dimensionless, rescale outputs
            tau_ref: 1.0,       // Geocentric reference
        }
    }

    /// Parameters for high-resolution ground-based clock networks.
    ///
    /// Uses 100m spacing for terrestrial experiments.
    #[inline]
    pub fn ground_clock_default() -> Self {
        Self {
            spacing_m: 100.0, // 100m spacing
            beta: 1.0,        // Dimensionless
            tau_ref: 1.0,     // Geocentric reference
        }
    }

    /// Compute the physical β coupling from lattice spacing.
    ///
    /// Returns β = c² / (G · a²), which is the chrono-gravitational stiffness.
    #[inline]
    pub fn physical_beta(spacing_m: f64) -> f64 {
        let c_sq = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
        let a_sq = spacing_m * spacing_m;
        c_sq / (GRAVITATIONAL_CONSTANT * a_sq)
    }

    /// Compute the gravitational potential scale φ_max = GM / (R · c²).
    ///
    /// This is the dimensionless gravitational parameter at Earth's surface.
    #[inline]
    pub fn potential_scale() -> f64 {
        EARTH_GM / (EARTH_RADIUS * SPEED_OF_LIGHT * SPEED_OF_LIGHT)
    }
}

impl<R: RealField> ChronoLatticeParams<R> {
    /// Creates new lattice parameters with custom values.
    #[inline]
    pub fn new(spacing_m: R, beta: R, tau_ref: R) -> Self {
        Self {
            spacing_m,
            beta,
            tau_ref,
        }
    }
}
