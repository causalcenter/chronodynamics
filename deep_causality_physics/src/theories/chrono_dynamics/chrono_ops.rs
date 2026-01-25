//! Chrono-gauge specific operations trait and implementation.
//!
//! This module provides the `ChronoGaugeOps` trait, which defines operations
//! specific to chrono-gauge lattice theory, including:
//!
//! - Mass/energy observables (via Wilson action)
//! - Vorticity and frame-dragging (via SU(2) plaquettes)
//! - Thermodynamic observables (Tolman temperature)
//! - Fleet coherence for gravitational observatory
//!
//! # Continuum Limit
//!
//! All observables should converge to continuum values as lattice spacing a → 0.
//! With standard Wilson action: O(a²) discretization errors.
//! With Symanzik improvement: O(a⁴) discretization errors.

use crate::SpaceTimeCoord;
use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;

/// Operations specific to chrono-gauge theory.
///
/// This trait extends `ChronoGauge` with physics-specific observables and
/// Monte Carlo methods for chrono-gravitational computations.
///
/// # Implementation Notes
///
/// All methods return `Result<T, TopologyError>` to propagate errors from
/// underlying lattice operations.
pub trait ChronoGaugeOps<R: RealField> {
    // =========================================================================
    // Mass/Energy Observables
    // =========================================================================

    /// Computes effective mass density from Wilson action.
    ///
    /// This is the lattice analog of the Laplacian operator in DEC:
    /// ρ ∝ ∆T where T is the chronometric field.
    ///
    /// # Continuum Limit
    ///
    /// As a → 0:
    /// $$\rho_{CGLT} \to \rho_{continuum} + O(a^2)$$
    ///
    /// With Symanzik improvement: O(a⁴)
    fn mass_density_action(&self) -> Result<R, TopologyError>;

    /// Computes virial ratio from action terms.
    ///
    /// The virial ratio measures the balance between kinetic and potential
    /// energy in the gravitational system.
    fn virial_ratio(&self) -> Result<R, TopologyError>;

    // =========================================================================
    // Angular Momentum / Frame-Dragging
    // =========================================================================

    /// Computes the vorticity tensor from SU(2) spatial plaquettes.
    ///
    /// Non-zero vorticity indicates gravitomagnetic (frame-dragging) effects.
    /// This is extracted from the SU(2) sector of the U(1) × SU(2) gauge group.
    fn vorticity_tensor(&self) -> Result<[[R; 3]; 3], TopologyError>;

    /// Computes the z-component of vorticity (parallel to Earth's rotation).
    ///
    /// For GNSS satellites, this measures the frame-dragging contribution
    /// from Earth's angular momentum.
    fn vorticity_z(&self) -> Result<R, TopologyError>;

    // =========================================================================
    // Thermodynamic Observables
    // =========================================================================

    /// Computes Tolman temperature from temporal Polyakov loop.
    ///
    /// The Tolman temperature relates local temperature to gravitational potential:
    /// T_local = T_∞ / √(g₀₀)
    ///
    /// On the lattice, this is extracted from the Polyakov loop magnitude.
    fn tolman_temperature(&self) -> Result<R, TopologyError>;

    /// Computes action-phase correlation for wave validation.
    fn action_phase_correlation(&self) -> Result<R, TopologyError>;

    // =========================================================================
    // Einstein Field Equation Inversion
    // =========================================================================

    /// Derives GM using Wilson action measurements across radial positions.
    /// Uses the relation: s(r) ∝ (GM/r²)²
    /// Extracts GM from sqrt(action) * r².
    ///
    /// # Workflow
    ///
    /// 1. The field must have source data attached via `with_source()`
    /// 2. Links must be populated via `populate_links_from_source()`
    /// 3. This method extracts GM from the populated lattice observables
    ///
    /// # Returns
    ///
    /// The gravitational parameter GM in m³/s².
    fn solve_gm_from_kinectic(&self) -> Result<R, TopologyError>;

    /// Derives GM using Wilson action measurements across radial positions.
    ///
    /// This determines GM from the field curvature (Wilson action density),
    /// realizing the Kinematic Inversion where matter properties are derived
    /// from field topology.
    ///
    /// Formula: GM ∝ r² ⋅ √Action
    fn solve_gm_from_action(&self) -> Result<R, TopologyError>;

    /// Inverts the Einstein field equation to compute gravity mass GM analytically.
    ///
    /// This method solves the Einstein field equations in reverse, deriving
    /// the gravitational parameter GM from observed curvature (clock effects)
    /// and kinetic energy differences between two space-time coordinates.
    ///
    /// # Formula
    ///
    /// $$GM = \frac{c^2(\dot{\tau}_b - \dot{\tau}_a) + \frac{1}{2}(v_b^2 - v_a^2)}{1/r_a - 1/r_b}$$
    fn solve_gm_analytical<C>(&self, coord_a: &C, coord_b: &C) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>;

    /// Computes J2 oblateness coefficient analytically (without the Gauge Field)
    /// from observed satellite time data.
    fn solve_j2_analytical<C>(&self, data: &[C]) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>;
}

/// Mutable operations for populating the gauge field from source data.
pub trait ChronoGaugeMutOps<R: RealField> {
    /// Populates temporal link variables from the internal source data.
    ///
    /// This method reads clock drift rates from `self.source()`, bins them
    /// by radius, and encodes the average drift into temporal link phases:
    ///
    /// $$U_0(x) = \exp(i \cdot (1 - \dot{\tau}) \cdot N_t)$$
    ///
    /// # Requirements
    ///
    /// The field must have source data attached (S = Vec<SpaceTimeCoordinate>).
    fn populate_links_from_source(&mut self) -> Result<(), TopologyError>;
}
