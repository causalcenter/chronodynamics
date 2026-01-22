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
    // Monte Carlo Methods
    // =========================================================================

    /// Thermalizes the gauge field with Metropolis sweeps.
    ///
    /// Performs `n_sweeps` Metropolis updates with step size `epsilon`.
    /// Returns the acceptance rate (should be 0.3-0.5 for efficiency).
    ///
    /// # Arguments
    ///
    /// * `n_sweeps` - Number of full lattice sweeps
    /// * `epsilon` - Step size for link updates
    /// * `rng` - Random number generator
    fn thermalize<Rng: deep_causality_rand::Rng>(
        &mut self,
        n_sweeps: usize,
        epsilon: R,
        rng: &mut Rng,
    ) -> Result<f64, TopologyError>;

    /// Computes an observable with jackknife error estimation.
    ///
    /// Performs `n_measurements` measurements with `skip` sweeps between each.
    /// Returns (mean, error) of the observable.
    ///
    /// # Arguments
    ///
    /// * `observable` - Function that computes the observable from the field
    /// * `n_measurements` - Number of measurements
    /// * `skip` - Number of sweeps between measurements (decorrelation)
    /// * `epsilon` - Step size for Metropolis updates
    /// * `rng` - Random number generator
    fn measure_with_error<F, Rng>(
        &mut self,
        observable: F,
        n_measurements: usize,
        skip: usize,
        epsilon: R,
        rng: &mut Rng,
    ) -> Result<(R, R), TopologyError>
    where
        F: Fn(&Self) -> Result<R, TopologyError>,
        Rng: deep_causality_rand::Rng;

    /// Automatically tunes the Metropolis step size `epsilon`.
    ///
    /// Iteratively adjusts epsilon to achieve a target acceptance rate.
    ///
    /// # Arguments
    ///
    /// * `initial_epsilon` - Starting step size
    /// * `target_acceptance` - Desired acceptance rate (e.g. 0.5)
    /// * `max_steps` - Maximum tuning iterations
    /// * `rng` - Random number generator
    fn auto_tune_metropolis<Rng: deep_causality_rand::Rng>(
        &mut self,
        initial_epsilon: R,
        target_acceptance: f64,
        max_steps: usize,
        rng: &mut Rng,
    ) -> Result<R, TopologyError>;
}
