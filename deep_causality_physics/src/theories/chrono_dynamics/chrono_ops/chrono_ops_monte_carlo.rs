use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;

pub trait ChronoGaugeMonteCarloOps<R: RealField> {
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
