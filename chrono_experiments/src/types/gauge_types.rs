//! Gauge experiment types for E02 chrono_gauge experiment.
//!
//! This module provides types for lattice gauge theory validation of the
//! chronometric field against theoretical predictions.

use deep_causality_num::{FromPrimitive, RealField};

/// Latitude bin for accumulating residuals.
#[derive(Clone)]
pub struct LatitudeBin<R: RealField> {
    pub residual_sum: R,
    pub count: usize,
}

impl<R: RealField> Default for LatitudeBin<R> {
    fn default() -> Self {
        Self {
            residual_sum: R::zero(),
            count: 0,
        }
    }
}
/// Result of J2 estimation.
#[derive(Debug, Clone, Copy)]
pub struct J2Result<R: RealField> {
    /// The derived J2 coefficient.
    pub derived_j2: R,
    /// The reference J2 value (JGM-3).
    pub reference_j2: R,
    /// The relative error in percent.
    pub error_percent: R,
    /// Number of data points used.
    pub data_points: usize,
    /// Number of valid bins used in regression.
    pub valid_bins: usize,
    /// The regression R² value (goodness of fit).
    pub r_squared: R,
}

/// Metrics computed per epoch from the ChronoGauge lattice.
#[derive(Debug, Clone, Default)]
pub struct EpochMetrics<R> {
    /// Mass density from Wilson action (LGT analog of Laplacian)
    pub mass_density: R,
    /// Frobenius norm of vorticity tensor (LGT analog of curl magnitude)
    pub curl_magnitude: R,
    /// U(1) link phase gradient magnitude (LGT analog of gradient)
    pub gradient_magnitude: R,
    /// SU(2) z-component of vorticity (Earth rotation axis alignment)
    pub vorticity_z: R,
    /// Time-velocity correlation (momentum validation)
    pub time_velocity_corr: R,
    /// Virial ratio |Kinetic| / |Potential| (energy validation)
    pub virial_ratio: R,
    /// Tolman temperature from Polyakov loop (thermodynamics validation)
    pub tolman_temp: R,
    /// Fleet coherence mean (Polyakov loop magnitude mean)
    pub coherence_mean: R,
    /// Fleet coherence variance (Polyakov loop magnitude variance)
    pub coherence_var: R,
    /// Derived Earth rotation rate from momentum optimization
    pub derived_rotation: R,
    /// Derived GM from clock rate
    pub derived_gm: R,
    /// Derived Mass from GM
    pub derived_mass: R,
    /// Derived Gravity from GM/r^2
    pub derived_gravity: R,
}

/// Result of LGT validation for one dataset.
#[derive(Debug, Clone, Default)]
pub struct GaugeValidationResult<R> {
    /// Dataset name (e.g., "gbm18710")
    pub dataset_name: String,
    /// Number of valid epochs processed
    pub num_epochs: usize,
    /// Lattice sites per epoch
    pub num_lattice_sites: usize,

    // === Core observables (LGT equivalents) ===
    /// Mean mass density from Wilson action
    pub mean_mass_density: R,
    /// Mean curl magnitude from vorticity tensor
    pub mean_curl_magnitude: R,
    /// Mean gradient magnitude from U(1) phases
    pub mean_gradient_magnitude: R,
    /// Mean vorticity z-component
    pub mean_vorticity_z: R,

    // === Physics correlations ===
    /// Time-velocity correlation (momentum)
    pub time_velocity_correlation: R,
    /// Virial ratio (energy)
    pub virial_ratio: R,
    /// Tolman consistency (thermodynamics)
    pub tolman_consistency: R,
    /// Derived Earth rotation rate
    pub derived_earth_rotation: R,
    /// Derived GM (mean)
    pub derived_gm: R,
    /// Derived Mass (mean)
    pub derived_mass: R,
    /// Derived Gravity (mean)
    pub derived_gravity: R,

    // === Fleet coherence (new in LGT) ===
    /// Fleet coherence (mean, variance)
    pub fleet_coherence: (R, R),
    /// Topological charge
    pub topological_charge: R,
}

impl<R: RealField + Default + Clone + FromPrimitive> GaugeValidationResult<R> {
    /// Creates a new validation result from epoch metrics.
    pub fn from_epochs(
        dataset_name: String,
        epochs: &[EpochMetrics<R>],
        num_lattice_sites: usize,
    ) -> Self {
        let n = epochs.len();
        if n == 0 {
            return Self {
                dataset_name,
                num_epochs: 0,
                num_lattice_sites,
                ..Default::default()
            };
        }

        let n_r = R::from_usize(n).unwrap_or_else(R::one);

        // Sum all metrics
        let mut sum = EpochMetrics::<R>::default();
        for e in epochs {
            sum.mass_density += e.mass_density;
            sum.curl_magnitude += e.curl_magnitude;
            sum.gradient_magnitude += e.gradient_magnitude;
            sum.vorticity_z += e.vorticity_z;
            sum.time_velocity_corr += e.time_velocity_corr;
            sum.virial_ratio += e.virial_ratio;
            sum.tolman_temp += e.tolman_temp;
            sum.coherence_mean += e.coherence_mean;
            sum.coherence_var += e.coherence_var;
            sum.derived_rotation += e.derived_rotation;
            sum.derived_gm += e.derived_gm;
            sum.derived_mass += e.derived_mass;
            sum.derived_gravity += e.derived_gravity;
        }

        Self {
            dataset_name,
            num_epochs: n,
            num_lattice_sites,
            mean_mass_density: sum.mass_density / n_r,
            mean_curl_magnitude: sum.curl_magnitude / n_r,
            mean_gradient_magnitude: sum.gradient_magnitude / n_r,
            mean_vorticity_z: sum.vorticity_z / n_r,
            time_velocity_correlation: sum.time_velocity_corr / n_r,
            virial_ratio: sum.virial_ratio / n_r,
            tolman_consistency: sum.tolman_temp / n_r,
            derived_earth_rotation: sum.derived_rotation / n_r,
            derived_gm: sum.derived_gm / n_r,
            derived_mass: sum.derived_mass / n_r,
            derived_gravity: sum.derived_gravity / n_r,
            fleet_coherence: (sum.coherence_mean / n_r, sum.coherence_var / n_r),
            topological_charge: R::default(),
        }
    }
}
