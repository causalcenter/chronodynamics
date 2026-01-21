//! Strong Gravity Binary Pulsar Types
//!
//! Data structures for analyzing binary pulsars from NANOGrav to validate
//! GQCD in strong gravitational fields (Φ/c² ~ 10⁻⁸ to 10⁻⁷).

use deep_causality_num::RealField;

/// Binary pulsar orbital parameters from .par file
#[derive(Debug, Clone)]
pub struct BinaryPulsarParams<R>
where
    R: RealField,
{
    /// Pulsar J-name (e.g., "J0437-4715")
    pub psrj: String,

    /// Binary model (DD, ELL1, BT)
    pub binary_model: String,

    /// Orbital period (days)
    pub pb_days: Option<R>,

    /// Eccentricity
    pub ecc: Option<R>,

    /// Semi-major axis (light-seconds)
    pub a1_lt_s: Option<R>,

    /// Epoch of periastron (MJD)
    pub t0_mjd: Option<R>,

    /// Longitude of periastron (degrees)
    pub om_deg: Option<R>,

    /// Pulsar mass (solar masses)
    pub m_pulsar: Option<R>,

    /// Companion mass (solar masses)
    pub m_companion: Option<R>,

    /// Gravitational redshift parameter (milliseconds)
    pub gamma_ms: Option<R>,

    /// Orbital decay (dimensionless)
    pub pbdot: Option<R>,

    /// ELL1 First Laplace-Lagrange parameter (e * sin(omega))
    pub eps1: Option<R>,

    /// ELL1 Second Laplace-Lagrange parameter (e * cos(omega))
    pub eps2: Option<R>,

    /// ELL1 Epoch of ascending node (MJD)
    pub tasc_mjd: Option<R>,
}

impl<R> BinaryPulsarParams<R>
where
    R: RealField,
{
    pub fn new(psrj: String) -> Self {
        Self {
            psrj,
            binary_model: String::new(),
            pb_days: None,
            ecc: None,
            a1_lt_s: None,
            t0_mjd: None,
            om_deg: None,
            m_pulsar: None,
            m_companion: None,
            gamma_ms: None,
            pbdot: None,
            eps1: None,
            eps2: None,
            tasc_mjd: None,
        }
    }

    pub fn is_usable(&self) -> bool {
        // Must have: binary model, orbital period, semi-major axis
        // And either (T0) OR (TASC for ELL1/ELL1H)
        !self.binary_model.is_empty()
            && self.pb_days.is_some()
            && self.a1_lt_s.is_some()
            && (self.t0_mjd.is_some() || self.tasc_mjd.is_some())
    }

    /// Check if masses are well-constrained
    /// Exclude if masses are missing (poorly constrained)
    pub fn has_well_constrained_masses(&self) -> bool {
        self.m_pulsar.is_some() && self.m_companion.is_some()
    }
}

/// Single epoch timing measurement for binary pulsar
#[derive(Debug, Clone)]
pub struct BinaryEpoch<R>
where
    R: RealField,
{
    /// Modified Julian Date
    pub mjd: R,

    /// Timing residual (nanoseconds)
    pub toa_residual_ns: R,

    /// Orbital phase (0 = periastron, 0.5 = apastron)
    pub orbital_phase: R,

    /// Pulsar-companion separation (meters)
    pub separation_m: R,

    /// Field strength Φ/c²
    pub field_strength_phi: R,
}

/// Results of analyzing a single binary system
#[derive(Debug, Clone)]
pub struct BinarySystemResult<R>
where
    R: RealField,
{
    pub pulsar_name: String,
    pub num_epochs: usize,
    pub field_strength_mean: R,
    pub field_strength_std: R,
    pub energy_balance_ratio: R,
    pub action_time_integral: R,
    pub action_orbital_computed: R,
    pub passes_criterion: bool,
    pub exclusion_reason: Option<String>,
}

impl<R> BinarySystemResult<R>
where
    R: RealField + Clone + From<f64>,
{
    /// Criterion for passing: |R - 1| < 0.1
    pub fn compute_passes(ratio: R) -> bool {
        let one = R::one();
        let threshold = R::from(0.1_f64);
        let diff = ratio - one;
        diff.abs() < threshold
    }
}

/// Population statistics across all analyzed systems
#[derive(Debug, Clone)]
pub struct PopulationStatistics<R>
where
    R: RealField,
{
    pub systems_analyzed: usize,
    pub systems_passing: usize,
    pub systems_excluded: usize,
    pub mean_ratio: R,
    pub std_ratio: R,
    pub median_ratio: R,
    pub field_strength_min: R,
    pub field_strength_max: R,
}

impl<R> PopulationStatistics<R>
where
    R: RealField + Clone + From<f64> + Into<f64>,
{
    pub fn new() -> Self {
        Self {
            systems_analyzed: 0,
            systems_passing: 0,
            systems_excluded: 0,
            mean_ratio: R::zero(),
            std_ratio: R::zero(),
            median_ratio: R::zero(),
            field_strength_min: R::from(1e308_f64),
            field_strength_max: R::from(-1e308_f64),
        }
    }

    pub fn compute_from_results(results: &[BinarySystemResult<R>]) -> Self {
        if results.is_empty() {
            return Self::new();
        }

        let ratios: Vec<R> = results.iter().map(|r| r.energy_balance_ratio).collect();
        let n = R::from(ratios.len() as f64);

        let sum: R = ratios.iter().cloned().fold(R::zero(), |acc, x| acc + x);
        let mean = sum / n;

        let variance_sum: R = ratios
            .iter()
            .cloned()
            .map(|r| {
                let diff = r - mean;
                diff * diff
            })
            .fold(R::zero(), |acc, x| acc + x);
        let variance = variance_sum / n;
        let std = variance.sqrt();

        let mut sorted_ratios = ratios.clone();
        sorted_ratios.sort_by(|a, b| {
            let a_f64: f64 = (*a).into();
            let b_f64: f64 = (*b).into();
            a_f64
                .partial_cmp(&b_f64)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let median = sorted_ratios[sorted_ratios.len() / 2];

        let field_min = results
            .iter()
            .map(|r| r.field_strength_mean)
            .fold(R::from(1e308_f64), |acc, x| if x < acc { x } else { acc });
        let field_max = results
            .iter()
            .map(|r| r.field_strength_mean)
            .fold(R::from(-1e308_f64), |acc, x| if x > acc { x } else { acc });

        let passing = results.iter().filter(|r| r.passes_criterion).count();

        Self {
            systems_analyzed: results.len(),
            systems_passing: passing,
            systems_excluded: 0, // Set externally
            mean_ratio: mean,
            std_ratio: std,
            median_ratio: median,
            field_strength_min: field_min,
            field_strength_max: field_max,
        }
    }
}

impl<R> Default for PopulationStatistics<R>
where
    R: RealField + Clone + From<f64> + Into<f64>,
{
    fn default() -> Self {
        Self::new()
    }
}

/// Timing residual from NANOGrav data
#[derive(Debug, Clone)]
pub struct TimingResidual<R>
where
    R: RealField,
{
    pub mjd: R,
    pub residual_ns: R,
    pub uncertainty_ns: R,
}
