use crate::observatory_utils::AnomalyTier;
use chrono_data_manager::SatelliteState;

#[derive(Debug, Clone)]
pub struct IncidentRecord {
    pub timestamp: i64,
    pub date_str: String,
    pub satellite: String,
    pub pc2_variance: f64,
    pub tier: AnomalyTier,
}

pub struct ObservatoryConfigParams<R> {
    pub g_ref: R,
    pub c: R,
    pub gm_earth: R,
    pub time_delay_s: u64,
}

/// Analysis results for 5T observatory measurements with GW detection metrics
#[derive(Debug, Clone)]
pub struct ObservatoryResults {
    pub total_epochs: usize,
    pub valid_5t_epochs: usize,
    // Differential observables
    pub mean_clock_diff_e14_e12: f64,
    pub std_clock_diff_e14_e12: f64,
    pub mean_altitude_diff_m: f64,
    // GR verification (raw - includes systematics)
    pub predicted_gr_shift_ns: f64,
    pub observed_shift_ns: f64,
    pub gr_agreement_pct: f64,
    // ISOLATED GR Analysis (systematics removed)
    pub clock_rate_altitude_correlation: f64,
    pub measured_gr_coefficient: f64,
    pub predicted_gr_coefficient: f64,
    pub gr_coefficient_agreement_pct: f64,
    pub detrended_clock_std_ns: f64,
    // Correlation metrics
    pub e12_e18_correlation: f64,
    pub temporal_autocorr_e14: f64,
    // GW Detection metrics (raw)
    pub gw_strain_upper_limit: f64,
    pub quadrupole_residual: f64,
    pub cross_correlation_coeff: f64,
    pub spectral_power_n_hz: f64,
    // Hellings-Downs analysis
    pub geodesic_strain: f64,
    pub angular_separation_deg: f64,
    pub expected_hd_correlation: f64,
    pub hd_significance_pct: f64,
    // GQCD Cleaning Pipeline (NEW)
    pub raw_rms_ns: f64,            // RMS of raw clock differences
    pub physics_rms_ns: f64,        // RMS after GR+SR subtraction
    pub cleaned_rms_ns: f64,        // RMS after polynomial detrend (GW candidate)
    pub cleaning_ratio: f64,        // raw_rms / cleaned_rms (improvement factor)
    pub red_noise_coefficient: f64, // Lag-1 autocorr (>0.5 = red noise = potential GW)
    pub drift_rate_ns_hr: f64,      // Hardware drift rate (from polynomial fit)
    // Frequency Domain Analysis (PSD)
    pub peak_12h_power: f64,       // Earth orbit
    pub peak_14h_power: f64,       // Galileo orbit
    pub peak_27d_power: f64,       // Moon
    pub nhz_background_power: f64, // Nanohertz hum
    // Swarm Analysis (7-Sat Covariance)
    pub swarm_quadrupole_power: f64, // Mode 3+ signal
    pub swarm_monopole_power: f64,   // Mode 1 noise
    // Anomaly Timing
    pub anomaly_start: String,
    pub anomaly_end: String,
}

/// 5T Field measurement for gravitational wave detection.
///
/// Represents a complete 5-temporal measurement configuration:
/// - T1: E12 (circular orbit) - Reference baseline
/// - T2: E18 (eccentric orbit) - Cross-correlation baseline  
/// - T3: E14 at t₀ (eccentric) - Dynamic clock
/// - T4: E14 at t₀+1h - Temporal offset 1
/// - T5: E14 at t₀+2h - Temporal offset 2
#[derive(Debug, Clone)]
pub struct FiveT {
    pub timestamp: i64,
    // Reference clocks (stable baselines)
    pub t1_e12: Option<SatelliteState>,
    pub t2_e18: Option<SatelliteState>,
    // Dynamic clock time series
    pub t3_e14_t0: Option<SatelliteState>,
    pub t4_e14_t1h: Option<SatelliteState>,
    pub t5_e14_t2h: Option<SatelliteState>,
}

impl FiveT {
    /// Create a new empty 5T field at the given timestamp.
    pub fn new(timestamp: i64) -> Self {
        Self {
            timestamp,
            t1_e12: None,
            t2_e18: None,
            t3_e14_t0: None,
            t4_e14_t1h: None,
            t5_e14_t2h: None,
        }
    }

    /// Check if the field has valid differential measurements (E14 and at least one reference).
    pub fn has_valid_differential(&self) -> bool {
        self.t3_e14_t0.is_some() && (self.t1_e12.is_some() || self.t2_e18.is_some())
    }

    /// Check if all 5T components are present.
    pub fn is_complete(&self) -> bool {
        self.t1_e12.is_some()
            && self.t2_e18.is_some()
            && self.t3_e14_t0.is_some()
            && self.t4_e14_t1h.is_some()
            && self.t5_e14_t2h.is_some()
    }
}

/// Observatory satellite configuration
/// (satellite_id, orbit_type, role)
#[derive(Debug, Clone, Copy)]
pub struct SatelliteConfig {
    pub id: &'static str,
    pub orbit_type: OrbitType,
    pub role: SatelliteRole,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitType {
    Circular,  // Stable orbit, constant altitude
    Eccentric, // Varying altitude (E14, E18 launched into wrong orbit)
}

#[derive(Debug, Clone, Copy)]
pub enum SatelliteRole {
    Reference, // Stable clock for baseline
    Dynamic,   // Clock experiencing varying gravitational potential
}

/// 7-Satellite Swarm Configuration
/// Order matches legacy for PCA loading consistency.
pub const OBSERVATORY_CONFIG: [SatelliteConfig; 7] = [
    // E14 (Eccentric - Dynamic Probe)
    SatelliteConfig {
        id: "E14",
        orbit_type: OrbitType::Eccentric,
        role: SatelliteRole::Dynamic,
    },
    // E18 (Eccentric - Dynamic Probe)
    SatelliteConfig {
        id: "E18",
        orbit_type: OrbitType::Eccentric,
        role: SatelliteRole::Dynamic,
    },
    // E12 (Circular - Reference Plane B)
    SatelliteConfig {
        id: "E12",
        orbit_type: OrbitType::Circular,
        role: SatelliteRole::Reference,
    },
    // E11 (Circular - Reference Plane A)
    SatelliteConfig {
        id: "E11",
        orbit_type: OrbitType::Circular,
        role: SatelliteRole::Reference,
    },
    // E19 (Circular - Reference Plane A)
    SatelliteConfig {
        id: "E19",
        orbit_type: OrbitType::Circular,
        role: SatelliteRole::Reference,
    },
    // E08 (Circular - Reference Plane C)
    SatelliteConfig {
        id: "E08",
        orbit_type: OrbitType::Circular,
        role: SatelliteRole::Reference,
    },
    // E09 (Circular - Reference Plane C)
    SatelliteConfig {
        id: "E09",
        orbit_type: OrbitType::Circular,
        role: SatelliteRole::Reference,
    },
];
