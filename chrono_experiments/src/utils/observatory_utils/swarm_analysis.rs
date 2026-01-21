/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Gravitational Wave Analysis utilities.
//!
//! Provides advanced GW detection methods:
//! - Geodesic deviation-based strain calculation
//! - Hellings-Downs angular correlation for GW background detection
//!   Ported with R: RealField generics.

use crate::utils::observatory_utils::linear_algebra::{
    DynamicMatrix, DynamicVector, Matrix, Vector, covariance_dynamic, jacobi_eigenvalues_dynamic,
};
use crate::utils::observatory_utils::residual_filter::apply_single_sat_cleaning;
use chrono_data_manager::ObservatoryEpoch;
use deep_causality_num::RealField;
use deep_causality_physics::GALILEO_NOMINAL_RADIUS_M;
use deep_causality_physics::SPEED_OF_LIGHT;

/// GW Analysis Results
#[derive(Debug, Clone)]
pub struct GwAnalysisResults<R> {
    /// Strain from geodesic deviation (h = Δa / a)
    pub geodesic_strain: R,
    /// Strain from timing residuals
    pub timing_strain: R,
    /// Hellings-Downs correlation coefficient
    pub hellings_downs_correlation: R,
    /// Angular separation between satellite pairs (radians)
    pub angular_separation_rad: R,
    /// Expected HD correlation for isotropic GW background
    pub expected_hd_correlation: R,
    /// Significance: how well measured matches HD expectation
    pub hd_significance: R,
}

/// Calculate GW strain from geodesic deviation principle.
/// h ≈ c × Δτ_rate / L where Δτ_rate is in seconds/second
pub fn calculate_geodesic_strain<R>(clock_rate_diff_ns_s: R, baseline_m: R) -> R
where
    R: RealField + Copy + From<f64>,
{
    if baseline_m <= R::zero() {
        return R::zero();
    }
    // h ≈ c × Δτ_rate / L where Δτ_rate is in seconds/second
    let clock_rate_s_s = clock_rate_diff_ns_s * R::from(1e-9);
    (R::from(SPEED_OF_LIGHT) * clock_rate_s_s) / baseline_m
}

/// Calculate the Hellings-Downs correlation function.
pub fn hellings_downs_correlation<R>(theta_rad: R) -> R
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    // Handle the special case of zero separation (auto-correlation)
    if theta_rad.abs() < R::from(1e-10) {
        return R::one();
    }

    let cos_theta = theta_rad.cos();
    let x = (R::one() - cos_theta) / R::from(2.0);

    if x <= R::zero() {
        return R::from(0.5); // θ = 0 case
    }

    // HD correlation: Γ(θ) = (3/2)x·ln(x) - (1/4)x + (1/2)
    R::from(1.5) * x * x.ln() - R::from(0.25) * x + R::from(0.5)
}

/// Calculate angular separation between two position vectors.
pub fn calculate_angular_separation<R>(pos1: [R; 3], pos2: [R; 3]) -> R
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let r1_sq = pos1[0] * pos1[0] + pos1[1] * pos1[1] + pos1[2] * pos1[2];
    let r2_sq = pos2[0] * pos2[0] + pos2[1] * pos2[1] + pos2[2] * pos2[2];

    if r1_sq < R::from(1e-10) || r2_sq < R::from(1e-10) {
        return R::zero();
    }

    let r1 = r1_sq.sqrt();
    let r2 = r2_sq.sqrt();

    // Dot product
    let dot = pos1[0] * pos2[0] + pos1[1] * pos2[1] + pos1[2] * pos2[2];

    // cos(θ) = (r1 · r2) / (|r1| × |r2|)
    let cos_theta = (dot / (r1 * r2)).clamp(R::from(-1.0), R::one());

    cos_theta.acos()
}

/// Perform complete GW analysis for satellite pair.
pub fn analyze_gw_signal<R>(
    clock_rate_residual_ns_s: R,
    cross_correlation: R,
    pos1: [R; 3],
    pos2: [R; 3],
    baseline_m: R,
) -> GwAnalysisResults<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    // 1. Calculate angular separation
    let angular_separation = calculate_angular_separation(pos1, pos2);

    // 2. Calculate expected HD correlation for this angle
    let expected_hd = hellings_downs_correlation(angular_separation);

    // 3. Calculate geodesic strain
    let geodesic_strain = calculate_geodesic_strain(clock_rate_residual_ns_s, baseline_m);

    // 4. Calculate timing strain (alternative method)
    let timing_residual_s = clock_rate_residual_ns_s * R::from(1e-9);
    let timing_strain = (R::from(SPEED_OF_LIGHT) * timing_residual_s) / baseline_m;

    // 5. Calculate HD significance
    let hd_diff = (cross_correlation - expected_hd).abs();
    let two = R::from(2.0);
    let hd_significance = if hd_diff < two {
        (R::one() - hd_diff / two) * R::from(100.0)
    } else {
        R::zero()
    };

    GwAnalysisResults {
        geodesic_strain,
        timing_strain,
        hellings_downs_correlation: cross_correlation,
        angular_separation_rad: angular_separation,
        expected_hd_correlation: expected_hd,
        hd_significance,
    }
}

/// Convert radians to degrees
pub fn rad_to_deg<R>(rad: R) -> R
where
    R: RealField + Copy + From<f64>,
{
    rad * R::from(180.0) / R::from(std::f64::consts::PI)
}

// ============================================================================
// FREQUENCY DOMAIN ANALYSIS (PSD)
// ============================================================================

/// A spectral component from PSD analysis
#[derive(Debug, Clone)]
pub struct SpectralComponent<R> {
    pub frequency_hz: R,
    pub power: R,
    pub period_hr: R,
}

/// Calculate Power Spectral Density (PSD) using a simple DFT.
pub fn calculate_psd<R>(time_values: &[R], signal_values: &[R]) -> Vec<SpectralComponent<R>>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = time_values.len();
    if n < 2 {
        return vec![];
    }

    // Determine frequency range
    let duration = *time_values.last().unwrap() - *time_values.first().unwrap();
    let min_freq = R::one() / duration; // Fundamental frequency

    // Average sampling rate for Nyquist
    let avg_dt = duration / (R::from(n as f64) - R::one());
    let max_freq = R::from(0.5) / avg_dt;

    // Generate frequency bins (logarithmic spacing for wide dynamic range)
    let mut freqs = Vec::new();
    let mut f = min_freq;
    let max_phys_freq = R::from(0.033); // ~30s period

    // Check if max_freq < max_phys_freq
    let limit_freq = if max_freq < max_phys_freq {
        max_freq
    } else {
        max_phys_freq
    };

    while f < limit_freq {
        freqs.push(f);
        f *= R::from(1.05); // 5% logarithmic steps
    }

    // Apply Hanning Window to reduce spectral leakage
    let two_pi = R::from(2.0 * std::f64::consts::PI);
    let n_minus_1 = R::from((n - 1) as f64);

    let windowed_signal: Vec<R> = signal_values
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            let cos_term = (two_pi * R::from(i as f64) / n_minus_1).cos();
            let w = R::from(0.5) * (R::one() - cos_term);
            x * w
        })
        .collect();

    let mut spectrum = Vec::with_capacity(freqs.len());
    let window_sum_sq: R = windowed_signal
        .iter()
        .fold(R::zero(), |acc, &w| acc + w * w);

    let norm_factor = if window_sum_sq > R::zero() {
        R::from(2.0) / window_sum_sq
    } else {
        R::one()
    };

    let n_real = R::from(n as f64);

    for &f in &freqs {
        let mut re = R::zero();
        let mut im = R::zero();
        for (i, &val) in windowed_signal.iter().enumerate() {
            let t = time_values[i];
            let angle = -two_pi * f * t;
            re += val * angle.cos();
            im += val * angle.sin();
        }

        let power = (re * re + im * im) * norm_factor / n_real;

        spectrum.push(SpectralComponent {
            frequency_hz: f,
            power,
            period_hr: (R::one() / f) / R::from(3600.0),
        });
    }

    spectrum
}

/// Detect significant spectral peaks in specific bands.
#[derive(Debug, Clone, Default)]
pub struct SpectralPeaks<R> {
    pub peak_12h_power: R,       // Earth orbit (Kepler)
    pub peak_14h_power: R,       // Galileo orbit
    pub peak_27d_power: R,       // Moon
    pub nhz_background_power: R, // Nanohertz hum (GW background)
}

pub fn analyze_spectral_peaks<R>(spectrum: &[SpectralComponent<R>]) -> SpectralPeaks<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let mut p12h = R::zero();
    let mut p14h = R::zero();
    let mut p27d = R::zero();
    let mut p_nhz = R::zero();

    let hr_11_5 = R::from(11.5);
    let hr_12_5 = R::from(12.5);
    let hr_13_5 = R::from(13.5);
    let hr_14_5 = R::from(14.5);
    let hr_600 = R::from(600.0);
    let hr_700 = R::from(700.0);
    let hz_1e9 = R::from(1e-9);
    let hz_1e8 = R::from(1e-8);

    for comp in spectrum {
        // 12-hour band (11.5 - 12.5h)
        if comp.period_hr >= hr_11_5 && comp.period_hr <= hr_12_5 && comp.power > p12h {
            p12h = comp.power;
        }
        // 14-hour band (13.5 - 14.5h)
        if comp.period_hr >= hr_13_5 && comp.period_hr <= hr_14_5 && comp.power > p14h {
            p14h = comp.power;
        }
        // 27-day band (25 - 29 days) -> 600 - 700 hours
        if comp.period_hr >= hr_600 && comp.period_hr <= hr_700 && comp.power > p27d {
            p27d = comp.power;
        }
        // Nanohertz band (1e-9 to 1e-8 Hz)
        if comp.frequency_hz >= hz_1e9 && comp.frequency_hz <= hz_1e8 {
            p_nhz += comp.power;
        }
    }

    SpectralPeaks {
        peak_12h_power: p12h,
        peak_14h_power: p14h,
        peak_27d_power: p27d,
        nhz_background_power: p_nhz,
    }
}

// ============================================================================
// SWARM COVARIANCE ANALYSIS (5-Satellite)
// ============================================================================

/// Results of the 5-satellite covariance analysis
#[derive(Debug, Clone)]
pub struct SwarmCovarianceResult<R> {
    pub eigenvalues: Vector<R>,
    pub eigenvectors: Matrix<R>,
    pub quadrupole_power: R, // Sum of specific quadrupole modes
    pub monopole_power: R,   // Common mode noise
}

/// Analyze the covariance of the 5-satellite swarm residuals.
pub fn analyze_swarm_vectors<R>(
    epochs: &[ObservatoryEpoch],
    target_sats: &[&str],
) -> DynamicSwarmResult<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let mut aligned_timestamps: Vec<i64> = Vec::new();
    let n_sats = target_sats.len();
    let mut sat_data_biases: Vec<Vec<R>> = vec![Vec::new(); n_sats];
    let mut sat_data_radii: Vec<Vec<R>> = vec![Vec::new(); n_sats];

    let mut sorted: Vec<&ObservatoryEpoch> = epochs.iter().collect();
    sorted.sort_by_key(|e| e.timestamp);

    for epoch in sorted {
        let mut present = true;
        for id in target_sats {
            if !epoch.satellites.contains_key(*id) {
                present = false;
                break;
            }
        }

        if present {
            aligned_timestamps.push(epoch.timestamp);
            for (i, id) in target_sats.iter().enumerate() {
                let s = &epoch.satellites[*id];
                sat_data_biases[i].push(R::from(s.clock_bias_ns));
                sat_data_radii[i].push(R::from(s.radius_m));
            }
        }
    }

    if aligned_timestamps.len() < 10 {
        return DynamicSwarmResult {
            eigenvalues: vec![],
            eigenvectors: vec![],
            quadrupole_power: R::zero(),
            monopole_power: R::zero(),
            num_satellites: 0,
        };
    }

    // Clean residuals using the residual filter
    let r_nominal = R::from(GALILEO_NOMINAL_RADIUS_M);
    let mut cleaned_vectors: Vec<Vec<R>> = Vec::new();

    for i in 0..n_sats {
        let cleaned = apply_single_sat_cleaning(
            &aligned_timestamps,
            &sat_data_biases[i],
            &sat_data_radii[i],
            r_nominal,
        );
        cleaned_vectors.push(cleaned.gw_candidates);
    }

    // Build residual rows for covariance analysis
    let n_epochs = cleaned_vectors[0].len();
    let mut residual_rows: Vec<DynamicVector<R>> = Vec::with_capacity(n_epochs);

    for t in 0..n_epochs {
        let mut row_vec: Vec<R> = Vec::with_capacity(n_sats);
        for sat_vec in &cleaned_vectors {
            row_vec.push(sat_vec[t]);
        }
        residual_rows.push(row_vec);
    }

    // Perform covariance analysis
    analyze_swarm_covariance_dynamic(&residual_rows)
}

// ============================================================================
// DYNAMIC SWARM COVARIANCE ANALYSIS (N-Satellite)
// ============================================================================

/// Results of the N-satellite covariance analysis
#[derive(Debug, Clone)]
pub struct DynamicSwarmResult<R> {
    pub eigenvalues: DynamicVector<R>,  // Length N
    pub eigenvectors: DynamicMatrix<R>, // N x N
    pub quadrupole_power: R,
    pub monopole_power: R,
    pub num_satellites: usize,
}

pub fn analyze_swarm_covariance_dynamic<R>(residuals: &[DynamicVector<R>]) -> DynamicSwarmResult<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    if residuals.is_empty() {
        return DynamicSwarmResult {
            eigenvalues: vec![],
            eigenvectors: vec![],
            quadrupole_power: R::zero(),
            monopole_power: R::zero(),
            num_satellites: 0,
        };
    }

    let cov = covariance_dynamic(residuals);
    let (eig_vals, eig_vecs) = jacobi_eigenvalues_dynamic(&cov);
    let dim = eig_vals.len();

    let mut monopole_pwr = R::zero();
    let mut quadrupole_pwr = R::zero();

    for i in 0..dim {
        let val = eig_vals[i];

        // Vector Mode Analysis
        let mut sum_components = R::zero();
        let mut abs_sum = R::zero();

        for item in eig_vecs.iter().take(dim) {
            let comp = item[i];
            sum_components += comp;
            abs_sum += comp.abs();
        }

        // Monopole Check (Coherent Sign)
        let eps = R::from(1e-9);
        let threshold = R::from(0.8);
        if abs_sum > eps && (sum_components.abs() / abs_sum) > threshold {
            monopole_pwr += val;
        } else {
            quadrupole_pwr += val;
        }
    }

    DynamicSwarmResult {
        eigenvalues: eig_vals,
        eigenvectors: eig_vecs,
        quadrupole_power: quadrupole_pwr,
        monopole_power: monopole_pwr,
        num_satellites: dim,
    }
}
