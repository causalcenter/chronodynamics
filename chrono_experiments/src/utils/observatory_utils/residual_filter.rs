/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! GQCD Residual Filter for GW Detection.
//!
//! Implements the 3-stage cleaning pipeline:
//! 1. Subtract known physics (GR + SR from theoretical model)
//! 2. Subtract hardware drift (polynomial detrending)
//! 3. Analyze cleaned residuals for GW signatures
use deep_causality_num::RealField;
use deep_causality_physics::{EARTH_GM, SPEED_OF_LIGHT};

/// Predicted relativistic clock difference based on known physics.
#[derive(Debug, Clone, Copy)]
pub struct TheoreticalPrediction<R> {
    /// GR gravitational time dilation (ns)
    pub gr_shift_ns: R,
    /// SR velocity time dilation (ns)
    pub sr_shift_ns: R,
    /// Total theoretical shift (ns)
    pub total_ns: R,
}

/// Calculate theoretical clock difference due to GR + SR.
pub fn calculate_theoretical_diff<R>(
    r1_m: R,
    r2_m: R,
    v1_ms: Option<R>,
    v2_ms: Option<R>,
    dt_s: R,
) -> TheoreticalPrediction<R>
where
    R: RealField + Copy + From<f64>,
{
    let c = R::from(SPEED_OF_LIGHT);
    let c_sq = c * c;
    let gm_earth = R::from(EARTH_GM);

    // GR: Δτ/τ = GM/c² × (1/r₂ - 1/r₁)
    // Clock at lower altitude (r1 < r2) runs slower, so bias (t1 - t2) decreases.
    // (1/r2 - 1/r1) is negative when r1 < r2.
    let r1_inv = R::one() / r1_m;
    let r2_inv = R::one() / r2_m;
    let gr_fractional = (gm_earth / c_sq) * (r2_inv - r1_inv);
    let gr_shift_ns = gr_fractional * dt_s * R::from(1e9); // Convert to ns

    // SR: Δτ/τ = -(v₁² - v₂²) / (2c²)
    let sr_fractional = match (v1_ms, v2_ms) {
        (Some(v1), Some(v2)) => -(v1 * v1 - v2 * v2) / (R::from(2.0) * c_sq),
        _ => R::zero(),
    };
    let sr_shift_ns = sr_fractional * dt_s * R::from(1e9);

    TheoreticalPrediction {
        gr_shift_ns,
        sr_shift_ns,
        total_ns: gr_shift_ns + sr_shift_ns,
    }
}

/// Fit a quadratic polynomial to data and return coefficients.
/// y = a*x² + b*x + c
pub fn fit_quadratic<R>(x: &[R], y: &[R]) -> (R, R, R)
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = x.len().min(y.len());
    if n < 3 {
        let sum_y: R = y.iter().fold(R::zero(), |acc, &val| acc + val);
        return (R::zero(), R::zero(), sum_y / R::from(n.max(1) as f64));
    }

    // Normal equations for least squares
    let mut sum_x = R::zero();
    let mut sum_x2 = R::zero();
    let mut sum_x3 = R::zero();
    let mut sum_x4 = R::zero();
    let mut sum_y = R::zero();
    let mut sum_xy = R::zero();
    let mut sum_x2y = R::zero();

    for i in 0..n {
        let xi = x[i];
        let yi = y[i];
        let x2 = xi * xi;
        let x3 = x2 * xi;
        let x4 = x3 * xi;

        sum_x += xi;
        sum_x2 += x2;
        sum_x3 += x3;
        sum_x4 += x4;
        sum_y += yi;
        sum_xy += xi * yi;
        sum_x2y += x2 * yi;
    }

    let nf = R::from(n as f64);

    // Solve 3x3 system using Cramer's rule
    let det = sum_x4 * (sum_x2 * nf - sum_x * sum_x) - sum_x3 * (sum_x3 * nf - sum_x * sum_x2)
        + sum_x2 * (sum_x3 * sum_x - sum_x2 * sum_x2);

    if det.abs() < R::from(1e-20) {
        // Nearly singular - fall back to linear or constant
        return (R::zero(), R::zero(), sum_y / nf);
    }

    let a = (sum_x2y * (sum_x2 * nf - sum_x * sum_x) - sum_xy * (sum_x3 * nf - sum_x * sum_x2)
        + sum_y * (sum_x3 * sum_x - sum_x2 * sum_x2))
        / det;

    let b = (sum_x4 * (sum_xy * nf - sum_y * sum_x) - sum_x3 * (sum_x2y * nf - sum_y * sum_x2)
        + sum_x2 * (sum_x2y * sum_x - sum_xy * sum_x2))
        / det;

    let c = (sum_x4 * (sum_x2 * sum_y - sum_x * sum_xy)
        - sum_x3 * (sum_x3 * sum_y - sum_x * sum_x2y)
        + sum_x2 * (sum_x3 * sum_xy - sum_x2 * sum_x2y))
        / det;

    (a, b, c)
}

/// Remove polynomial trend from data.
/// Fits a quadratic and subtracts it, returning the residuals.
pub fn detrend_polynomial<R>(timestamps: &[i64], values: &[R]) -> Vec<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = timestamps.len().min(values.len());
    if n == 0 {
        return vec![];
    }

    // Normalize timestamps to avoid numerical issues
    let t0 = timestamps[0] as f64;
    let x: Vec<R> = timestamps
        .iter()
        .take(n)
        .map(|&t| R::from((t as f64 - t0) / 3600.0))
        .collect(); // Hours

    let values_slice = &values[..n];
    let (a, b, c) = fit_quadratic(&x, values_slice);

    // Subtract trend
    x.iter()
        .zip(values_slice.iter())
        .map(|(&xi, &yi)| yi - (a * xi * xi + b * xi + c))
        .collect()
}

/// Results from the data cleaning pipeline
#[derive(Debug, Clone)]
pub struct CleanedResiduals<R> {
    /// Timestamps
    pub timestamps: Vec<i64>,
    /// Raw observed clock differences (ns)
    pub raw_diffs: Vec<R>,
    /// Physics residuals (raw - GR - SR model)
    pub physics_residuals: Vec<R>,
    /// Cleaned residuals (physics - polynomial drift)
    pub gw_candidates: Vec<R>,
    /// Polynomial coefficients (a, b, c)
    pub drift_coefficients: (R, R, R),
    /// RMS of raw diffs
    pub raw_rms_ns: R,
    /// RMS of physics residuals
    pub physics_rms_ns: R,
    /// RMS of GW candidates (should be much smaller)
    pub gw_rms_ns: R,
}

/// Apply the full data cleaning pipeline.
pub fn apply_gqcd_cleaning<R>(
    timestamps: &[i64],
    raw_diffs: &[R],
    r1_m: R,
    r2_m: R,
    v1_ms: Option<R>,
    v2_ms: Option<R>,
) -> CleanedResiduals<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = timestamps.len().min(raw_diffs.len());

    // Step 1: Subtract known physics
    let mut physics_residuals = Vec::with_capacity(n);

    for i in 0..n {
        // Calculate time interval from previous measurement
        let dt_s = if i > 0 {
            R::from((timestamps[i] - timestamps[i - 1]) as f64)
        } else {
            R::from(30.0) // Default 30 second epochs
        };

        let prediction = calculate_theoretical_diff(r1_m, r2_m, v1_ms, v2_ms, dt_s);

        // For accumulated bias, we need the cumulative effect
        let cumulative_theoretical = prediction.total_ns * R::from(i as f64);

        physics_residuals.push(raw_diffs[i] - cumulative_theoretical);
    }

    // Step 2: Subtract hardware drift (polynomial)
    let gw_candidates = detrend_polynomial(timestamps, &physics_residuals);

    let time_vec: Vec<R> = timestamps
        .iter()
        .map(|&t| R::from((t - timestamps[0]) as f64 / 3600.0))
        .collect();

    let drift_coefficients = fit_quadratic(&time_vec, &physics_residuals);

    // RMS Calculation
    let zero = R::zero();
    let n_f64 = R::from(n.max(1) as f64);

    let raw_sum_sq = raw_diffs.iter().fold(zero, |acc, &x| acc + x * x);
    let raw_rms_ns = (raw_sum_sq / n_f64).sqrt();

    let phys_sum_sq = physics_residuals.iter().fold(zero, |acc, &x| acc + x * x);
    let physics_rms_ns = (phys_sum_sq / n_f64).sqrt();

    let gw_len_f64 = R::from(gw_candidates.len().max(1) as f64);
    let gw_sum_sq = gw_candidates.iter().fold(zero, |acc, &x| acc + x * x);
    let gw_rms_ns = (gw_sum_sq / gw_len_f64).sqrt();

    CleanedResiduals {
        timestamps: timestamps.to_vec(),
        raw_diffs: raw_diffs.to_vec(),
        physics_residuals,
        gw_candidates,
        drift_coefficients,
        raw_rms_ns,
        physics_rms_ns,
        gw_rms_ns,
    }
}

/// Apply GQCD cleaning to a single satellite's clock data.
pub fn apply_single_sat_cleaning<R>(
    timestamps: &[i64],
    clock_values: &[R],
    radii_m: &[R],
    r_nominal_m: R,
) -> CleanedResiduals<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = timestamps.len().min(clock_values.len());
    let c = R::from(SPEED_OF_LIGHT);
    let c_sq = c * c;
    let gm_earth = R::from(EARTH_GM);

    // Step 1: Subtract Physics (GR variation)
    let mut physics_residuals = Vec::with_capacity(n);
    let mut cumulative_physics = R::zero();

    for i in 0..n {
        let dt_s = if i > 0 {
            R::from((timestamps[i] - timestamps[i - 1]) as f64)
        } else {
            R::from(30.0)
        };

        let r_actual = radii_m[i];

        // GR: Δτ/τ = GM/c² × (1/r_nominal - 1/r_actual)
        let r_nominal_inv = R::one() / r_nominal_m;
        let r_actual_inv = R::one() / r_actual;

        let gr_fractional = (gm_earth / c_sq) * (r_nominal_inv - r_actual_inv);
        let gr_shift_ns = gr_fractional * dt_s * R::from(1e9);

        cumulative_physics += gr_shift_ns;

        physics_residuals.push(clock_values[i] - cumulative_physics);
    }

    // Step 2: Polynomial Detrend (Drift + IGS Reference Trend)
    let gw_candidates = detrend_polynomial(timestamps, &physics_residuals);

    let time_vec: Vec<R> = timestamps
        .iter()
        .map(|&t| R::from((t - timestamps[0]) as f64 / 3600.0))
        .collect();

    let drift_coefficients = fit_quadratic(&time_vec, &physics_residuals);

    // RMS Calculation
    let zero = R::zero();
    let n_f64 = R::from(n.max(1) as f64);

    let raw_sum_sq = clock_values.iter().fold(zero, |acc, &x| acc + x * x);
    let raw_rms_ns = (raw_sum_sq / n_f64).sqrt();

    let phys_sum_sq = physics_residuals.iter().fold(zero, |acc, &x| acc + x * x);
    let physics_rms_ns = (phys_sum_sq / n_f64).sqrt();

    let gw_sum_sq = gw_candidates.iter().fold(zero, |acc, &x| acc + x * x);
    let gw_rms_ns = (gw_sum_sq / n_f64).sqrt();

    CleanedResiduals {
        timestamps: timestamps.to_vec(),
        raw_diffs: clock_values.to_vec(),
        physics_residuals,
        gw_candidates,
        drift_coefficients,
        raw_rms_ns,
        physics_rms_ns,
        gw_rms_ns,
    }
}

/// Check if residuals show "red noise" (correlated low-frequency wobble).
pub fn check_for_red_noise<R>(residuals: &[R]) -> R
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    if residuals.len() < 3 {
        return R::zero();
    }

    let n = residuals.len();
    let n_f64 = R::from(n as f64);
    let sum = residuals.iter().fold(R::zero(), |acc, &x| acc + x);
    let mean = sum / n_f64;

    let mut cov_0 = R::zero(); // Variance
    let mut cov_1 = R::zero(); // Lag-1 covariance

    for i in 0..n {
        let xi = residuals[i] - mean;
        cov_0 += xi * xi;
        if i < n - 1 {
            let xi1 = residuals[i + 1] - mean;
            cov_1 += xi * xi1;
        }
    }

    if cov_0 < R::from(1e-20) {
        return R::zero();
    }

    cov_1 / cov_0
}
