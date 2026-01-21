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

use deep_causality_physics::{EARTH_GM, SPEED_OF_LIGHT};

/// Predicted relativistic clock difference based on known physics.
#[derive(Debug, Clone, Copy)]
pub struct TheoreticalPrediction {
    /// GR gravitational time dilation (ns)
    pub gr_shift_ns: f64,
    /// SR velocity time dilation (ns)
    pub sr_shift_ns: f64,
    /// Total theoretical shift (ns)
    pub total_ns: f64,
}

/// Calculate theoretical clock difference due to GR + SR.
///
/// Uses the full relativistic formula for gravitational and velocity time dilation.
///
/// # Arguments
/// * `r1_m` - Radius of clock 1 from Earth center (meters)
/// * `r2_m` - Radius of clock 2 from Earth center (meters)  
/// * `v1_ms` - Velocity of clock 1 (m/s), if available
/// * `v2_ms` - Velocity of clock 2 (m/s), if available
/// * `dt_s` - Time interval for integration (seconds)
///
/// # Returns
/// Theoretical clock difference prediction
pub fn calculate_theoretical_diff(
    r1_m: f64,
    r2_m: f64,
    v1_ms: Option<f64>,
    v2_ms: Option<f64>,
    dt_s: f64,
) -> TheoreticalPrediction {
    let c = SPEED_OF_LIGHT;
    let c_sq = c * c;

    // GR: Δτ/τ = GM/c² × (1/r₂ - 1/r₁)
    // Clock at lower altitude (r1 < r2) runs slower, so bias (t1 - t2) decreases.
    // (1/r2 - 1/r1) is negative when r1 < r2.
    let gr_fractional = (EARTH_GM / c_sq) * (1.0 / r2_m - 1.0 / r1_m);
    let gr_shift_ns = gr_fractional * dt_s * 1e9; // Convert to ns

    // SR: Δτ/τ = -(v₁² - v₂²) / (2c²)
    // Faster clock runs slower (velocity time dilation)
    let sr_fractional = match (v1_ms, v2_ms) {
        (Some(v1), Some(v2)) => -(v1 * v1 - v2 * v2) / (2.0 * c_sq),
        _ => 0.0,
    };
    let sr_shift_ns = sr_fractional * dt_s * 1e9;

    TheoreticalPrediction {
        gr_shift_ns,
        sr_shift_ns,
        total_ns: gr_shift_ns + sr_shift_ns,
    }
}

/// Fit a quadratic polynomial to data and return coefficients.
/// y = a*x² + b*x + c
///
/// Uses simple least squares.
///
/// # Returns
/// (a, b, c) coefficients
pub fn fit_quadratic(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let n = x.len().min(y.len());
    if n < 3 {
        return (0.0, 0.0, y.iter().sum::<f64>() / n.max(1) as f64);
    }

    // Normal equations for least squares
    // [Σx⁴  Σx³  Σx²] [a]   [Σx²y]
    // [Σx³  Σx²  Σx ] [b] = [Σxy ]
    // [Σx²  Σx   n  ] [c]   [Σy  ]

    let mut sum_x = 0.0f64;
    let mut sum_x2 = 0.0f64;
    let mut sum_x3 = 0.0f64;
    let mut sum_x4 = 0.0f64;
    let mut sum_y = 0.0f64;
    let mut sum_xy = 0.0f64;
    let mut sum_x2y = 0.0f64;

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

    let nf = n as f64;

    // Solve 3x3 system using Cramer's rule
    let det = sum_x4 * (sum_x2 * nf - sum_x * sum_x) - sum_x3 * (sum_x3 * nf - sum_x * sum_x2)
        + sum_x2 * (sum_x3 * sum_x - sum_x2 * sum_x2);

    if det.abs() < 1e-20 {
        // Nearly singular - fall back to linear or constant
        return (0.0, 0.0, sum_y / nf);
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
///
/// Fits a quadratic and subtracts it, returning the residuals.
pub fn detrend_polynomial(timestamps: &[i64], values: &[f64]) -> Vec<f64> {
    let n = timestamps.len().min(values.len());
    if n == 0 {
        return vec![];
    }

    // Normalize timestamps to avoid numerical issues
    let t0 = timestamps[0] as f64;
    let x: Vec<f64> = timestamps
        .iter()
        .take(n)
        .map(|&t| (t as f64 - t0) / 3600.0)
        .collect(); // Hours

    let (a, b, c) = fit_quadratic(&x, &values[..n]);

    // Subtract trend
    x.iter()
        .zip(values.iter().take(n))
        .map(|(&xi, &yi)| yi - (a * xi * xi + b * xi + c))
        .collect()
}

/// Results from the GQCD cleaning pipeline
#[derive(Debug, Clone)]
pub struct CleanedResiduals {
    /// Timestamps
    pub timestamps: Vec<i64>,
    /// Raw observed clock differences (ns)
    pub raw_diffs: Vec<f64>,
    /// Physics residuals (raw - GR - SR model)
    pub physics_residuals: Vec<f64>,
    /// Cleaned residuals (physics - polynomial drift)
    pub gw_candidates: Vec<f64>,
    /// Polynomial coefficients (a, b, c)
    pub drift_coefficients: (f64, f64, f64),
    /// RMS of raw diffs
    pub raw_rms_ns: f64,
    /// RMS of physics residuals
    pub physics_rms_ns: f64,
    /// RMS of GW candidates (should be much smaller)
    pub gw_rms_ns: f64,
}

/// Apply the full GQCD cleaning pipeline.
///
/// 1. Subtract known physics (GR + SR)
/// 2. Subtract hardware drift (quadratic polynomial)
/// 3. Return cleaned residuals
///
/// # Arguments
/// * `timestamps` - Unix timestamps
/// * `raw_diffs` - Raw clock differences (ns)
/// * `r1_m` - Average radius of clock 1 (meters)
/// * `r2_m` - Average radius of clock 2 (meters)
/// * `v1_ms` - Average velocity of clock 1 (m/s), optional
/// * `v2_ms` - Average velocity of clock 2 (m/s), optional
pub fn apply_gqcd_cleaning(
    timestamps: &[i64],
    raw_diffs: &[f64],
    r1_m: f64,
    r2_m: f64,
    v1_ms: Option<f64>,
    v2_ms: Option<f64>,
) -> CleanedResiduals {
    let n = timestamps.len().min(raw_diffs.len());

    // Step 1: Subtract known physics
    let mut physics_residuals = Vec::with_capacity(n);

    for i in 0..n {
        // Calculate time interval from previous measurement
        let dt_s = if i > 0 {
            (timestamps[i] - timestamps[i - 1]) as f64
        } else {
            30.0 // Default 30 second epochs
        };

        let prediction = calculate_theoretical_diff(r1_m, r2_m, v1_ms, v2_ms, dt_s);

        // For accumulated bias, we need the cumulative effect
        // This is a simplification - real analysis would integrate the rate
        let cumulative_theoretical = prediction.total_ns * (i as f64);

        // Actually, for clock BIAS (not rate), we compare positions, not rates
        // The theoretical bias grows linearly with the rate difference
        physics_residuals.push(raw_diffs[i] - cumulative_theoretical);
    }

    // Step 2: Subtract hardware drift (polynomial)
    let gw_candidates = detrend_polynomial(timestamps, &physics_residuals);
    let drift_coefficients = fit_quadratic(
        &timestamps
            .iter()
            .map(|&t| (t - timestamps[0]) as f64 / 3600.0)
            .collect::<Vec<_>>(),
        &physics_residuals,
    );

    // Calculate RMS values
    let raw_rms_ns = (raw_diffs.iter().map(|x| x * x).sum::<f64>() / n.max(1) as f64).sqrt();
    let physics_rms_ns =
        (physics_residuals.iter().map(|x| x * x).sum::<f64>() / n.max(1) as f64).sqrt();
    let gw_rms_ns = (gw_candidates.iter().map(|x| x * x).sum::<f64>()
        / gw_candidates.len().max(1) as f64)
        .sqrt();

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
///
/// Removes the relativistic effect of orbital eccentricity (variation from nominal radius).
///
/// # Arguments
/// * `timestamps` - Unix timestamps
/// * `clock_values` - Raw clock biases (ns)
/// * `radii_m` - Radius at each timestamp (meters)
/// * `r_nominal_m` - Nominal semi-major axis (meters)
pub fn apply_single_sat_cleaning(
    timestamps: &[i64],
    clock_values: &[f64],
    radii_m: &[f64],
    r_nominal_m: f64,
) -> CleanedResiduals {
    let n = timestamps.len().min(clock_values.len());
    let c = SPEED_OF_LIGHT;
    let c_sq = c * c;

    // Step 1: Subtract Physics (GR variation)
    let mut physics_residuals = Vec::with_capacity(n);
    let mut cumulative_physics = 0.0;

    for i in 0..n {
        let dt_s = if i > 0 {
            (timestamps[i] - timestamps[i - 1]) as f64
        } else {
            30.0
        };

        let r_actual = radii_m[i];

        // GR: Δτ/τ = GM/c² × (1/r_nominal - 1/r_actual)
        // If r_actual < r_nominal, term is negative (clock slow)
        let gr_fractional = (EARTH_GM / c_sq) * (1.0 / r_nominal_m - 1.0 / r_actual);
        let gr_shift_ns = gr_fractional * dt_s * 1e9;

        cumulative_physics += gr_shift_ns;

        // Residual = Measured - Physics
        // Note: measured bias accumulates the physics effect naturally.
        physics_residuals.push(clock_values[i] - cumulative_physics);
    }

    // Step 2: Polynomial Detrend (Drift + IGS Reference Trend)
    // We remove the quadratic drift of the clock itself.
    // WARNING: This also removes the Monopole (IGS noise) if it looks quadratic!
    // However, GW signals in the ~hour band are higher frequency than the daily drift.
    // So extracting residuals after quadratic detrend preserves the relevant signal band.
    let gw_candidates = detrend_polynomial(timestamps, &physics_residuals);

    let drift_coefficients = fit_quadratic(
        &timestamps
            .iter()
            .map(|&t| (t - timestamps[0]) as f64 / 3600.0)
            .collect::<Vec<_>>(),
        &physics_residuals,
    );

    // RMS Calculation
    let raw_rms_ns = (clock_values.iter().map(|x| x * x).sum::<f64>() / n.max(1) as f64).sqrt();
    let physics_rms_ns =
        (physics_residuals.iter().map(|x| x * x).sum::<f64>() / n.max(1) as f64).sqrt();
    let gw_rms_ns = (gw_candidates.iter().map(|x| x * x).sum::<f64>() / n.max(1) as f64).sqrt();

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
///
/// Returns the lag-1 autocorrelation. Values close to 1 indicate red noise.
/// Values close to 0 indicate white noise.
pub fn check_for_red_noise(residuals: &[f64]) -> f64 {
    if residuals.len() < 3 {
        return 0.0;
    }

    let n = residuals.len();
    let mean: f64 = residuals.iter().sum::<f64>() / n as f64;

    let mut cov_0 = 0.0f64; // Variance
    let mut cov_1 = 0.0f64; // Lag-1 covariance

    for i in 0..n {
        let xi = residuals[i] - mean;
        cov_0 += xi * xi;
        if i < n - 1 {
            let xi1 = residuals[i + 1] - mean;
            cov_1 += xi * xi1;
        }
    }

    if cov_0 < 1e-20 {
        return 0.0;
    }

    cov_1 / cov_0
}
