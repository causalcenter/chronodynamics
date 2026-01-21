/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use deep_causality_num::Float;

/// Apply MAD (Median Absolute Deviation) filtering to remove outliers.
///
/// Generic over any `Float` type to maintain precision (f64, DoubleFloat, etc.).
///
/// # Statistical Background
///
/// The constant 1.4826 ≈ 1/Φ⁻¹(3/4) where Φ⁻¹ is the inverse normal CDF.
/// This makes MAD a consistent estimator of σ for Gaussian-distributed data:
/// σ ≈ 1.4826 × MAD
pub fn apply_mad_filter<T: Float>(data: &[T], outlier_sigma: T) -> Vec<T> {
    if data.is_empty() {
        return Vec::new();
    }

    // Sort for median calculation (preserves original order in `data` for final filter)
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // True median: average of middle two elements for even-length arrays
    let median = if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2]
    } else {
        let mid = sorted.len() / 2;
        let two = T::from(2.0).unwrap_or_else(|| T::one() + T::one());
        (sorted[mid - 1] + sorted[mid]) / two
    };

    // Calculate MAD (Median Absolute Deviation)
    let residuals: Vec<T> = sorted.iter().map(|x| (*x - median).abs()).collect();
    let mut sorted_residuals = residuals;
    sorted_residuals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // True median of residuals
    let mad = if sorted_residuals.len() % 2 == 1 {
        sorted_residuals[sorted_residuals.len() / 2]
    } else {
        let mid = sorted_residuals.len() / 2;
        let two = T::from(2.0).unwrap_or_else(|| T::one() + T::one());
        (sorted_residuals[mid - 1] + sorted_residuals[mid]) / two
    };

    // Estimate sigma from MAD (σ ≈ 1.4826 × MAD)
    // Fallback chain: 1.4826 → 1.5 → 1.0
    let scale_factor = T::from(1.4826)
        .or_else(|| T::from(1.5))
        .unwrap_or_else(T::one);
    let sigma_est = scale_factor * mad;

    if sigma_est > T::zero() {
        data.iter()
            .copied()
            .filter(|&g| (g - median).abs() <= outlier_sigma * sigma_est)
            .collect()
    } else {
        data.to_vec()
    }
}
