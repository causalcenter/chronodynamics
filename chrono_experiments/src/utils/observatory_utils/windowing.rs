/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::utils::observatory_utils::linear_algebra::calculate_pca;
use crate::utils::statistics_utils::NullDistribution;
use deep_causality_num::RealField;

#[derive(Debug, Clone)]
pub struct AnomalyWindow<R> {
    pub start_index: usize,
    pub end_index: usize,
    pub pc2_variance: R,
    pub p_value: f64,
    pub tier: AnomalyTier,
    pub dominant_sat: String, // Satellite with max PC2 loading
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnomalyTier {
    Normal,
    Elevated,    // > 90%
    Significant, // > 99%
    Extreme,     // > 99.9%
}

pub struct DetectionResults<R>
where
    R: RealField,
{
    pub baseline_windows: usize,
    pub test_windows: usize,
    pub null_dist: NullDistribution<R>,
    pub anomalies: Vec<AnomalyWindow<R>>,
}

/// Runs a sliding window anomaly scan over the fleet data matrix.
///
/// # Arguments
/// * `matrix` - [Sat][Epoch] matrix (Rows=Sats, Cols=Epochs)
/// * `satellite_names` - Names of M satellites (Rows)
/// * `window_size` - Number of epochs per window
/// * `step_size` - Step between windows
pub fn scan_temporal_anomalies<R>(
    matrix: &[Vec<R>],
    satellite_names: &[String],
    window_size: usize,
    step_size: usize,
) -> DetectionResults<R>
where
    R: RealField + Clone + From<f64> + Into<f64> + PartialOrd,
{
    // 1. Validate Input
    let n_sats = matrix.len();
    if n_sats == 0 {
        return DetectionResults {
            baseline_windows: 0,
            test_windows: 0,
            null_dist: NullDistribution::from_samples(vec![]),
            anomalies: vec![],
        };
    }
    let n_epochs = matrix[0].len();
    if n_epochs < 100 {
        return DetectionResults {
            baseline_windows: 0,
            test_windows: 0,
            null_dist: NullDistribution::from_samples(vec![]),
            anomalies: vec![],
        };
    }

    // =========================================================
    // STEP 1: Build Null Distribution (Truncated Baseline)
    // =========================================================
    // Legacy logic: Baseline ends at 40% of total epochs.
    // Windows in baseline are TRUNCATED at the baseline boundary.
    let baseline_epoch_limit = (n_epochs as f64 * 0.40).floor() as usize;
    let mut null_samples: Vec<R> = Vec::new();
    let mut baseline_count = 0;

    let mut start = 0;
    while start < baseline_epoch_limit {
        let end = (start + window_size).min(baseline_epoch_limit);
        if end - start >= 100 {
            // Extract window data
            let mut window_data: Vec<Vec<R>> = Vec::with_capacity(n_sats);
            for row in matrix.iter().take(n_sats) {
                window_data.push(row[start..end].to_vec());
            }

            // PCA
            let pca = calculate_pca(&window_data, 3);
            if pca.variance_explained.len() > 1 {
                null_samples.push(pca.variance_explained[1]);
                baseline_count += 1;
            }
        }
        start += step_size;
    }

    let null_dist = NullDistribution::from_samples(null_samples);

    // Calculate thresholds for Tiering (Interpolated)
    let thresh_p90 = null_dist.percentile(0.90);
    let thresh_p99 = null_dist.percentile(0.99);
    let thresh_p999 = null_dist.percentile(0.999);

    // =========================================================
    // 2. Full Scan (Sliding Window)
    // =========================================================
    let mut anomalies = Vec::new();
    let mut test_count = 0;

    start = 0;
    while start < n_epochs {
        let end = (start + window_size).min(n_epochs);
        if end - start >= 100 {
            test_count += 1;

            // Extract window data
            let mut window_data: Vec<Vec<R>> = Vec::with_capacity(n_sats);
            for row in matrix.iter().take(n_sats) {
                window_data.push(row[start..end].to_vec());
            }

            let pca = calculate_pca(&window_data, 3);
            let pc2_var = if pca.variance_explained.len() > 1 {
                pca.variance_explained[1]
            } else {
                R::zero()
            };

            let dominant_sat = if pca.loadings.len() > 1 && !pca.loadings[1].is_empty() {
                find_dominant_satellite(&pca.loadings[1], satellite_names)
            } else {
                "???".to_string()
            };

            // Check against Null Dist (P-Value for display)
            let p_val = null_dist.p_value(&pc2_var);

            // Use Interpolated Thresholds for Classification (Matches Legacy logic)
            let tier = if pc2_var > thresh_p999 {
                AnomalyTier::Extreme
            } else if pc2_var > thresh_p99 {
                AnomalyTier::Significant
            } else if pc2_var > thresh_p90 {
                AnomalyTier::Elevated
            } else {
                AnomalyTier::Normal
            };

            // Legacy stores ALL events? No, only anomalies.
            // But my main.rs expects anomalies to be filtered?
            // "DETECTED ANOMALIES" table usually lists them.
            // But "Total Windows Tested" counts all.

            if tier != AnomalyTier::Normal {
                anomalies.push(AnomalyWindow {
                    start_index: start,
                    end_index: end,
                    pc2_variance: pc2_var,
                    p_value: p_val,
                    tier,
                    dominant_sat,
                });
            }
        }
        start += step_size;
    }

    DetectionResults {
        baseline_windows: baseline_count,
        test_windows: test_count,
        null_dist,
        anomalies,
    }
}

/// Find satellites with significant loadings (>= 75% of max loading).
fn find_dominant_satellite<R>(loadings: &[R], satellite_names: &[String]) -> String
where
    R: RealField + Clone + Into<f64>,
{
    if loadings.len() != satellite_names.len() || loadings.is_empty() {
        return "???".to_string();
    }

    // 1. Find Max Loading
    let mut max_abs_loading: f64 = 0.0;
    for loading in loadings {
        let val: f64 = (*loading).into();
        if val.abs() > max_abs_loading {
            max_abs_loading = val.abs();
        }
    }

    if max_abs_loading < 1e-9 {
        return "None".to_string();
    }

    // 2. Collect all satellites within 75% of max loading
    let threshold = max_abs_loading * 0.75;
    let mut contributors: Vec<String> = Vec::new();

    for (idx, loading) in loadings.iter().enumerate() {
        let val: f64 = (*loading).into();
        if val.abs() >= threshold {
            contributors.push(satellite_names[idx].clone());
        }
    }

    contributors.join("/")
}
