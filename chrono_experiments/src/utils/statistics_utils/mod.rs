/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! # Generic Statistics Utilities
//!
//! Statistics infrastructure for GQCD experiments, generic over `R: RealField`.

use deep_causality_num::RealField;

/// Generic statistics results for a dataset.
#[derive(Debug, Clone)]
pub struct Statistics<R>
where
    R: RealField,
{
    pub count: usize,
    pub mean: R,
    pub median: R,
    pub variance: R,
    pub std_dev: R,
    pub std_error: R,
    pub skewness: R,
}

impl<R> Statistics<R>
where
    R: RealField + Clone + From<f64>,
{
    /// Creates a new Statistics struct with all fields set to zero.
    pub fn empty() -> Self {
        Self {
            count: 0,
            mean: R::zero(),
            median: R::zero(),
            variance: R::zero(),
            std_dev: R::zero(),
            std_error: R::zero(),
            skewness: R::zero(),
        }
    }

    /// Computes the error percentage relative to a reference value.
    pub fn error_percentage(&self, reference: f64) -> f64
    where
        R: Into<f64>,
    {
        let mean_f64: f64 = self.mean.into();
        if reference.abs() > 1e-15 {
            ((mean_f64 - reference) / reference).abs() * 100.0
        } else {
            0.0
        }
    }
}

/// PCA result containing eigenvalues, eigenvectors, and variance explained.
#[derive(Debug, Clone)]
pub struct PCAResult<R>
where
    R: RealField,
{
    pub eigenvalues: Vec<R>,
    pub eigenvectors: Vec<Vec<R>>, // Loadings (Components)
    pub variance_explained: Vec<R>,
    pub total_variance: R,
}

/// Null distribution for anomaly detection.
#[derive(Debug, Clone)]
pub struct NullDistribution<R>
where
    R: RealField,
{
    pub count: usize,
    pub mean: R,
    pub std_dev: R,
    sorted_samples: Vec<R>,
}

impl<R> NullDistribution<R>
where
    R: RealField + Clone + From<f64> + Into<f64>,
{
    /// Create from samples.
    pub fn from_samples(mut samples: Vec<R>) -> Self {
        if samples.is_empty() {
            return Self {
                count: 0,
                mean: R::zero(),
                std_dev: R::zero(),
                sorted_samples: Vec::new(),
            };
        }

        // Calculate mean
        let n_f64 = samples.len() as f64;
        let sum: R = samples.iter().cloned().fold(R::zero(), |a, b| a + b);
        let mean = sum / R::from(n_f64);

        // Calculate std_dev
        let var_sum: R = samples
            .iter()
            .cloned()
            .map(|x| {
                let d = x - mean;
                d * d
            })
            .fold(R::zero(), |a, b| a + b);
        let std_dev = (var_sum / R::from(n_f64)).sqrt();

        // Sort for percentile calculations
        samples.sort_by(|a, b| {
            let a_f64: f64 = (*a).into();
            let b_f64: f64 = (*b).into();
            a_f64
                .partial_cmp(&b_f64)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        Self {
            count: samples.len(),
            mean,
            std_dev,
            sorted_samples: samples,
        }
    }

    /// Get the p-value percentile for a test statistic.
    pub fn p_value(&self, test_stat: &R) -> f64 {
        if self.sorted_samples.is_empty() {
            return 1.0;
        }
        let test_f64: f64 = (*test_stat).into();
        let exceeding = self
            .sorted_samples
            .iter()
            .filter(|s| {
                let s_f64: f64 = (*(*s)).into();
                s_f64 >= test_f64
            })
            .count();
        exceeding as f64 / self.count as f64
    }

    /// Calculate the value at a specific percentile (0.0 - 1.0).
    /// e.g. 0.99 returns the value greater than 99% of samples.
    pub fn percentile(&self, p: f64) -> R {
        if self.sorted_samples.is_empty() {
            return R::zero();
        }
        let n = self.sorted_samples.len();
        // Rank from 0 to n-1
        let rank = p * (n - 1) as f64;
        let idx = rank.floor() as usize;
        let frac = rank - rank.floor();

        if idx >= n - 1 {
            return self.sorted_samples[n - 1];
        }

        let a = self.sorted_samples[idx];
        let b = self.sorted_samples[idx + 1];

        // Linear interpolation: a + (b - a) * frac
        // R needs to support arithmetic
        a * R::from(1.0 - frac) + b * R::from(frac)
    }
}

/// Calculates comprehensive statistics from a slice of generic data.
pub fn calculate_statistics<R>(data: &[R]) -> Statistics<R>
where
    R: RealField + Clone + From<f64> + Into<f64>,
{
    if data.is_empty() {
        return Statistics::empty();
    }

    let count = data.len();
    let n = R::from(count as f64);

    // Single pass for mean
    let sum: R = data.iter().cloned().fold(R::zero(), |a, b| a + b);
    let mean = sum / n;

    // Second pass for variance
    let sum_sq_diff: R = data
        .iter()
        .cloned()
        .map(|x| {
            let diff = x - mean;
            diff * diff
        })
        .fold(R::zero(), |a, b| a + b);

    let variance = sum_sq_diff / n;
    let std_dev = variance.sqrt();
    let std_error = std_dev / n.sqrt();

    // Skewness
    let skewness = if std_dev > R::zero() {
        let sum_cubed_diff: R = data
            .iter()
            .cloned()
            .map(|x| {
                let diff = x - mean;
                diff * diff * diff
            })
            .fold(R::zero(), |a, b| a + b);
        (sum_cubed_diff / n) / (std_dev * std_dev * std_dev)
    } else {
        R::zero()
    };

    // Median (requires sorting)
    let mut sorted: Vec<R> = data.to_vec();
    sorted.sort_by(|a, b| {
        let a_f64: f64 = (*a).into();
        let b_f64: f64 = (*b).into();
        a_f64
            .partial_cmp(&b_f64)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let median = if count.is_multiple_of(2) {
        (sorted[count / 2 - 1] + sorted[count / 2]) / R::from(2.0_f64)
    } else {
        sorted[count / 2]
    };

    Statistics {
        count,
        mean,
        median,
        variance,
        std_dev,
        std_error,
        skewness,
    }
}

/// Calculate Pearson correlation between two series.
pub fn calculate_correlation<R>(x: &[R], y: &[R]) -> R
where
    R: RealField + Clone + From<f64>,
{
    if x.len() != y.len() || x.is_empty() {
        return R::zero();
    }

    let n = R::from(x.len() as f64);

    // Means
    let mean_x: R = x.iter().cloned().fold(R::zero(), |a, b| a + b) / n;
    let mean_y: R = y.iter().cloned().fold(R::zero(), |a, b| a + b) / n;

    // Covariance and variances
    let mut cov = R::zero();
    let mut var_x = R::zero();
    let mut var_y = R::zero();

    for (xi, yi) in x.iter().zip(y.iter()) {
        let dx = *xi - mean_x;
        let dy = *yi - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    let denom = (var_x * var_y).sqrt();
    if denom == R::zero() {
        R::zero()
    } else {
        cov / denom
    }
}

/// Calculate PCA using power iteration method.
///
/// # Arguments
/// * `data` - Matrix as Vec of rows, each row is a sample with features as columns
/// * `num_components` - Number of principal components to extract
///
/// # Returns
/// PCA result with eigenvalues, eigenvectors, and variance explained
pub fn calculate_pca<R>(data: &[Vec<R>], num_components: usize) -> PCAResult<R>
where
    R: RealField + Clone + Into<f64> + From<f64>,
{
    if data.is_empty() || data[0].is_empty() {
        return PCAResult {
            eigenvalues: Vec::new(),
            eigenvectors: Vec::new(),
            variance_explained: Vec::new(),
            total_variance: R::zero(),
        };
    }

    let n_samples = data.len();
    let n_features = data[0].len();
    let n_f64 = n_samples as f64;

    // 1. Standardize the data (z-score: mean=0, std=1) - matches legacy
    let mut means = vec![R::zero(); n_features];
    let mut std_devs = vec![R::zero(); n_features];

    // Calculate means
    for row in data {
        for (j, val) in row.iter().enumerate() {
            means[j] += *val;
        }
    }
    for m in &mut means {
        *m /= R::from(n_f64);
    }

    // Calculate std devs
    for row in data {
        for (j, val) in row.iter().enumerate() {
            let diff = *val - means[j];
            std_devs[j] += diff * diff;
        }
    }
    for s in &mut std_devs {
        *s = (*s / R::from(n_f64)).sqrt();
    }

    // Standardize (z-score normalization)
    let standardized: Vec<Vec<R>> = data
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(j, v)| {
                    if std_devs[j] > R::from(1e-10) {
                        (*v - means[j]) / std_devs[j]
                    } else {
                        R::zero()
                    }
                })
                .collect()
        })
        .collect();

    // 2. Compute covariance matrix (n_features x n_features) with sample covariance (n-1)
    let mut cov = vec![vec![R::zero(); n_features]; n_features];
    for row in &standardized {
        for i in 0..n_features {
            for j in 0..n_features {
                cov[i][j] += row[i] * row[j];
            }
        }
    }
    // Use (n-1) for sample covariance - matches legacy
    let divisor = R::from((n_samples - 1).max(1) as f64);
    for row in cov.iter_mut().take(n_features) {
        for val in row.iter_mut().take(n_features) {
            *val /= divisor;
        }
    }

    // 3. Power iteration for eigenvalues
    let mut eigenvalues = Vec::with_capacity(num_components);
    let mut eigenvectors = Vec::with_capacity(num_components);
    let num_to_extract = num_components.min(n_features);

    for _ in 0..num_to_extract {
        // Initialize uniform vector (matches legacy: 1.0 / sqrt(n))
        let init_val = R::from(1.0 / (n_features as f64).sqrt());
        let mut v: Vec<R> = vec![init_val; n_features];

        // Normalize
        let norm: R = v
            .iter()
            .map(|x| *x * *x)
            .fold(R::zero(), |a, b| a + b)
            .sqrt();
        for x in &mut v {
            *x /= norm;
        }

        // Power iteration (100 iterations)
        for _ in 0..100 {
            // Matrix-vector multiply: v_new = cov * v
            let mut v_new = vec![R::zero(); n_features];
            for i in 0..n_features {
                for j in 0..n_features {
                    v_new[i] += cov[i][j] * v[j];
                }
            }

            // Normalize
            let norm: R = v_new
                .iter()
                .map(|x| *x * *x)
                .fold(R::zero(), |a, b| a + b)
                .sqrt();
            if norm == R::zero() {
                break;
            }
            for x in &mut v_new {
                *x /= norm;
            }
            v = v_new;
        }

        // Eigenvalue = v^T * cov * v (Rayleigh quotient)
        let mut cov_v = vec![R::zero(); n_features];
        for i in 0..n_features {
            for j in 0..n_features {
                cov_v[i] += cov[i][j] * v[j];
            }
        }
        let eigenvalue: R = v
            .iter()
            .zip(cov_v.iter())
            .map(|(vi, ci)| *vi * *ci)
            .fold(R::zero(), |a, b| a + b);

        eigenvalues.push(eigenvalue);
        eigenvectors.push(v.clone());

        // Deflate: cov = cov - eigenvalue * v * v^T
        for i in 0..n_features {
            for j in 0..n_features {
                cov[i][j] -= eigenvalue * v[i] * v[j];
            }
        }
    }

    // 4. Calculate total variance and variance explained
    let total_variance: R = eigenvalues.iter().cloned().fold(R::zero(), |a, b| a + b);
    let variance_explained: Vec<R> = if total_variance == R::zero() {
        eigenvalues.iter().map(|_| R::zero()).collect()
    } else {
        eigenvalues.iter().map(|e| *e / total_variance).collect()
    };

    PCAResult {
        eigenvalues,
        eigenvectors,
        variance_explained,
        total_variance,
    }
}
