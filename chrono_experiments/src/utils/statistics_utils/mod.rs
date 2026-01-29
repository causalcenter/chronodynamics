//! # Generic Statistics Utilities
//!
//! Statistics infrastructure for GQCD experiments, generic over `R: RealField`.

use deep_causality_num::RealField;

pub mod error_correlation;
mod pca;

pub use error_correlation::*;
pub use pca::calculate_pca;

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
