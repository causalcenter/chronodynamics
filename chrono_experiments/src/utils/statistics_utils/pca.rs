use crate::statistics_utils::PCAResult;
use deep_causality_num::RealField;

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
