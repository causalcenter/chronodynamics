/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Minimal Linear Algebra module for 7-satellite constellation analysis.
//! Implements 7x7 Matrix operations, Covariance, and Jacobi Eigenvalue solver.
//!   Ported with R: RealField generics.

use deep_causality_num::RealField;

/// 5x5 Matrix type (for 5-satellite covariance)
pub type Matrix<R> = [[R; 5]; 5];

/// 5-element Vector type (for 5-satellite residuals)
pub type Vector<R> = [R; 5];

/// Dynamic Matrix type (Vec of Rows)
pub type DynamicMatrix<R> = Vec<Vec<R>>;

/// Dynamic Vector type
pub type DynamicVector<R> = Vec<R>;

/// Compute the covariance matrix of a set of 5-element vectors.
pub fn covariance<R>(data: &[Vector<R>]) -> Matrix<R>
where
    R: RealField + Copy + From<f64>,
{
    let n = data.len();
    if n < 2 {
        return [[R::zero(); 5]; 5];
    }

    // 1. Compute Mean Vector
    let mut mean = [R::zero(); 5];
    for row in data {
        for i in 0..5 {
            mean[i] += row[i];
        }
    }
    for x in &mut mean {
        *x /= R::from(n as f64);
    }

    // 2. Compute Covariance
    // Cov(X, Y) = E[(X - E[X])(Y - E[Y])]
    let mut cov = [[R::zero(); 5]; 5];

    for row in data {
        // Center the data
        let mut centered = [R::zero(); 5];
        for i in 0..5 {
            centered[i] = row[i] - mean[i];
        }

        // Outer product accumulator
        for i in 0..5 {
            for j in 0..5 {
                cov[i][j] += centered[i] * centered[j];
            }
        }
    }

    // Normalize by (N-1) for sample covariance
    for row in &mut cov {
        for x in row {
            *x /= R::from((n - 1) as f64);
        }
    }

    cov
}

/// Compute Covariance Matrix for Dynamic Vectors (N-dimensional)
pub fn covariance_dynamic<R>(data: &[DynamicVector<R>]) -> DynamicMatrix<R>
where
    R: RealField + Copy + From<f64>,
{
    let n_samples = data.len();
    if n_samples < 2 {
        return Vec::new(); // Return empty or handle error
    }
    let dim = data[0].len();

    // 1. Compute Mean
    let mut mean = vec![R::zero(); dim];
    for row in data {
        if row.len() != dim {
            continue;
        } // Skip mismatched rows (should not happen if aligned)
        for i in 0..dim {
            mean[i] += row[i];
        }
    }
    for x in &mut mean {
        *x /= R::from(n_samples as f64);
    }

    // 2. Compute Covariance
    let mut cov = vec![vec![R::zero(); dim]; dim];

    for row in data {
        if row.len() != dim {
            continue;
        }
        // Center
        let mut centered = vec![R::zero(); dim];
        for i in 0..dim {
            centered[i] = row[i] - mean[i];
        }

        // Accumulate
        for i in 0..dim {
            for j in 0..dim {
                cov[i][j] += centered[i] * centered[j];
            }
        }
    }

    // Normalize
    for item in cov.iter_mut().take(dim) {
        for val in item.iter_mut().take(dim) {
            *val /= R::from(n_samples as f64);
        }
    }

    cov
}

/// Compute Eigenvalues and Eigenvectors of a symmetric matrix using Jacobi method.
/// Returns (eigenvalues, eigenvectors).
/// Eigenvectors are the columns of the returned matrix.
pub fn jacobi_eigenvalues<R>(matrix: &Matrix<R>) -> (Vector<R>, Matrix<R>)
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = 5;
    let max_iter = 100;
    let tolerance = R::from(1e-10);

    let mut a = *matrix; // Copy of matrix to destroy
    let mut v = identity_matrix::<R>(); // Accumulator for eigenvectors
    let mut d = [R::zero(); 5]; // Eigenvalues (diagonal)

    // Initialize d to diagonal of A
    for i in 0..n {
        d[i] = a[i][i];
    }

    // Jacobi Rotation Loop
    for _iter in 0..max_iter {
        // Find largest off-diagonal element
        let mut max_off_diag = R::zero();
        let mut p = 0;
        let mut q = 0;

        for (i, row) in a.iter().enumerate().take(n) {
            for (j, val) in row.iter().enumerate().skip(i + 1).take(n) {
                let val = val.abs();
                if val > max_off_diag {
                    max_off_diag = val;
                    p = i;
                    q = j;
                }
            }
        }

        if max_off_diag < tolerance {
            break; // Converged
        }

        // Compute rotation angle (p, q)
        let diff = d[q] - d[p];
        let phi = R::from(0.5) * (R::from(2.0) * a[p][q]).atan2(diff);

        let c = phi.cos();
        let s = phi.sin();

        // Update diagonal elements
        let d_p_new = c * c * d[p] + s * s * d[q] - R::from(2.0) * s * c * a[p][q];
        let d_q_new = s * s * d[p] + c * c * d[q] + R::from(2.0) * s * c * a[p][q];
        d[p] = d_p_new;
        d[q] = d_q_new;

        // Zero out the off-diagonal element (conceptually)
        a[p][q] = R::zero();

        // Update off-diagonal elements
        // #[allow(clippy::needless_range_loop)]
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            if i != p && i != q {
                let a_ip = a[i][p];
                let a_iq = a[i][q];
                a[i][p] = c * a_ip - s * a_iq;
                a[i][q] = s * a_ip + c * a_iq;
                // Symmetry
                a[p][i] = a[i][p];
                a[q][i] = a[i][q];
            }
        }

        // Update Eigenvectors matrix V
        // V_new = V * Rotation_matrix
        for row in v.iter_mut().take(n) {
            let v_ip = row[p];
            let v_iq = row[q];
            row[p] = c * v_ip - s * v_iq;
            row[q] = s * v_ip + c * v_iq;
        }
    }

    (d, v)
}

/// Dynamic Version of Jacobi Eigenvalue Solver
pub fn jacobi_eigenvalues_dynamic<R>(
    matrix: &DynamicMatrix<R>,
) -> (DynamicVector<R>, DynamicMatrix<R>)
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = matrix.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    let max_iter = n * n * 20; // Scale iterations with size
    let tolerance = R::from(1e-10);

    let mut a = matrix.clone();
    let mut v = identity_matrix_dynamic(n);
    let mut d = vec![R::zero(); n];

    for i in 0..n {
        d[i] = a[i][i];
    }

    for _iter in 0..max_iter {
        let mut max_off_diag = R::zero();
        let mut p = 0;
        let mut q = 0;

        for (i, row) in a.iter().enumerate().take(n) {
            for (j, val) in row.iter().enumerate().take(n).skip(i + 1) {
                let val = val.abs();
                if val > max_off_diag {
                    max_off_diag = val;
                    p = i;
                    q = j;
                }
            }
        }

        if max_off_diag < tolerance {
            break;
        }

        let diff = d[q] - d[p];
        let phi = R::from(0.5) * (R::from(2.0) * a[p][q]).atan2(diff);
        let c = phi.cos();
        let s = phi.sin();

        let d_p_new = c * c * d[p] + s * s * d[q] - R::from(2.0) * s * c * a[p][q];
        let d_q_new = s * s * d[p] + c * c * d[q] + R::from(2.0) * s * c * a[p][q];
        d[p] = d_p_new;
        d[q] = d_q_new;
        a[p][q] = R::zero();

        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            if i != p && i != q {
                let a_ip = a[i][p];
                let a_iq = a[i][q];
                a[i][p] = c * a_ip - s * a_iq;
                a[i][q] = s * a_ip + c * a_iq;
                a[p][i] = a[i][p];
                a[q][i] = a[i][q];
            }
        }

        for row in v.iter_mut() {
            let v_ip = row[p];
            let v_iq = row[q];
            row[p] = c * v_ip - s * v_iq;
            row[q] = s * v_ip + c * v_iq;
        }
    }
    (d, v)
}

fn identity_matrix<R>() -> Matrix<R>
where
    R: RealField + Copy + From<f64>,
{
    let mut m = [[R::zero(); 5]; 5];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = R::one();
    }
    m
}

pub fn identity_matrix_dynamic<R>(n: usize) -> DynamicMatrix<R>
where
    R: RealField + Copy + From<f64>,
{
    let mut m = vec![vec![R::zero(); n]; n];
    for (i, row) in m.iter_mut().enumerate().take(n) {
        row[i] = R::one();
    }
    m
}

/// Transpose a matrix
pub fn transpose<R>(m: &DynamicMatrix<R>) -> DynamicMatrix<R>
where
    R: RealField + Copy,
{
    if m.is_empty() {
        return Vec::new();
    }
    let rows = m.len();
    let cols = m[0].len();
    let mut res = vec![vec![R::zero(); rows]; cols];
    for i in 0..rows {
        for j in 0..cols {
            res[j][i] = m[i][j];
        }
    }
    res
}

/// Multiply two matrices (A * B)
pub fn matrix_multiply<R>(a: &DynamicMatrix<R>, b: &DynamicMatrix<R>) -> DynamicMatrix<R>
where
    R: RealField + Copy,
{
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let rows_a = a.len();
    let cols_a = a[0].len();
    let rows_b = b.len();
    let cols_b = b[0].len();

    if cols_a != rows_b {
        return Vec::new();
    }

    let mut res = vec![vec![R::zero(); cols_b]; rows_a];
    for i in 0..rows_a {
        for j in 0..cols_b {
            let mut sum = R::zero();
            for k in 0..cols_a {
                sum += a[i][k] * b[k][j];
            }
            res[i][j] = sum;
        }
    }
    res
}

/// Solve linear system Ax = b using Gaussian elimination with pivoting
/// Returns x vector
pub fn solve_linear_system<R>(
    a: &DynamicMatrix<R>,
    b: &DynamicVector<R>,
) -> Option<DynamicVector<R>>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = a.len();
    if n == 0 || b.len() != n {
        return None;
    }

    // Augmented matrix [A | b]
    let mut mat = vec![vec![R::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            mat[i][j] = a[i][j];
        }
        mat[i][n] = b[i];
    }

    for i in 0..n {
        // Pivot
        let mut max = i;
        for j in i + 1..n {
            if mat[j][i].abs() > mat[max][i].abs() {
                max = j;
            }
        }
        mat.swap(i, max);

        // Singular Check
        if mat[i][i].abs() < R::from(1e-12) {
            return None;
        }

        // Elimination
        for j in i + 1..n {
            let factor = mat[j][i] / mat[i][i];

            // Use split_at_mut to get both rows safely
            let (before, after) = mat.split_at_mut(j);
            let row_i = &before[i];
            let row_j = &mut after[0];

            for (vj, vi) in row_j
                .iter_mut()
                .take(n + 1)
                .zip(row_i.iter().take(n + 1))
                .skip(i)
            {
                *vj -= factor * *vi;
            }
        }
    }

    // Back substitution
    let mut x = vec![R::zero(); n];
    for i in (0..n).rev() {
        let mut sum = R::zero();
        for j in i + 1..n {
            sum += mat[i][j] * x[j];
        }
        x[i] = (mat[i][n] - sum) / mat[i][i];
    }

    Some(x)
}

// Structure for PCA result
#[derive(Debug, Clone)]
pub struct PcaResult<R> {
    pub eigenvalues: Vec<R>,
    pub eigenvectors: Vec<Vec<R>>,
    pub variance_explained: Vec<R>,
    pub loadings: Vec<Vec<R>>, // Loadings (Components)
}

/// Calculate PCA of a dataset.
pub fn calculate_pca<R>(data: &[Vec<R>], components: usize) -> PcaResult<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    // data is [Variable][Observation]
    let num_vars = data.len();
    let num_obs = data[0].len();

    // 1. Compute Means
    // 1. Compute Means and Standard Deviations (Population)
    let mut means = vec![R::zero(); num_vars];
    let mut std_devs = vec![R::zero(); num_vars];
    let n_f64 = R::from(num_obs as f64);

    for i in 0..num_vars {
        let sum = data[i].iter().fold(R::zero(), |acc, &x| acc + x);
        let mean = sum / n_f64;
        means[i] = mean;

        let sum_sq_diff = data[i].iter().fold(R::zero(), |acc, &x| {
            let diff = x - mean;
            acc + diff * diff
        });

        // Use Population Standard Deviation (N) to match legacy's calculate_statistics
        std_devs[i] = (sum_sq_diff / n_f64).sqrt();
    }

    // 2. Compute Covariance (Sample) and Normalize to Correlation
    let mut cov_matrix = vec![vec![R::zero(); num_vars]; num_vars];
    let n_minus_1 = R::from((num_obs - 1).max(1) as f64);

    for i in 0..num_vars {
        for j in i..num_vars {
            let mut sum = R::zero();
            for (val_i_raw, val_j_raw) in data[i].iter().zip(data[j].iter()).take(num_obs) {
                let val_i = *val_i_raw - means[i];
                let val_j = *val_j_raw - means[j];
                sum += val_i * val_j;
            }
            let cov = sum / n_minus_1;

            // Normalize by Population StdDevs to match Legacy logic
            // (effectively Correlation Matrix * (N / N-1))
            let denom = std_devs[i] * std_devs[j];
            let val = if denom.abs() > R::from(1e-12) {
                cov / denom
            } else {
                R::zero()
            };

            cov_matrix[i][j] = val;
            cov_matrix[j][i] = val;
        }
    }

    // 2. Eigen Decomposition
    let (eig_vals, eig_vecs) = jacobi_eigenvalues_dynamic(&cov_matrix);

    // 3. Sort by Eigenvalue descending
    // We need to pair them to sort
    let dim = eig_vals.len();
    let mut indices: Vec<usize> = (0..dim).collect();
    indices.sort_by(|&i, &j| {
        eig_vals[j]
            .partial_cmp(&eig_vals[i])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let sorted_vals: Vec<R> = indices.iter().map(|&i| eig_vals[i]).collect();

    let mut sorted_vecs = vec![vec![R::zero(); dim]; dim];
    for r in 0..dim {
        for (new_c, &old_c) in indices.iter().enumerate() {
            sorted_vecs[r][new_c] = eig_vecs[r][old_c];
        }
    }

    // 4. Calculate Variance Explained (Only from top k components - matches legacy power iteration)
    let k = components.min(dim);
    let total_var: R = sorted_vals
        .iter()
        .take(k)
        .fold(R::zero(), |acc, &v| acc + v);
    let var_explained: Vec<R> = sorted_vals
        .iter()
        .take(k)
        .map(|&v| {
            if total_var > R::zero() {
                v / total_var
            } else {
                R::zero()
            }
        })
        .collect();

    // 5. Loadings

    let mut loadings = vec![vec![R::zero(); dim]; k];

    for i in 0..k {
        // PC index
        for j in 0..dim {
            // Variable index
            loadings[i][j] = sorted_vecs[j][i];
        }
    }

    PcaResult {
        eigenvalues: sorted_vals.into_iter().take(k).collect(),
        eigenvectors: sorted_vecs, // Full matrix sorted
        variance_explained: var_explained,
        loadings,
    }
}
