/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::utils::statistics_utils::{PCAResult, calculate_pca};
use deep_causality_num::RealField;

pub struct FleetCoherenceMetrics<R>
where
    R: RealField,
{
    pub pc1_variance: R,
    pub pc2_variance: R,
    pub quadrupole_power: R, // Differential (PC2+)
    pub monopole_power: R,   // Common mode (PC1)
    pub is_coherent: bool,
    pub pca_result: PCAResult<R>,
}

/// Computes gauge theoretic coherence metrics for the fleet.
pub fn compute_fleet_coherence<R>(matrix: &[Vec<R>]) -> FleetCoherenceMetrics<R>
where
    R: RealField + Clone + From<f64> + Into<f64>,
{
    // Perform PCA on the [N x 7] matrix
    // N epochs, 7 satellites (features)
    let pca = calculate_pca(matrix, 5);

    if pca.variance_explained.len() < 2 {
        return FleetCoherenceMetrics {
            pc1_variance: R::zero(),
            pc2_variance: R::zero(),
            quadrupole_power: R::zero(),
            monopole_power: R::zero(),
            is_coherent: false,
            pca_result: pca,
        };
    }

    let pc1_var = pca.variance_explained[0];
    let pc2_var = pca.variance_explained[1];

    // Interpreting Gauge Theory analogues:
    // Monopole Power (Scalar field) ~ PC1 (Common Mode)
    // Quadrupole Power (Tensor field) ~ Sum of PC2..PC5 (Differential Modes)

    let monopole_power = pca.eigenvalues[0];

    let mut quadrupole_power = R::zero();
    for i in 1..pca.eigenvalues.len() {
        quadrupole_power += pca.eigenvalues[i];
    }

    // Coherence Check:
    // System is "Coherent" if Monopole (Common Mode) dominates > 50% variance
    let threshold: R = R::from(0.5);
    let is_coherent = pc1_var > threshold;

    FleetCoherenceMetrics {
        pc1_variance: pc1_var,
        pc2_variance: pc2_var,
        quadrupole_power,
        monopole_power,
        is_coherent,
        pca_result: pca,
    }
}
