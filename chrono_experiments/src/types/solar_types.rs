/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

#[derive(Debug, Clone)]
pub struct SolarVectorResult {
    pub timestamp: i64,
    pub gradient_vector: [f64; 3], // Measured G-vector (pointing TO source)
    pub optical_sun_vector: [f64; 3], // Expected Optical
    pub true_sun_vector: [f64; 3], // Expected True
    pub angle_optical_deg: f64,
    pub angle_true_deg: f64,
    pub scalar_magnitude_ns: f64,
}

impl SolarVectorResult {
    pub fn new(
        timestamp: i64,
        gradient_vector: [f64; 3],
        optical_sun_vector: [f64; 3],
        true_sun_vector: [f64; 3],
        angle_optical_deg: f64,
        angle_true_deg: f64,
        scalar_magnitude_ns: f64,
    ) -> Self {
        Self {
            timestamp,
            gradient_vector,
            optical_sun_vector,
            true_sun_vector,
            angle_optical_deg,
            angle_true_deg,
            scalar_magnitude_ns,
        }
    }
}
