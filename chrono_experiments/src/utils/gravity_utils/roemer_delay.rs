/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Römer Delay Calculations for Binary Pulsars
//!
//! Computes light travel time variations due to orbital motion

/// Compute Römer delay for a binary pulsar
///
/// Römer delay is the light travel time variation as the pulsar orbits:
/// Δt = (a₁ sin i / c) × [sin(ω + ν) + e sin ω]
///
/// where:
/// - a₁: projected semi-major axis (light-seconds)
/// - i: orbital inclination (we use sin i ≈ 1 for edge-on, typical for MSPs)
/// - ω: longitude of periastron
/// - ν: true anomaly
/// - e: eccentricity
///
/// Returns: Römer delay in seconds
pub fn compute_roemer_delay(a1_lt_s: f64, ecc: f64, om_deg: f64, true_anomaly: f64) -> f64 {
    let om_rad = om_deg.to_radians();

    // For MSPs, typically sin i ≈ 1 (edge-on orbits observed via timing)
    let sin_i = 1.0;

    // Römer delay formula

    a1_lt_s * sin_i * ((om_rad + true_anomaly).sin() + ecc * om_rad.sin())
}

/// Compute true anomaly from eccentric anomaly and eccentricity
pub fn eccentric_to_true_anomaly(ecc_anom: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 + ecc).sqrt() * (ecc_anom / 2.0).tan() / (1.0 - ecc).sqrt()).atan()
}
