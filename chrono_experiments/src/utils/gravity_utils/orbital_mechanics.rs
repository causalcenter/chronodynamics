/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Orbital Mechanics for Binary Pulsars
//!
//! Keplerian orbit calculations for computing orbital phase, separation,
//! and gravitational field strength.

use deep_causality_physics::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT};
use std::f64::consts::PI;

const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
const MSUN: f64 = 1.98898e30; // kg (solar mass)

/// Default companion mass for systems without measured mass
/// (Typical white dwarf: 0.5 solar masses)
pub const DEFAULT_COMPANION_MASS: f64 = 0.5;

/// Default pulsar mass (typical neutron star)
pub const DEFAULT_PULSAR_MASS: f64 = 1.4;

/// Compute gravitational field strength Φ/c²
///
/// Φ/c² = GM / (r c²)
///
/// Arguments:
/// - separation_m: Pulsar-companion separation (meters)
/// - companion_mass_msun: Companion mass (solar masses)
///
/// Returns: Dimensionless field strength
pub fn compute_field_strength(separation_m: f64, companion_mass_msun: f64) -> f64 {
    let m_kg = companion_mass_msun * MSUN;
    let phi = G * m_kg / separation_m;
    phi / (SPEED_OF_LIGHT * SPEED_OF_LIGHT)
}

/// Compute orbital phase from MJD
///
/// Returns phase in range [0, 1] where:
/// - 0.0 = periastron (closest approach)
/// - 0.5 = apastron (farthest point)
pub fn compute_orbital_phase(mjd: f64, t0_mjd: f64, pb_days: f64) -> f64 {
    let phase = ((mjd - t0_mjd) / pb_days) % 1.0;
    if phase < 0.0 { phase + 1.0 } else { phase }
}

/// Compute pulsar-companion separation
///
/// Uses Keplerian orbit formulation:
/// r = a(1 - e²) / (1 + e cos(true_anomaly))
///
/// Arguments:
/// - phase: Orbital phase [0,1]
/// - a1_lt_s: Semi-major axis in light-seconds
/// - ecc: Eccentricity
/// - om_deg: Longitude of periastron (degrees)
///
/// Returns: Separation in meters
pub fn compute_separation(phase: f64, a1_lt_s: f64, ecc: f64, om_deg: f64) -> f64 {
    // Convert phase to mean anomaly
    let mean_anomaly = 2.0 * PI * phase;

    // Solve Kepler's equation for eccentric anomaly
    let ecc_anomaly = solve_kepler_equation(mean_anomaly, ecc);

    // Compute true anomaly
    let true_anomaly =
        2.0 * ((1.0 + ecc).sqrt() * (ecc_anomaly / 2.0).tan() / ((1.0 - ecc).sqrt())).atan();

    // Add longitude of periastron
    let total_angle = true_anomaly + om_deg.to_radians();

    // Compute separation using orbit equation
    let a_meters = a1_lt_s * SPEED_OF_LIGHT;

    a_meters * (1.0 - ecc * ecc) / (1.0 + ecc * total_angle.cos())
}

/// Solve Kepler's equation: E - e sin(E) = M
///
/// Uses Newton-Raphson iteration to find eccentric anomaly E
/// given mean anomaly M and eccentricity e.
pub fn solve_kepler_equation(mean_anom: f64, ecc: f64) -> f64 {
    let mut e_anom = mean_anom; // Initial guess

    // Newton-Raphson: E_new = E - f(E)/f'(E)
    // where f(E) = E - e sin(E) - M
    // and f'(E) = 1 - e cos(E)
    for _ in 0..10 {
        let f = e_anom - ecc * e_anom.sin() - mean_anom;
        let fp = 1.0 - ecc * e_anom.cos();
        e_anom -= f / fp;

        if f.abs() < 1e-10 {
            break;
        }
    }

    e_anom
}
