/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use deep_causality_physics::NEWTONIAN_CONSTANT_OF_GRAVITATION;

/// Calculates the Mass of a planetary body from a derived Gravitational Constant ($G$).
///
/// This assumes the Universality of Free Fall (Equivalence Principle). If we measure
/// a specific $G_{derived}$ from the time dilation field of a planet, we can solve
/// for the Mass $M$ that would generate that field, assuming the standard value
/// of $G_{ref}$ is the true universal constant.
///
/// # Formula
/// $$ M_{derived} = \frac{G_{derived} \times M_{ref}}{G_{ref}} $$
///
/// # Arguments
/// * `derived_g` - The Gravitational Constant calculated from clock drift ($m^3 kg^{-1} s^{-2}$).
/// * `reference_mass` - The known reference mass of the body (kg) used in the initial G derivation.
///
/// # Returns
/// * `f64` - The derived mass (kg).
pub fn calculate_derived_mass(derived_g: f64, reference_mass: f64) -> f64 {
    (derived_g * reference_mass) / NEWTONIAN_CONSTANT_OF_GRAVITATION
}

/// Calculates the percentage error between a derived value and a reference value.
pub fn calculate_error_percentage(derived: f64, reference: f64) -> f64 {
    (derived - reference).abs() / reference * 100.0
}
