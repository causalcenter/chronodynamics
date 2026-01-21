/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use deep_causality_num::RealField;
use deep_causality_physics::{EARTH_MASS_KG, NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT};

/// Results from Lunar Inversion
#[derive(Debug, Clone, Copy)]
pub struct LunarParameters<T> {
    pub mass_kg: T,
    pub distance_m: T,
    pub mass_error_pct: T,
    pub distance_error_pct: T,
}

/// Solve for Moon Mass and Distance using the 27-day clock signal.
///
/// Physics:
/// 1. Period (T) gives Distance (d) via Kepler's 3rd Law: d = cbrt(GM_earth * T^2 / 4pi^2)
/// 2. Amplitude (A) gives Mass (M_moon) via Tidal Potential: A = (GM_moon * r_sat / (c^2 * d^2)) * (T / 2pi) ?
///    Actually, let's use the potential difference formula derived in plan:
///    Amplitude (ns) approx = [ (GM_moon * r_sat) / (c^2 * d^2) ] * (T / 2pi) * 1e9 ??
///    
///    Wait, let's re-derive carefully.
///    Potenital at satellite due to Moon: U = -GM_m / |d - r|
///    Expansion: U approx -GM_m/d - (GM_m/d^2)*r*cos(theta) ...
///    The first term is constant (offset).
///    The second term varies with orbital angle theta (period T).
///    Delta U = (GM_m * r_sat / d^2) * cos(theta).
///    Clock Rate Shift: df/f = Delta U / c^2.
///    Time Bias Accumulation: Integral(df/f dt).
///    Int(cos(2pi*t/T)) = (T/2pi) * sin(...).
///    So Amplitude (Peak) of Bias = (Peak Rate Shift) * (T / 2pi).
///    Peak Rate Shift = (GM_m * r_sat) / (c^2 * d^2).
///    Amplitude A (sec) = [ (GM_m * r_sat) / (c^2 * d^2) ] * (T / 2pi).
///
///    So: A_ns * 1e-9 = [ (G * M_m * r) / (c^2 * d^2) ] * (T / 2pi).
///    Solve for M_m:
///    M_m = (A_ns * 1e-9 * c^2 * d^2 * 2pi) / (G * r * T).
///
pub fn solve_lunar_parameters<T>(period_s: T, amplitude_ns: T, r_sat_m: T) -> LunarParameters<T>
where
    T: RealField + From<f64> + Copy + deep_causality_num::Float,
{
    let g = <T as From<f64>>::from(NEWTONIAN_CONSTANT_OF_GRAVITATION);
    let c = <T as From<f64>>::from(SPEED_OF_LIGHT);
    let m_earth = <T as From<f64>>::from(EARTH_MASS_KG);
    let pi = <T as From<f64>>::from(std::f64::consts::PI);
    let four = <T as From<f64>>::from(4.0);
    let two = <T as From<f64>>::from(2.0);

    // 1. Solve for Distance (d) using Period (Kepler's 3rd Law for Moon orbiting Earth)
    // d = (G * M_earth * T^2 / 4pi^2)^(1/3)
    let d_cubed = (g * m_earth * period_s.powi(2)) / (four * pi.powi(2));
    let distance_m = d_cubed.cbrt();

    // 2. Solve for Mass (M_moon) using Tidal Amplitude
    let threshold = <T as From<f64>>::from(1e-6);
    let mass_kg = if amplitude_ns > threshold {
        (amplitude_ns * <T as From<f64>>::from(1e-9) * c.powi(2) * distance_m.powi(3) * two * pi)
            / (g * r_sat_m.powi(2) * period_s)
    } else {
        T::zero()
    };

    // Reference Values (for error calc)
    let ref_mass = <T as From<f64>>::from(7.342e22); // Moon Mass
    let ref_dist = <T as From<f64>>::from(3.844e8); // Semi-major axis
    let hundred = <T as From<f64>>::from(100.0);

    LunarParameters {
        mass_kg,
        distance_m,
        mass_error_pct: RealField::abs((mass_kg - ref_mass) / ref_mass) * hundred,
        distance_error_pct: RealField::abs((distance_m - ref_dist) / ref_dist) * hundred,
    }
}
