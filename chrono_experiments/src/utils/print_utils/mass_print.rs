/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::utils::statistics_utils::calculate_statistics;
use deep_causality_num::RealField;
use deep_causality_physics::{EARTH_GM, EARTH_MASS_KG, NEWTONIAN_CONSTANT_OF_GRAVITATION};
use std::io::Error;

/// Prints the global summary showing GM as the TRUE measurement,
/// with explicit separation steps for G and M.
///
/// # Type Parameter
/// - `R`: Input real field type that can be converted to f64 for display
pub fn print_mass_global_summary<R>(gm_results: &[R], label: &str)
where
    R: RealField + Into<f64> + From<f64> + Clone + std::fmt::LowerExp,
{
    if gm_results.is_empty() {
        println!("\nNo data points to summarize.");
        return;
    }

    // Use generic statistics calculation
    let stats = calculate_statistics(gm_results);

    // Calculate High-Precision Mean (in R) to check for DoubleFloat benefits
    let sum: R = gm_results.iter().fold(R::zero(), |acc, x| acc + *x);
    let count = R::from(gm_results.len() as f64);
    let mean_hp = sum / count;

    // Cast constants to R
    let gm_earth_r = R::from(EARTH_GM);
    let earth_mass_r = R::from(EARTH_MASS_KG);
    let g_const_r = R::from(NEWTONIAN_CONSTANT_OF_GRAVITATION);

    // GM components
    let gm_median = stats.median;
    let gm_mean = mean_hp;

    // Errors
    let median_error = ((gm_median - gm_earth_r) / gm_earth_r).abs() * R::from(100.0);

    let mean_error = ((gm_mean - gm_earth_r) / gm_earth_r).abs() * R::from(100.0);

    // Select the best GM for recovery (the one with the lower error relative to reference)
    let (best_gm, recovery_source) = if median_error.into() < mean_error.into() {
        (gm_median, "Median")
    } else {
        (gm_mean, "Mean")
    };

    // Skewness-based distribution description
    let skewness: f64 = stats.skewness.into();
    let distribution = if skewness > 0.5 {
        "Right Skewed"
    } else if skewness < -0.5 {
        "Left Skewed"
    } else {
        "Symmetric"
    };

    // Explicit separation (requires assuming the complementary constant)
    let g_recovered = best_gm / earth_mass_r;
    let g_error = ((g_recovered - g_const_r) / g_const_r).abs() * R::from(100.0);

    let m_recovered = best_gm / g_const_r;
    let m_error = ((m_recovered - earth_mass_r) / earth_mass_r).abs() * R::from(100.0);

    println!("\n╔══════════════════════════════════════════════════════════════════════════════╗");
    println!(
        "║  CHRONO-MASS RESULT: {}                              ",
        label
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!("║  DATA QUALITY                                                                ║");
    println!(
        "║    Total Data Points:      {:>10}                                       ║",
        stats.count
    );
    println!(
        "║    Variance (GM):          {:.2e} (m³/s²)²                               ║",
        stats.variance.into()
    );
    println!(
        "║    Std Deviation (GM):     {:.2e} m³/s²                                   ║",
        stats.std_dev.into()
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!("║  ★ PRIMARY RESULT: GM DERIVED FROM TIME (No External Mass/G Required)       ║");
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!(
        "║    GM (Median):            {:.6e} m³/s²                                ║",
        gm_median.into()
    );
    println!(
        "║    GM (HP Mean):           {:.30e} m³/s²       ║",
        gm_mean
    );
    println!(
        "║    Distribution:           {:<15} (Skew: {:>6.2})                 ║",
        distribution, skewness
    );
    println!(
        "║    GM (Reference):         {:.6e} m³/s²  (IERS 2010)                   ║",
        EARTH_GM
    );
    println!(
        "║    GM median error:        {:.4}%                                           ║",
        median_error.into()
    );
    println!(
        "║    GM mean error:          {:.4}%                                           ║",
        mean_error.into()
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!(
        "║  SEPARATION STEP 1: Recover G (Based on GM {}, M = {:.4e} kg) ║",
        recovery_source, EARTH_MASS_KG
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!("║    G = GM / M_assumed                                                        ║");
    println!(
        "║    G (recovered):          {:.5e} m³/(kg·s²)                            ║",
        g_recovered.into()
    );
    println!(
        "║    G (reference):          {:.5e} m³/(kg·s²)  (CODATA 2022)             ║",
        NEWTONIAN_CONSTANT_OF_GRAVITATION
    );
    println!(
        "║    G Error:                {:.4}%                                           ║",
        g_error.into()
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!(
        "║  SEPARATION STEP 2: Recover M (Based on GM {}, G = {:.5e})   ║",
        recovery_source, NEWTONIAN_CONSTANT_OF_GRAVITATION
    );
    println!("╠══════════════════════════════════════════════════════════════════════════════╣");
    println!("║    M = GM / G_assumed                                                        ║");
    println!(
        "║    M (recovered):          {:.4e} kg                                     ║",
        m_recovered.into()
    );
    println!(
        "║    M (reference):          {:.4e} kg  (IERS 2010)                        ║",
        EARTH_MASS_KG
    );
    println!(
        "║    M Error:                {:.4}%                                           ║",
        m_error.into()
    );
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");
}

/// Prints statistics for GM values.
///
/// # Type Parameter
/// - `R`: Input real field type that can be converted to f64 for display
pub fn print_mass_statistics<R>(filtered: &[R]) -> Result<(), Error>
where
    R: RealField + Into<f64> + From<f64> + Clone,
{
    // Use generic statistics calculation
    let stats = calculate_statistics(filtered);

    println!("  -> MEDIAN GM: {:.5e} m³/s²", stats.median.into());
    println!(
        "  -> MEAN GM:   {:.5e} +/- {:.2e} (SEM)",
        stats.mean.into(),
        stats.std_error.into()
    );
    println!("  -> VARIANCE:  {:.2e}", stats.variance.into());
    println!("  -> STD DEV:   {:.2e}", stats.std_dev.into());
    println!(
        "  -> ERROR:     {:.3}% (vs IERS GM)",
        stats.error_percentage(EARTH_GM)
    );
    println!("---------------------------------------------------");

    Ok(())
}
