use crate::types::gauge_types::GaugeValidationResult;
use deep_causality_num::{FromPrimitive, RealField, ToPrimitive};
use deep_causality_physics::{EARTH_GM, EARTH_MASS_KG, EARTH_ROTATION_RATE};
use std::io;

/// Prints the experiment header.
pub fn print_gauge_header() {
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║               E02 CHRONO-GAUGE VALIDATION EXPERIMENT                         ║");
    println!("║               Lattice Gauge Theory Cross-Validation                          ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝\n");

    println!(
        "🔬 Float Type: {}",
        std::any::type_name::<deep_causality_num::Float106>()
    );
}

/// Prints summary statistics for a set of validation results.
pub fn print_gauge_summary<R: RealField + Clone + Into<f64> + std::fmt::Display>(
    results: &[GaugeValidationResult<R>],
    label: &str,
) -> io::Result<()> {
    if results.is_empty() {
        println!("  No results to summarize for {}", label);
        return Ok(());
    }

    let total_epochs: usize = results.iter().map(|r| r.num_epochs).sum();
    let total_sites: usize = results.iter().map(|r| r.num_lattice_sites).sum();

    // Sum observables
    let mut mass_sum = 0.0f64;
    let mut curl_sum = 0.0f64;
    let mut grad_sum = 0.0f64;
    let mut vort_z_sum = 0.0f64;
    let mut corr_sum = 0.0f64;
    let mut virial_sum = 0.0f64;
    let mut tolman_sum = 0.0f64;
    let mut rot_sum = 0.0f64;

    for r in results {
        mass_sum += r.mean_mass_density.into();
        curl_sum += r.mean_curl_magnitude.into();
        grad_sum += r.mean_gradient_magnitude.into();
        vort_z_sum += r.mean_vorticity_z.into();
        corr_sum += r.time_velocity_correlation.into();
        virial_sum += r.virial_ratio.into();
        tolman_sum += r.tolman_consistency.into();
        rot_sum += r.derived_earth_rotation.into();
    }

    let n = results.len() as f64;

    println!("\n───────────────────────────────────────────────────");
    println!("  {} SUMMARY", label);
    println!("───────────────────────────────────────────────────");
    println!("  Datasets: {:>10}", results.len());
    println!("  Epochs:   {:>10}", total_epochs);
    println!("  Sites:    {:>10}", total_sites);
    println!();
    println!("  📊 OBSERVABLES (mean across datasets):");
    println!("  Mass Density:      {:.6e}", mass_sum / n);
    println!("  Curl Magnitude:    {:.6e}", curl_sum / n);
    println!("  Gradient Mag:      {:.6e}", grad_sum / n);
    println!("  Vorticity-Z:       {:.6e}", vort_z_sum / n);
    println!("  Derived Rotation:  {:.6e}", rot_sum / n);
    println!();
    println!("  📐 PHYSICS VALIDATIONS:");
    println!("  Time-Velocity Corr: {:.6}", corr_sum / n);
    println!("  Virial Ratio:       {:.6}", virial_sum / n);
    println!("  Tolman Consistency: {:.6}", tolman_sum / n);

    Ok(())
}

/// Prints cross-validation table with an optional global correlation override.
pub fn print_cross_validation_table<
    R: RealField + Clone + Copy + std::fmt::Display + std::fmt::LowerExp + FromPrimitive + ToPrimitive,
>(
    results: &[GaugeValidationResult<R>],
    global_corr: Option<R>,
) {
    if results.is_empty() {
        return;
    }

    // Aggregate means
    let n = R::from_usize(results.len()).unwrap();
    let zero = R::zero();
    let one = R::one();

    let sum_accumulator = |acc: R, x: R| acc + x;

    let mass_mean: R = results
        .iter()
        .map(|r| r.mean_mass_density)
        .fold(zero, sum_accumulator)
        / n;

    let curl_mean: R = results
        .iter()
        .map(|r| r.mean_curl_magnitude)
        .fold(zero, sum_accumulator)
        / n;

    let vort_z_mean: R = results
        .iter()
        .map(|r| r.mean_vorticity_z)
        .fold(zero, sum_accumulator)
        / n;

    let tolman_mean: R = results
        .iter()
        .map(|r| r.tolman_consistency)
        .fold(zero, sum_accumulator)
        / n;

    let rot_mean: R = results
        .iter()
        .map(|r| r.derived_earth_rotation)
        .fold(zero, sum_accumulator)
        / n;

    // Filter out zero mass results (datasets where E14 wasn't present)
    let mass_results: Vec<_> = results.iter().filter(|r| r.derived_mass > zero).collect();
    let n_mass = R::from_usize(mass_results.len()).unwrap();

    let gm_mean: R = if n_mass > zero {
        mass_results
            .iter()
            .map(|r| r.derived_gm)
            .fold(zero, sum_accumulator)
            / n_mass
    } else {
        zero
    };

    let mass_val: R = if n_mass > zero {
        mass_results
            .iter()
            .map(|r| r.derived_mass)
            .fold(zero, sum_accumulator)
            / n_mass
    } else {
        zero
    };

    // Use global correlation if provided, otherwise average from results
    let corr_mean: R = global_corr.unwrap_or_else(|| {
        results
            .iter()
            .map(|r| r.time_velocity_correlation)
            .fold(zero, sum_accumulator)
            / n
    });

    println!("\n╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║                     CHRONO-GAUGE VALIDATION SUMMARY                          ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝\n");

    println!("🔬 CROSS-VALIDATION WITH REFERENCE VALUES");
    println!("  ┌──────────────────────────┬───────────┬─────────────┬─────────────┐");
    println!("  │ Experiment               │ Reference │ Gauge Value │ Status      │");
    println!("  ├──────────────────────────┼───────────┼─────────────┼─────────────┤");

    // Constants for validation
    let tolerance_small = R::from_f64(1e-6).unwrap();
    let tolerance_corr = R::from_f64(0.95).unwrap();
    let tolerance_pct = R::from_f64(0.01).unwrap();

    // Tidal (Laplacian) - Reference: 0.00
    let mass_status = if mass_mean.abs() < tolerance_small {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Tidal (Laplacian)        │ 0.00      │ {:>11} │ {:11} │",
        format_near_zero(mass_mean),
        mass_status
    );

    // Frame-Dragging (Curl) - Reference: 0.00
    let curl_status = if curl_mean.abs() < tolerance_small {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Frame-Dragging (Curl)    │ 0.00      │ {:>11} │ {:11} │",
        format_near_zero(curl_mean),
        curl_status
    );

    // Spin (Vorticity-Z) - Reference: 0.00
    let vort_status = if vort_z_mean.abs() < tolerance_small {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Spin (Vorticity-Z)       │ 0.00      │ {:>11} │ {:11} │",
        format_near_zero(vort_z_mean),
        vort_status
    );

    // Momentum (Correlation) - Reference: 1.00
    let corr_status = if corr_mean > tolerance_corr {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Momentum (Correlation)   │ 1.00      │ {:>11.2} │ {:11} │",
        corr_mean.to_f64().unwrap(),
        corr_status
    );

    // Thermodynamics (Tolman) - Reference: 1.00
    let tolman_status = if (tolman_mean - one).abs() < tolerance_pct {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Thermodynamics (Tolman)  │ 1.00      │ {:>11.2} │ {:11} │",
        tolman_mean.to_f64().unwrap(),
        tolman_status
    );

    // Earth Rotation - Reference: ~7.29e-5
    let ref_rot = R::from_f64(EARTH_ROTATION_RATE).unwrap();
    let rot_diff = (rot_mean - ref_rot).abs() / ref_rot;
    let rot_status = if rot_diff < tolerance_pct {
        "✅ PASS"
    } else {
        "⚠ DEVIATES"
    };
    // Precision format for rot: .2e
    println!(
        "  │ Earth Rotation (rad/s)   │ {:>1.2e}   │ {:>11.2e} │ {:11} │",
        EARTH_ROTATION_RATE, rot_mean, rot_status
    );

    // GM - Reference: 3.986e14
    let ref_gm = R::from_f64(EARTH_GM).unwrap();
    let gm_diff = if gm_mean > zero {
        (gm_mean - ref_gm).abs() / ref_gm
    } else {
        one
    };
    let gm_status = if gm_mean > zero && gm_diff < tolerance_pct {
        "✅ PASS"
    } else if gm_mean == zero {
        "SKIP"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Derived GM (m^3/s^2)     │ {:>1.2e}   │ {:>11.2e} │ {:11} │",
        EARTH_GM, gm_mean, gm_status
    );

    // Mass - Reference: 5.972e24
    let ref_mass = R::from_f64(EARTH_MASS_KG).unwrap();
    let mass_diff = if mass_val > zero {
        (mass_val - ref_mass).abs() / ref_mass
    } else {
        one
    };
    let mass_status = if mass_val > zero && mass_diff < tolerance_pct {
        "✅ PASS"
    } else if mass_val == zero {
        "SKIP"
    } else {
        "⚠ DEVIATES"
    };
    println!(
        "  │ Derived Mass (kg)        │ {:>1.2e}   │ {:>11.2e} │ {:11} │",
        EARTH_MASS_KG, mass_val, mass_status
    );

    println!("  └──────────────────────────┴───────────┴─────────────┴─────────────┘");
}

/// Formats a value as "0.00" if it's essentially zero (< 1e-10), otherwise shows scientific notation.
fn format_near_zero<R: RealField + std::fmt::Display + std::fmt::LowerExp + FromPrimitive>(
    value: R,
) -> String {
    let threshold = R::from_f64(1e-10).unwrap();
    if value.abs() < threshold {
        "0.00".to_string()
    } else {
        format!("{:.2e}", value)
    }
}
