/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use deep_causality_physics::{EARTH_GM, EARTH_MASS_KG, NEWTONIAN_CONSTANT_OF_GRAVITATION};
use std::fs::{self, File};
use std::io::{self, Write};

/// Saves global aggregated GM results to a CSV file.
/// Now correctly saves GM (geocentric gravitational constant) with proper error calculations.
pub fn save_global_results(out_path: &str, year: &str, gm_results: &[f64]) -> io::Result<()> {
    if gm_results.is_empty() {
        return Ok(());
    }

    let year_out_path = format!("{}/{}", out_path, year);
    fs::create_dir_all(&year_out_path)?;

    let file_path = format!("{}/{}_global_results_aggregated.csv", year_out_path, year);
    let mut file = File::create(&file_path)?;

    // Calculate statistics for GM values
    let mut sorted = gm_results.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_gm = sorted[sorted.len() / 2];
    let mean_gm: f64 = sorted.iter().sum::<f64>() / sorted.len() as f64;
    let variance: f64 = sorted
        .iter()
        .map(|value| {
            let diff = mean_gm - (*value);
            diff * diff
        })
        .sum::<f64>()
        / sorted.len() as f64;
    let std_dev = variance.sqrt();

    // GM error against IERS 2010 reference
    let gm_error = ((median_gm - EARTH_GM) / EARTH_GM).abs() * 100.0;

    // Separation: Recover G (assuming M)
    let g_recovered = median_gm / EARTH_MASS_KG;
    let g_error = ((g_recovered - NEWTONIAN_CONSTANT_OF_GRAVITATION)
        / NEWTONIAN_CONSTANT_OF_GRAVITATION)
        .abs()
        * 100.0;

    // Separation: Recover M (assuming G)
    let m_recovered = median_gm / NEWTONIAN_CONSTANT_OF_GRAVITATION;
    let m_error = ((m_recovered - EARTH_MASS_KG) / EARTH_MASS_KG).abs() * 100.0;

    // Write header and data
    writeln!(file, "Metric,Value")?;
    writeln!(file, "Year,{}", year)?;
    writeln!(file, "Total_Data_Points,{}", gm_results.len())?;

    // Primary result: GM
    writeln!(file, "Global_Median_GM,{:.10e}", median_gm)?;
    writeln!(file, "Global_Mean_GM,{:.10e}", mean_gm)?;
    writeln!(file, "Std_Deviation_GM,{:.10e}", std_dev)?;
    writeln!(file, "Reference_GM,{:.10e}", EARTH_GM)?;
    writeln!(file, "Error_GM_Pct,{:.6}", gm_error)?;

    // Separation Step 1: G (assuming M)
    writeln!(file, "Recovered_G,{:.10e}", g_recovered)?;
    writeln!(
        file,
        "Reference_G,{:.10e}",
        NEWTONIAN_CONSTANT_OF_GRAVITATION
    )?;
    writeln!(file, "Assumed_M_kg,{:.10e}", EARTH_MASS_KG)?;
    writeln!(file, "Error_G_Pct,{:.6}", g_error)?;

    // Separation Step 2: M (assuming G)
    writeln!(file, "Recovered_M_kg,{:.10e}", m_recovered)?;
    writeln!(file, "Reference_M_kg,{:.10e}", EARTH_MASS_KG)?;
    writeln!(file, "Assumed_G,{:.10e}", NEWTONIAN_CONSTANT_OF_GRAVITATION)?;
    writeln!(file, "Error_M_Pct,{:.6}", m_error)?;

    println!("Saved global results to: {}", file_path);

    Ok(())
}
