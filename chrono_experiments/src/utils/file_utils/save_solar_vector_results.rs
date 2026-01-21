/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::types::solar_types::SolarVectorResult;
use std::fs::File;
use std::io::{self, Write};

pub fn save_solar_vector_results(path: &str, results: &[SolarVectorResult]) -> io::Result<()> {
    let filename = format!("{}/solar_vector_analysis.csv", path);
    let mut file = File::create(&filename)?;

    // Write header
    writeln!(
        file,
        "Timestamp,G_Diff_ns,Angle_Sep_Deg,Opt_Sun_X,Opt_Sun_Y,Opt_Sun_Z,True_Sun_X,True_Sun_Y,True_Sun_Z,Grad_X,Grad_Y,Grad_Z"
    )?;

    // Write data
    for r in results {
        writeln!(
            file,
            "{},{:.6e},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.6e},{:.6e},{:.6e}",
            r.timestamp,
            r.scalar_magnitude_ns,
            r.angle_optical_deg,
            r.optical_sun_vector[0],
            r.optical_sun_vector[1],
            r.optical_sun_vector[2],
            r.true_sun_vector[0],
            r.true_sun_vector[1],
            r.true_sun_vector[2],
            r.gradient_vector[0],
            r.gradient_vector[1],
            r.gradient_vector[2]
        )?;
    }

    println!("  → Saved: {}", filename);
    Ok(())
}
