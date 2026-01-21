/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */
use deep_causality_num::RealField;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Loads GM (geocentric gravitational constant) values from a CSV file.
pub fn load_gm_values_from_csv<R>(path: &Path) -> io::Result<Vec<R>>
where
    R: RealField + From<f64>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut gm_values = Vec::new();

    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        // Skip header
        if idx == 0 {
            continue;
        }

        // CSV format: Timestamp,R_A,R_B,Delta_H_km,GM,Error_Pct,Term_Time,Term_Kinetic,Term_Potential
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 5
            && let Ok(gm_f64) = parts[4].parse::<f64>()
        {
            let gm: R = R::from(gm_f64);
            gm_values.push(gm);
        }
    }

    Ok(gm_values)
}
