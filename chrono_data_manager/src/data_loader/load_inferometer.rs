/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Data loaders for atom interferometry datasets.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use deep_causality_num::RealField;

use crate::types::interferometry_types::{EinsteinElevatorData, RamanShot};

/// Get the base path for Zenodo data
pub fn get_zenodo_data_path() -> std::path::PathBuf {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    std::path::Path::new(manifest_dir)
        .join("data")
        .join("zenodo")
}

/// Load Raman interferometry data from arbitrary orientations dataset
pub fn load_arb_rot_data<R>(base_path: &Path) -> io::Result<Vec<RamanShot<R>>>
where
    R: RealField + From<f64>,
{
    let mut all_shots = Vec::new();

    // Look for data in Figure 3 & 5 subdirectories
    let fig_path = base_path.join("Figure 3 & 5");
    if !fig_path.exists() {
        return Ok(all_shots);
    }

    // Iterate through date directories
    for date_entry in std::fs::read_dir(&fig_path)? {
        let date_entry = date_entry?;
        if !date_entry.file_type()?.is_dir() {
            continue;
        }

        let raman_path = date_entry.path().join("Raman");
        if !raman_path.exists() {
            continue;
        }

        // Iterate through Run directories
        for run_entry in std::fs::read_dir(&raman_path)? {
            let run_entry = run_entry?;
            if !run_entry.file_type()?.is_dir() {
                continue;
            }

            // Look for AvgRatios-kD.txt files
            for file_entry in std::fs::read_dir(run_entry.path())? {
                let file_entry = file_entry?;
                let file_name = file_entry.file_name();
                let file_name_str = file_name.to_string_lossy();

                if file_name_str.contains("AvgRatios-kD.txt") {
                    let shots = parse_raman_file(&file_entry.path())?;
                    all_shots.extend(shots);
                }
            }
        }
    }

    println!(
        "  Loaded {} Raman shots from arb_rot dataset",
        all_shots.len()
    );
    Ok(all_shots)
}

/// Parse a single Raman AvgRatios file
fn parse_raman_file<R>(path: &Path) -> io::Result<Vec<RamanShot<R>>>
where
    R: RealField + From<f64>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut shots = Vec::new();

    for line in reader.lines() {
        let line = line?;

        // Skip header (starts with #)
        if line.starts_with('#') || line.is_empty() {
            continue;
        }

        // Parse tab-separated values
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 8
            && let (Ok(iteration), Ok(phase_f64), Ok(ratio_mean_f64), Ok(ratio_std_dev_f64)) = (
                parts[0].parse::<u32>(),
                parts[4].parse::<f64>(),
                parts[6].parse::<f64>(),
                parts[7].parse::<f64>(),
            )
        {
            let timestamp = format!("{} {}", parts[1], parts[2]);
            let phase: R = R::from(phase_f64);
            let ratio_mean: R = R::from(ratio_mean_f64);
            let ratio_std_dev: R = R::from(ratio_std_dev_f64);
            shots.push(RamanShot::new(
                iteration,
                timestamp,
                phase,
                ratio_mean,
                ratio_std_dev,
            ));
        }
    }

    Ok(shots)
}

/// Load Einstein Elevator data from MATLAB .m file
pub fn load_einstein_elevator_data<R>(path: &Path) -> io::Result<EinsteinElevatorData<R>>
where
    R: RealField + From<f64>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut freq_data_f64 = Vec::new();
    let mut ratio_data_f64 = Vec::new();
    let mut pe_simulated_f64 = Vec::new();
    let mut residual_f64 = Vec::new();

    let mut temperature: f64 = 7.5e-6;

    let mut current_array: Option<&str> = None;
    let mut collecting_array = false;
    let mut array_buffer = String::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Parse simple assignments
        if trimmed.starts_with("Temp =")
            && let Some(val) = extract_matlab_value(trimmed, "Temp =")
        {
            temperature = val;
        }

        // Detect array start
        if trimmed.starts_with("Freq_data =") {
            current_array = Some("freq");
            collecting_array = true;
            array_buffer.clear();
        } else if trimmed.starts_with("Ratio_data =") {
            current_array = Some("ratio");
            collecting_array = true;
            array_buffer.clear();
        } else if trimmed.starts_with("Pe =") && !trimmed.starts_with("PEpoch") {
            current_array = Some("pe");
            collecting_array = true;
            array_buffer.clear();
        } else if trimmed.starts_with("Residual =") {
            current_array = Some("residual");
            collecting_array = true;
            array_buffer.clear();
        }

        // Collect array content
        if collecting_array {
            array_buffer.push_str(&line);
            array_buffer.push(' ');

            // Check for array end (]; or just ])
            if trimmed.contains("];") || (trimmed.ends_with(']') && !trimmed.contains('[')) {
                // Parse the collected array
                let values = parse_matlab_array(&array_buffer);
                match current_array {
                    Some("freq") => freq_data_f64 = values,
                    Some("ratio") => ratio_data_f64 = values,
                    Some("pe") => pe_simulated_f64 = values,
                    Some("residual") => residual_f64 = values,
                    _ => {}
                }
                collecting_array = false;
                current_array = None;
                array_buffer.clear();
            }
        }
    }

    // Einstein Elevator uses 480ms interrogation time (from Figure 6 folder name)
    let interrogation_time: f64 = 0.480;

    // Convert all f64 to R
    let freq_data: Vec<R> = freq_data_f64.into_iter().map(R::from).collect();
    let ratio_data: Vec<R> = ratio_data_f64.into_iter().map(R::from).collect();
    let pe_simulated: Vec<R> = pe_simulated_f64.into_iter().map(R::from).collect();
    let residual: Vec<R> = residual_f64.into_iter().map(R::from).collect();

    Ok(EinsteinElevatorData::new(
        freq_data,
        ratio_data,
        pe_simulated,
        residual,
        R::from(interrogation_time),
        R::from(temperature),
    ))
}

/// Extract a simple numeric value from MATLAB assignment
fn extract_matlab_value(line: &str, prefix: &str) -> Option<f64> {
    let after = line.strip_prefix(prefix)?;
    let cleaned: String = after
        .chars()
        .filter(|c| {
            c.is_ascii_digit() || *c == '.' || *c == '-' || *c == 'E' || *c == 'e' || *c == '+'
        })
        .collect();
    cleaned.parse().ok()
}

/// Parse MATLAB array notation [a b c ...] or [a, b, c, ...]
fn parse_matlab_array(content: &str) -> Vec<f64> {
    let mut values = Vec::new();

    // Remove brackets and continuation markers
    let cleaned = content
        .replace(['[', ']', ';'], " ")
        .replace("...", " ")
        .replace(['\n', '\r'], " ");

    // Split by whitespace or commas
    for token in cleaned.split(|c: char| c.is_whitespace() || c == ',') {
        let token = token.trim();
        if token.is_empty() || token.starts_with('%') || token.starts_with('#') {
            continue;
        }
        if let Ok(val) = token.parse::<f64>() {
            values.push(val);
        }
    }

    values
}
