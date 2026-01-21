/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::{PulsarDataPoint, PulsarObservation};
use deep_causality_num::RealField;
use std::fs;
use std::io::{self, BufRead};
use std::path::{Path, PathBuf};

/// Get the absolute path to the root data directory (CARGO_MANIFEST_DIR/data)
pub fn get_data_root_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data");
    path
}

/// Get the absolute path to the IPTA data directory (data/ipta/DR2-master/release/VersionB)
pub fn get_ipta_data_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data/ipta/DR2-master/release/VersionB");
    path
}

/// Get NANOGrav data path
pub fn get_nanograv_data_path() -> PathBuf {
    get_data_root_path().join("zenodo/nanograv")
}

/// Recursively loads data points from a .tim file, handling INCLUDE directives.
fn load_tim_data<R>(tim_path: &Path, data_points: &mut Vec<PulsarDataPoint<R>>) -> io::Result<()>
where
    R: RealField + From<f64>,
{
    let file = fs::File::open(tim_path)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Handle INCLUDE directive
        if trimmed.starts_with("INCLUDE") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 2 {
                let include_rel_path = parts[1];
                if let Some(parent) = tim_path.parent() {
                    let include_path = parent.join(include_rel_path);
                    if include_path.exists() {
                        load_tim_data(&include_path, data_points)?;
                    } else {
                        println!("Warning: Included file not found: {:?}", include_path);
                    }
                }
            }
            continue;
        }

        if trimmed.is_empty()
            || trimmed.starts_with("C")
            || trimmed.starts_with("#")
            || trimmed.starts_with("FORMAT")
            || trimmed.starts_with("MODE")
        {
            continue;
        }

        // Parse TEMPO2 format line
        // File Freq MJD Error Site ...
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }

        let freq_res = parts[1].parse::<f64>();
        let mjd_res = parts[2].parse::<f64>();
        let err_res = parts[3].parse::<f64>();

        if let (Ok(freq_f64), Ok(mjd_f64), Ok(err_f64)) = (freq_res, mjd_res, err_res) {
            let unix_seconds = (mjd_f64 - 40587.0) * 86400.0;
            data_points.push(PulsarDataPoint {
                timestamp_unix: unix_seconds as i64,
                mjd: R::from(mjd_f64),
                error_us: R::from(err_f64),
                freq_mhz: R::from(freq_f64),
                residual_s: R::zero(), // Needs detrending later
            });
        }
    }
    Ok(())
}

/// Loads a single pulsar's data from .par and .tim files.
/// Returns None if the pulsar is a Binary (excluded from analysis for now).
pub fn load_pulsar<R>(par_path: &str, tim_path: &str) -> io::Result<Option<PulsarObservation<R>>>
where
    R: RealField + From<f64>,
{
    // 1. Check .par file for BINARY model
    let par_content = fs::read_to_string(par_path)?;
    if par_content.contains("BINARY") {
        return Ok(None); // Skip binary pulsars
    }

    // Extract name from par file
    let mut name = String::new();
    for line in par_content.lines() {
        if line.starts_with("PSRJ") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                name = parts[1].to_string();
            }
        }
    }
    if name.is_empty()
        && let Some(n) = Path::new(par_path).file_stem()
    {
        name = n.to_string_lossy().to_string();
    }

    // 2. Load .tim file (recursively)
    let mut data_points = Vec::new();
    load_tim_data(Path::new(tim_path), &mut data_points)?;

    // Sort by timestamp
    data_points.sort_by(|a, b| a.timestamp_unix.cmp(&b.timestamp_unix));

    if data_points.is_empty() {
        return Ok(None);
    }

    Ok(Some(PulsarObservation {
        name,
        observations: data_points,
    }))
}

/// Recursively find all .tim files in the directory
fn find_tim_files(dir: &Path) -> io::Result<Vec<PathBuf>> {
    let mut tim_files = Vec::new();
    if dir.is_dir() {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                tim_files.extend(find_tim_files(&path)?);
            } else if let Some(name) = path.file_name().and_then(|s| s.to_str()) {
                // IPTA DR2 specific filtering: prefer the main combined tim files
                if name.ends_with(".IPTADR2.tim") {
                    tim_files.push(path);
                }
            }
        }
    }
    Ok(tim_files)
}

pub fn load_all_pta<R>(dir_path: &str) -> io::Result<Vec<PulsarObservation<R>>>
where
    R: RealField + From<f64>,
{
    let mut pulsars = Vec::new();
    let root = Path::new(dir_path);

    // Recursively collect all .tim files
    let mut tim_files = find_tim_files(root)?;

    // Sort for consistent order
    tim_files.sort();

    println!(
        "Found {} IPTA DR2 .tim files in {}",
        tim_files.len(),
        dir_path
    );

    for tim_path in tim_files {
        // Look for corresponding .par file
        // IPTA Pattern: J1234+5678.IPTADR2.tim -> J1234+5678.IPTADR2.par
        let par_path = tim_path.with_extension("par");

        if par_path.exists() {
            println!("Loading PTA Pulsar: {:?}", tim_path.file_name().unwrap());
            match load_pulsar(par_path.to_str().unwrap(), tim_path.to_str().unwrap()) {
                Ok(Some(p)) => {
                    println!("  -> Loaded {} ({} obs)", p.name, p.observations.len());
                    pulsars.push(p);
                }
                Ok(None) => {
                    println!("  -> Skipped (Binary or Empty)");
                }
                Err(e) => {
                    println!("  -> Error: {}", e);
                }
            }
        } else {
            println!(
                "  -> Skipped (No .par file found for {:?})",
                tim_path.file_name()
            );
        }
    }

    Ok(pulsars)
}
