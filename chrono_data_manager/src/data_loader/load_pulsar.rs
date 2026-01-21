/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! NANOGrav .par File Parser
//!
//! Parses tempo2 format .par files containing binary pulsar orbital parameters.
use crate::types::strong_gravity_types::{BinaryPulsarParams, TimingResidual};
use deep_causality_num::RealField;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Raw Time of Arrival measurement from .tim file
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields reserved for future TOA analysis
pub struct RawTOA<R>
where
    R: RealField,
{
    pub mjd: R,
    pub frequency_mhz: R,
    pub error_us: R,
}

/// Parse a tempo2 .par file
///
/// Format: PARAM value [uncertainty] [fit_flag]
/// Example: PB 5.74104238 1 1.81e-10
pub fn parse_par_file<R>(path: &Path) -> io::Result<BinaryPulsarParams<R>>
where
    R: RealField + From<f64> + Into<f64>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut params: BinaryPulsarParams<R> = BinaryPulsarParams::new(String::from("UNKNOWN"));

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }

        match tokens[0] {
            "PSR" | "PSRJ" => {
                if tokens.len() > 1 {
                    params.psrj = tokens[1].to_string();
                }
            }
            "BINARY" => {
                if tokens.len() > 1 {
                    params.binary_model = tokens[1].to_string();
                }
            }
            "PB" => {
                params.pb_days = parse_value_generic(tokens.get(1));
            }
            "ECC" => {
                params.ecc = parse_value_generic(tokens.get(1));
            }
            "A1" => {
                params.a1_lt_s = parse_value_generic(tokens.get(1));
            }
            "T0" => {
                params.t0_mjd = parse_value_generic(tokens.get(1));
            }
            "OM" => {
                params.om_deg = parse_value_generic(tokens.get(1));
            }
            "M2" => {
                params.m_companion = parse_value_generic(tokens.get(1));
            }
            "MTOT" => {
                // If we have total mass and companion mass, can infer pulsar mass
                if let Some(mtot) = parse_value_f64(tokens.get(1))
                    && let Some(m2_ref) = &params.m_companion
                {
                    let m2_f64: f64 = (*m2_ref).into();
                    params.m_pulsar = Some(R::from(mtot - m2_f64));
                }
            }
            "M1" | "MASS1" => {
                params.m_pulsar = parse_value_generic(tokens.get(1));
            }
            "GAMMA" => {
                params.gamma_ms = parse_value_generic(tokens.get(1));
            }
            "PBDOT" => {
                params.pbdot = parse_value_generic(tokens.get(1));
            }
            "EPS1" => {
                params.eps1 = parse_value_generic(tokens.get(1));
            }
            "EPS2" => {
                params.eps2 = parse_value_generic(tokens.get(1));
            }
            "TASC" => {
                params.tasc_mjd = parse_value_generic(tokens.get(1));
            }
            _ => {}
        }
    }

    // NANOGrav MSP binaries typically don't have measured masses
    // Use defaults: 1.4 M☉ pulsar (canonical NS), 0.5 M☉ companion (typical WD)
    if params.m_pulsar.is_none() {
        params.m_pulsar = Some(R::from(1.4));
    }
    if params.m_companion.is_none() {
        params.m_companion = Some(R::from(0.5));
    }

    Ok(params)
}

/// Parse scientific notation value (handles D-format: 1.5D-10) to f64
fn parse_value_f64(token: Option<&&str>) -> Option<f64> {
    token.and_then(|s| {
        // Replace 'D' with 'E' for Fortran-style scientific notation
        let normalized = s.replace('D', "E").replace('d', "e");
        normalized.parse::<f64>().ok()
    })
}

/// Parse scientific notation value to generic R
fn parse_value_generic<R>(token: Option<&&str>) -> Option<R>
where
    R: RealField + From<f64>,
{
    parse_value_f64(token).map(R::from)
}

/// Discover all binary pulsars in NANOGrav dataset
pub fn discover_binary_systems<R>(nanograv_dir: &Path) -> io::Result<Vec<BinaryPulsarParams<R>>>
where
    R: RealField + From<f64> + Into<f64>,
{
    let par_dir = nanograv_dir.join("narrowband/par");

    if !par_dir.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Par directory not found: {}", par_dir.display()),
        ));
    }

    let mut binaries = Vec::new();

    for entry in std::fs::read_dir(par_dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.extension().and_then(|s| s.to_str()) == Some("par") {
            match parse_par_file(&path) {
                Ok(params) => {
                    // Only include if it's a binary with sufficient parameters
                    if params.is_usable() {
                        binaries.push(params);
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Failed to parse {}: {}", path.display(), e);
                }
            }
        }
    }

    Ok(binaries)
}

/// Load timing residuals from NANOGrav residuals file
///
/// Format: MJD  residual_us  uncertainty_us
#[allow(dead_code)] // Reserved for future timing residual analysis
pub fn load_timing_residuals<R>(
    nanograv_dir: &Path,
    pulsar_name: &str,
) -> io::Result<Vec<TimingResidual<R>>>
where
    R: RealField + From<f64>,
{
    // NANOGrav format: residuals/{PSR}_NG15yr_nb.full.res
    let residual_path = nanograv_dir
        .join("residuals")
        .join(format!("{}_NG15yr_nb.full.res", pulsar_name));

    if !residual_path.exists() {
        return Ok(Vec::new()); // No residuals available for this pulsar
    }

    let file = File::open(&residual_path)?;
    let reader = BufReader::new(file);

    let mut residuals = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() >= 3
            && let (Some(mjd_f64), Some(res_us_f64), Some(unc_us_f64)) = (
                tokens[0].parse::<f64>().ok(),
                tokens[1].parse::<f64>().ok(),
                tokens[2].parse::<f64>().ok(),
            )
        {
            residuals.push(TimingResidual {
                mjd: R::from(mjd_f64),
                residual_ns: R::from(res_us_f64 * 1000.0), // Convert μs → ns
                uncertainty_ns: R::from(unc_us_f64 * 1000.0),
            });
        }
    }

    Ok(residuals)
}

/// Parse NANOGrav .tim file to extract raw TOAs
///
/// Format: filename frequency MJD error site [flags...]
/// Example: file.fits 1400.5 58000.123456 0.5 ao -format Tempo2 ...
pub fn parse_tim_file<R>(tim_path: &Path) -> io::Result<Vec<RawTOA<R>>>
where
    R: RealField + From<f64>,
{
    let file = File::open(tim_path)?;
    let reader = BufReader::new(file);

    let mut toas = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Skip FORMAT, comments, empty lines
        if trimmed.is_empty() || trimmed.starts_with("C ") || trimmed.starts_with("FORMAT") {
            continue;
        }

        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() < 4 {
            continue;
        }

        // Parse: filename freq MJD error site ...
        if let (Some(freq_f64), Some(mjd_f64), Some(err_f64)) = (
            tokens[1].parse::<f64>().ok(),
            tokens[2].parse::<f64>().ok(),
            tokens[3].parse::<f64>().ok(),
        ) {
            toas.push(RawTOA {
                mjd: R::from(mjd_f64),
                frequency_mhz: R::from(freq_f64),
                error_us: R::from(err_f64),
            });
        }
    }

    Ok(toas)
}

/// Load raw TOAs for a pulsar from NANOGrav .tim file
pub fn load_raw_toas_data<R>(nanograv_dir: &Path, pulsar_name: &str) -> io::Result<Vec<RawTOA<R>>>
where
    R: RealField + From<f64>,
{
    // NANOGrav format: narrowband/tim/{PSR}_PINT_*.nb.tim
    let tim_pattern = format!("{}_PINT_", pulsar_name);
    let tim_dir = nanograv_dir.join("narrowband/tim");

    if !tim_dir.exists() {
        return Ok(Vec::new());
    }

    // Find matching .tim file
    for entry in std::fs::read_dir(tim_dir)? {
        let entry = entry?;
        let path = entry.path();

        if let Some(name) = path.file_name().and_then(|n| n.to_str())
            && name.starts_with(&tim_pattern)
            && name.ends_with(".tim")
        {
            return parse_tim_file(&path);
        }
    }

    Ok(Vec::new())
}
