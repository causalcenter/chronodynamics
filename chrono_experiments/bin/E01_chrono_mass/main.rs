/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */
use chrono_data_manager::{
    ANOMALOUS_WEEKS, DataManager, YEARS, get_gnss_data_input_path, get_year_datasets,
};
use chrono_experiments::print_utils::{print_mass_global_summary, print_mass_statistics};
use chrono_experiments::proces_utils::{apply_mad_filter, interpolate_space_time};
use chrono_experiments::{AnalysisConfig, SpaceTimeCoordinate};
use deep_causality_num::{Float106, RealField};
use deep_causality_physics::{ChronoGauge, ChronoGaugeOps};
use deep_causality_topology::Lattice;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;
use std::io;
use std::sync::Arc;
use std::sync::Mutex;

/// Change this to `f64` for standard precision or `Float106` for quad-precision.
type FloatType = Float106;

/// Satellite ID for this experiment (Galileo E14).
pub const SAT_ID: &str = "E14";
pub const DBG: bool = false;

/// Macro to convert f64 literals to FloatType.
macro_rules! flt {
    ($val:expr) => {
        FloatType::from($val)
    };
}

// =============================================================================
// SETUP & INITIALIZATION
// =============================================================================

fn main() -> io::Result<()> {
    println!("╔══════════════════════════════════════════════════════════════════════╗");
    println!("║  E00: CHRONO-MASS EXPERIMENT                                         ║");
    println!("║  GM Derivation from Time Dilation via ChronoGauge                    ║");
    println!("╠══════════════════════════════════════════════════════════════════════╣");
    println!("║  Float Type: {:<56}║", std::any::type_name::<FloatType>());
    println!("║  Satellite:  E14 (Galileo)                                           ║");
    println!("║  Theory:     Chrono-Gauge Lattice (U(1) × SU(2)                      ║");
    println!("╚══════════════════════════════════════════════════════════════════════╝\n");

    // Verify data path
    let data_path_buf = get_gnss_data_input_path();
    let data_path = data_path_buf.to_str().unwrap();

    // Create the ChronoGauge field
    let gauge_field = create_chrono_gauge();
    // Run the experiment
    run_experiment(&gauge_field, data_path)?;

    Ok(())
}

fn create_chrono_gauge() -> ChronoGauge<FloatType> {
    // 4D lattice, size 32^4, periodic boundaries
    let shape = [32, 32, 32, 32];
    let periodic = [true, true, true, true];
    let lattice = Arc::new(Lattice::new(shape, periodic));

    // Initialize with identity configuration
    // beta = 1.0 (coupling constant parameter, irrelevant for static analysis)
    ChronoGauge::identity(lattice, flt!(1.0))
}

// =============================================================================
// EXPERIMENT EXECUTION
// =============================================================================

/// Runs the chrono-mass experiment on real GNSS data.
fn run_experiment(gauge_field: &ChronoGauge<FloatType>, data_path: &str) -> io::Result<()> {
    let mut yearly_results: Vec<(String, Vec<FloatType>)> = Vec::new();

    for year in YEARS {
        let results = run_year_analysis(year, gauge_field, data_path)?;
        yearly_results.push((year.to_string(), results));
    }

    let mut master_results: Vec<FloatType> = Vec::new();
    for (year, results) in &yearly_results {
        print_mass_global_summary(results, year);
        master_results.extend(results.iter().cloned());
    }

    print_mass_global_summary(&master_results, "TOTAL (All Years)");

    Ok(())
}

/// Analyze a single year of data.
fn run_year_analysis(
    year: &str,
    gauge_field: &ChronoGauge<FloatType>,
    data_path: &str,
) -> io::Result<Vec<FloatType>> {
    let datasets = get_year_datasets(year);
    let global_results = Mutex::new(Vec::new());

    // Process datasets in parallel
    datasets.par_iter().for_each(|dataset| {
        // Extract GPS dataset ID from filename (e.g., "gbm19670" -> 19670)
        if let Some(dataset_id) = extract_gps_dataset_id(dataset)
            && ANOMALOUS_WEEKS.contains(&dataset_id)
        {
            if DBG {
                println!(
                    "[{}] Skipping anomalous dataset {} (2017/2018 Data Crisis)",
                    dataset, dataset_id
                );
            }
            return;
        }

        let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
        let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);

        match process_dataset(dataset, &clk_path, &sp3_path, gauge_field) {
            Ok(results) => {
                let mut acc = global_results.lock().unwrap();
                acc.extend(results);
            }
            Err(e) => {
                eprintln!("Failed to process dataset {}: {}", dataset, e);
            }
        }
    });

    let final_results = global_results.into_inner().unwrap();
    Ok(final_results)
}

/// Extract GPS dataset identifier from dataset name (e.g., "gbm19670" -> 19670)
fn extract_gps_dataset_id(dataset: &str) -> Option<u32> {
    // The dataset string is typically just the basename without extension, e.g., "gbm19670"
    // "gbm" (3 chars) + Week (4 chars) + Day (1 char) = 5-digit identifier
    if dataset.len() >= 8 && dataset.starts_with("gbm") {
        let id_str = &dataset[3..8];
        id_str.parse::<u32>().ok()
    } else {
        None
    }
}

/// Processes a single dataset file pair.
fn process_dataset(
    dataset_name: &str,
    clk_path: &str,
    sp3_path: &str,
    gauge_field: &ChronoGauge<FloatType>,
) -> Result<Vec<FloatType>, io::Error> {
    let config = AnalysisConfig::default();

    // Load GNSS data
    let dm = DataManager::default();
    let (clocks, orbits) = dm.load_gnss_single_satellite(clk_path, sp3_path, SAT_ID)?;

    // Interpolate orbits to clock timestamps using 10th-order Lagrange polynomial
    let data: Vec<SpaceTimeCoordinate<FloatType>> = interpolate_space_time(&clocks, &orbits);

    if DBG {
        println!(
            "[{}] Interpolated {} orbit × {} clock → {} space-time coords",
            dataset_name,
            orbits.len(),
            clocks.len(),
            data.len()
        );
    }

    // Skip datasets with insufficient data
    if data.len() <= config.window_size_indices {
        if DBG {
            println!("  → Insufficient data points for window size");
        }
        return Ok(Vec::new());
    }

    // Derive GM values using ChronoGaugeWitness::source()
    let mut raw_gm_values: Vec<FloatType> = Vec::new();

    let mut i = 0;
    while i < data.len() - config.window_size_indices {
        let idx_a = i;
        let idx_b = i + config.window_size_indices;

        let r_a: FloatType = data[idx_a].r_m;
        let r_b: FloatType = data[idx_b].r_m;
        let d_h: FloatType = (r_a - r_b).abs();

        // Skip if height difference is too small
        if d_h < config.min_height_diff_m {
            i += 1;
            continue;
        }

        // Use source() to invert the Einstein field equation
        if let Ok(gm) = gauge_field.source(&data[idx_a], &data[idx_b]) {
            raw_gm_values.push(gm);
        }

        i += config.step_size;
    }

    if raw_gm_values.is_empty() {
        println!("  → No valid GM derivations");
        return Ok(Vec::new());
    }

    // Apply MAD filter for outlier rejection (now generic over Float)
    let filtered = apply_mad_filter(&raw_gm_values, flt!(config.outlier_sigma));

    if DBG {
        println!(
            "  → {} raw → {} filtered GM values",
            raw_gm_values.len(),
            filtered.len()
        );
        print_mass_statistics(&filtered)?;
    }

    Ok(filtered)
}
