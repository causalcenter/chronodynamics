/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! E02 Chrono-Gauge Experiment
//!
//!
//! | DEC Observable | LGT Method | Reference Value |
//! |----------------|------------|-----------------|
//! | Laplacian      | Wilson action (`mass_density_action()`) | 0.00 |
//! | Curl           | SU(2) plaquettes (`vorticity_tensor()`) | 0.00 |
//! | Gradient       | U(1) link phases | ~1e-10 |
//! | Vorticity-Z    | `vorticity_z()` | 0.00 |
//! | Time-Vel Corr  | Correlations | +1.00 |
//! | Tolman         | `tolman_temperature()` | 1.00 |

use chrono_data_manager::{
    ANOMALOUS_WEEKS, DataManager, YEARS, get_gnss_data_input_path, get_year_datasets,
};
use chrono_experiments::gauge_utils::{
    compute_dataset_time_velocity_correlation, compute_epoch_metrics_fast,
};
use chrono_experiments::{
    EpochMetrics, GaugeValidationResult, SpaceTimeCoordinate, folder_utils, gauge_utils,
    gravity_utils::derive_mass_from_dataset,
    print_utils::{print_cross_validation_table, print_gauge_header, print_gauge_summary},
    proces_utils::interpolate_space_time_single_pass,
    statistics_utils::calculate_pca,
};
use deep_causality_num::{Float106, RealField};
use deep_causality_physics::{ChronoGauge, ChronoGaugeOps, EARTH_J2};
use deep_causality_topology::{Lattice, LatticeGaugeField};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io;
use std::sync::Arc;
use std::sync::Mutex;
// =============================================================================
// CONFIGURATION
// =============================================================================

/// Float type for calculations. Use `Float106` for quad-precision.
type FloatType = Float106;

/// Minimum satellites needed per epoch for LGT construction.
const MIN_SATELLITES: usize = 6;

// =============================================================================
// MAIN
// =============================================================================

fn main() -> io::Result<()> {
    print_gauge_header();

    let data_path_buf = get_gnss_data_input_path();
    let data_path = data_path_buf.to_str().unwrap();

    println!("📊 Data path: {}", data_path);
    folder_utils::check_folder(data_path)?;
    println!();

    // Collect all validation results
    // Create the ChronoGauge field (Identity configuration for reference)
    let gauge_field = create_chrono_gauge();

    let mut all_results: Vec<GaugeValidationResult<FloatType>> = Vec::new();

    for year in YEARS {
        let year_results = process_year(year, data_path, &gauge_field)?;
        print_gauge_summary(&year_results, year)?;
        all_results.extend(year_results);
    }

    // PCA Analysis
    run_pca_analysis(&all_results);

    // J2 Analysis
    // run_j2_analysis(data_path, &gauge_field);

    // Pass None for global correlation override since we now compute it per-dataset correctly
    print_cross_validation_table(&all_results, None);

    println!("\n✅ E02 Chrono-Gauge experiment complete.");

    Ok(())
}

fn run_j2_analysis(data_path: &str, gauge_field: &ChronoGauge<FloatType>) {
    // J2 Oblateness Estimation
    println!("\n📐 Running J2 Oblateness Estimation (ChronoGaugeOps::solve_j2)...");
    let all_coords =
        collect_all_coordinates(data_path).expect("couldn't collect all coordinates data");
    println!(
        "   Collected {} coordinates for J2 analysis.",
        all_coords.len()
    );

    if !all_coords.is_empty() {
        // Use the Gauge Field Witness to solve for J2
        // This leverages the HKT implementation in deep_causality_physics
        match gauge_field.solve_j2(&all_coords) {
            Ok(derived_j2) => {
                // Reference J2 (JGM-3)
                let reference_j2 = FloatType::from(EARTH_J2);
                let error_percent =
                    ((derived_j2 - reference_j2) / reference_j2).abs() * FloatType::from(100.0);

                println!(
                    "╔══════════════════════════════════════════════════════════════════════════════╗"
                );
                println!(
                    "║  J2 ESTIMATION RESULT (via ChronoGaugeWitness HKT)                          ║"
                );
                println!(
                    "╠══════════════════════════════════════════════════════════════════════════════╣"
                );
                println!(
                    "║    Data Points:       {:>10}                                              ║",
                    all_coords.len()
                );
                println!(
                    "╠══════════════════════════════════════════════════════════════════════════════╣"
                );
                println!(
                    "║    Derived J2:        {:>12.6e}                                          ║",
                    derived_j2
                );
                println!(
                    "║    Reference J2:      {:>12.6e}  (JGM-3)                                  ║",
                    reference_j2
                );
                println!(
                    "║    Error:              {:.4e}%                                              ║",
                    error_percent
                );
                println!(
                    "╚══════════════════════════════════════════════════════════════════════════════╝"
                );
            }
            Err(e) => {
                println!("⚠️ J2 gauge calculation failed: {}", e);
            }
        }
    } else {
        println!("⚠️ No data collected for J2 estimation.");
    }
}

// =============================================================================
// YEAR PROCESSING
// =============================================================================

fn process_year(
    year: &str,
    data_path: &str,
    gauge_field: &ChronoGauge<FloatType>,
) -> io::Result<Vec<GaugeValidationResult<FloatType>>> {
    let datasets = get_year_datasets(year);

    let results: Mutex<Vec<GaugeValidationResult<FloatType>>> = Mutex::new(Vec::new());

    datasets.par_iter().for_each(|dataset| {
        let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
        let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);

        match process_dataset(dataset, &clk_path, &sp3_path, gauge_field) {
            Ok(Some(result)) => {
                if result.num_epochs > 0 {
                    results.lock().unwrap().push(result);
                }
            }
            Ok(None) => {
                // No data or insufficient data for this dataset, skip.
            }
            Err(e) => {
                eprintln!("[{}] Error: {}", dataset, e);
            }
        }
    });

    Ok(results.into_inner().unwrap())
}

// =============================================================================
// DATASET PROCESSING
// =============================================================================

fn process_dataset(
    dataset_name: &str,
    clk_path: &str,
    sp3_path: &str,
    gauge_field: &ChronoGauge<FloatType>,
) -> Result<Option<GaugeValidationResult<FloatType>>, io::Error> {
    // Load GNSS data
    let dm = DataManager::default();
    let (clocks, orbits) = dm.load_gnss_all_satellites(clk_path, sp3_path)?;

    if clocks.is_empty() || orbits.is_empty() {
        return Ok(Some(GaugeValidationResult::default()));
    }

    // Get satellite IDs and interpolate
    let sat_ids: Vec<u32> = clocks
        .iter()
        .map(|c| c.sat_id().as_num())
        .chain(orbits.iter().map(|o| o.sat_id().as_num()))
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();

    let mut all_coordinates: Vec<SpaceTimeCoordinate<FloatType>> = Vec::new();

    for sat_num in sat_ids {
        let sat_clocks: Vec<_> = clocks
            .iter()
            .filter(|c| c.sat_id().as_num() == sat_num)
            .cloned()
            .collect();

        let mut sat_orbits: Vec<_> = orbits
            .iter()
            .filter(|o| o.sat_id().as_num() == sat_num)
            .cloned()
            .collect();

        if sat_clocks.is_empty() || sat_orbits.len() < 10 {
            continue;
        }

        sat_orbits.sort_by_key(|a| a.timestamp());
        let coords = interpolate_space_time_single_pass(&sat_clocks, &sat_orbits);
        all_coordinates.extend(coords);
    }

    // Group by epoch
    let mut epoch_groups: HashMap<i64, Vec<SpaceTimeCoordinate<FloatType>>> = HashMap::new();
    for coord in &all_coordinates {
        let ts = coord.timestamp as i64;
        epoch_groups.entry(ts).or_default().push(*coord);
    }

    // Process each epoch - FAST: no lattice creation
    let mut epoch_metrics: Vec<EpochMetrics<FloatType>> = Vec::new();

    for (_ts, satellites) in epoch_groups.iter() {
        if satellites.len() < MIN_SATELLITES {
            continue;
        }

        // Compute observables directly from satellite data (no lattice!)
        let metrics = compute_epoch_metrics_fast(satellites);
        epoch_metrics.push(metrics);
    }

    let mut result = GaugeValidationResult::from_epochs(
        dataset_name.to_string(),
        &epoch_metrics,
        epoch_metrics.len(),
    );

    // FIX: Recompute Momentum (Correlation) using ONLY coordinates from valid epochs (>= MIN_SATELLITES)
    // This matches legacy logic which only processed valid epochs.
    let mut valid_coordinates: Vec<SpaceTimeCoordinate<FloatType>> = Vec::new();
    for (_ts, satellites) in epoch_groups.iter() {
        if satellites.len() >= MIN_SATELLITES {
            valid_coordinates.extend(satellites.iter().cloned());
        }
    }

    let dataset_corr = compute_dataset_time_velocity_correlation::<FloatType>(&valid_coordinates);

    // DEBUG: Print datasets with significant E14/E18 data but poor correlation
    // This helps identify if specific files are corrupt or "weird"
    let corr_f64: f64 = dataset_corr.into();
    if corr_f64 < 0.9 {
        let e14_18_count = valid_coordinates
            .iter()
            .filter(|s| s.sat_id == 14 || s.sat_id == 18)
            .count();

        if e14_18_count > 50 {
            println!(
                "⚠️  [{}]: Data Points: {}, Correlation: {:.4} (LOW)",
                dataset_name, e14_18_count, corr_f64
            );
        }
    }

    result.time_velocity_correlation = dataset_corr;

    // Derived Rotation using ONLY coordinates from valid epochs
    let derived_rot = gauge_utils::estimate_earth_rotation::<FloatType>(&valid_coordinates);
    result.derived_earth_rotation = derived_rot;

    // We construct paths assuming standard layout (which we have in arguments)
    // Only derive for E14 to match E01 methodology, or if specific satellite requested.
    // The E02 experiment iterates ALL datasets. We should try to derive mass for this dataset.
    // derive_mass_from_dataset handles file loading internally to ensure fresh two-pass interpolation.
    if let Ok(Some(mass_res)) = derive_mass_from_dataset(
        dataset_name,
        clk_path,
        sp3_path,
        "E14", // Try E14 first
        gauge_field,
        ANOMALOUS_WEEKS,
    ) {
        result.derived_gm = mass_res.derived_gm;
        result.derived_mass = mass_res.derived_mass;
        result.derived_gravity = mass_res.derived_gravity;
    } else {
        // If E14 failed or not present, try "E18" or fallback to primary if possible?
        // For now, we stick to E14 as the "Gold Standard" mass witness.
        // If satellite is not in file, it returns None.
    }

    Ok(Some(result))
}

fn create_chrono_gauge() -> ChronoGauge<FloatType> {
    // 4D lattice, size 32^4, periodic boundaries (same as E01)
    let shape = [32, 32, 32, 32];
    let periodic = [true, true, true, true];
    let lattice = Arc::new(Lattice::new(shape, periodic));

    // Initialize with identity configuration
    // beta = 1.0 (coupling constant parameter, irrelevant for static analysis)
    let one = <FloatType as From<f64>>::from(1.0);
    let lattice_gauge: ChronoGauge<FloatType> = LatticeGaugeField::identity(lattice, one);

    lattice_gauge
}

// =============================================================================
// COORDINATE COLLECTION FOR J2 ESTIMATION
// =============================================================================

/// Collects all SpaceTimeCoordinates from all datasets for J2 analysis.
/// Uses Rayon for parallel processing across datasets.
fn collect_all_coordinates(data_path: &str) -> io::Result<Vec<SpaceTimeCoordinate<FloatType>>> {
    // Collect all dataset paths first
    let mut all_dataset_paths: Vec<(String, String, String)> = Vec::new();

    for year in YEARS {
        let datasets = get_year_datasets(year);
        for dataset in datasets {
            let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
            let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);
            all_dataset_paths.push((year.to_string(), clk_path, sp3_path));
        }
    }

    // Process datasets in parallel
    let all_coords: Mutex<Vec<SpaceTimeCoordinate<FloatType>>> = Mutex::new(Vec::new());

    all_dataset_paths
        .par_iter()
        .for_each(|(_year, clk_path, sp3_path)| {
            if !std::path::Path::new(clk_path).exists() {
                return;
            }

            let dm = DataManager::default();
            let Ok((clocks, orbits)) = dm.load_gnss_all_satellites(clk_path, sp3_path) else {
                return;
            };

            if clocks.is_empty() || orbits.is_empty() {
                return;
            }

            // Get satellite IDs
            let sat_ids: Vec<u32> = clocks
                .iter()
                .map(|c| c.sat_id().as_num())
                .chain(orbits.iter().map(|o| o.sat_id().as_num()))
                .collect::<std::collections::HashSet<_>>()
                .into_iter()
                .collect();

            let mut local_coords: Vec<SpaceTimeCoordinate<FloatType>> = Vec::new();

            for sat_num in sat_ids {
                // FILTER: Exclude E14 (14) and E18 (18) - eccentric satellites
                // These have high eccentricity (~0.16) which creates massive altitude
                // variations that contaminate the J2 regression with monopole effects.
                if sat_num == 14 || sat_num == 18 {
                    continue;
                }
                let sat_clocks: Vec<_> = clocks
                    .iter()
                    .filter(|c| c.sat_id().as_num() == sat_num)
                    .cloned()
                    .collect();

                let mut sat_orbits: Vec<_> = orbits
                    .iter()
                    .filter(|o| o.sat_id().as_num() == sat_num)
                    .cloned()
                    .collect();

                if sat_clocks.is_empty() || sat_orbits.len() < 10 {
                    continue;
                }

                sat_orbits.sort_by_key(|a| a.timestamp());
                let coords = interpolate_space_time_single_pass(&sat_clocks, &sat_orbits);
                local_coords.extend(coords);
            }

            // Merge into global collection
            all_coords.lock().unwrap().extend(local_coords);
        });

    Ok(all_coords.into_inner().unwrap())
}

// =============================================================================
// PCA ANALYSIS
// =============================================================================

fn run_pca_analysis(results: &[GaugeValidationResult<FloatType>]) {
    if results.len() < 3 {
        println!("\n⚠ Insufficient data for PCA analysis (need at least 3 observations)\n");
        return;
    }

    // Extract 7 observables for PCA
    // Construct PCA matrix (observations x variables)
    // Format: Each row is a sample (dataset), each column is a feature (metric)
    let pca_matrix: Vec<Vec<FloatType>> = results
        .iter()
        .map(|r| {
            vec![
                r.mean_mass_density,
                r.mean_curl_magnitude,
                r.mean_gradient_magnitude,
                r.mean_vorticity_z,
                r.time_velocity_correlation,
                r.virial_ratio,
                r.tolman_consistency,
                r.derived_earth_rotation,
            ]
        })
        .collect();

    let variable_names = [
        "Mass/Tidal (Laplacian)",
        "Frame-Dragging (Curl)",
        "Geoid (Gradient)",
        "Spin (Vorticity-Z)",
        "Momentum (Time-Vel Corr)",
        "Energy (Virial Ratio)",
        "Thermodynamics (Tolman)",
        "Earth Rotation (Derived)",
    ];

    // Perform PCA
    let pca = calculate_pca(&pca_matrix, 3);

    println!("\n╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║               PRINCIPAL COMPONENT ANALYSIS                                   ║");
    println!("║               Dimensionality Reduction of Validation Metrics                 ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝\n");

    let num_obs = results.len();
    let num_vars = 8;

    println!(
        "📊 PCA Results ({} observations, {} variables)",
        num_obs, num_vars
    );
    println!("───────────────────────────────────────────────────────────────────────────────\n");

    // Compute cumulative variance from variance_explained
    let mut cumulative = 0.0f64;
    let cumulative_variances: Vec<f64> = pca
        .variance_explained
        .iter()
        .map(|&v| {
            let v_f64: f64 = v.into();
            cumulative += v_f64;
            cumulative
        })
        .collect();

    // Variance explained
    println!("Variance Explained by Principal Components:");
    for (i, var_exp) in pca.variance_explained.iter().enumerate() {
        let var_pct: f64 = (*var_exp).into();
        let cum_pct = cumulative_variances.get(i).copied().unwrap_or(var_pct);
        println!(
            "  PC{}: {:.2}% (Cumulative: {:.2}%)",
            i + 1,
            var_pct * 100.0,
            cum_pct * 100.0
        );
    }

    println!("\n───────────────────────────────────────────────────────────────────────────────");
    println!("Principal Component Loadings (Contribution of each variable):\n");

    // eigenvectors contains the loadings
    for (pc_idx, loadings) in pca.eigenvectors.iter().enumerate() {
        let var_pct: f64 = pca
            .variance_explained
            .get(pc_idx)
            .copied()
            .map(|v| v.into())
            .unwrap_or(0.0);
        println!(
            "PC{} (explains {:.1}% of variance):",
            pc_idx + 1,
            var_pct * 100.0
        );
        for var_idx in 0..loadings.len().min(variable_names.len()) {
            let loading_f64: f64 = loadings[var_idx].into();
            let abs_loading = loading_f64.abs();
            let bar_length = (abs_loading * 40.0) as usize;
            let bar = "█".repeat(bar_length);
            let sign = if loading_f64 >= 0.0 { "+" } else { "-" };
            println!(
                "  {:30} {} {:.3} {}",
                variable_names[var_idx], sign, abs_loading, bar
            );
        }
        println!();
    }

    println!("───────────────────────────────────────────────────────────────────────────────");
    println!("Interpretation Guide:");
    println!("  • PC1 typically represents the dominant validation pattern");
    println!("  • High loadings (>0.5) indicate strong contribution to that component");
    println!("  • Components with >70% cumulative variance capture most information");
    println!("  • Variables with similar loadings are correlated");
    println!("═══════════════════════════════════════════════════════════════════════════════\n");
}
