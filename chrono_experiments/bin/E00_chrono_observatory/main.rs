/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */
// Physics constants from deep_causality_physics crate
use chrono_experiments::utils::observatory_utils::chrono_gravity_utils::{
    TIME_DELAY_S, analyze_5t_fields, analyze_swarm_vectors, build_5t_fields,
};
use chrono_experiments::utils::observatory_utils::data_loader::merge_satellite_data;
use chrono_experiments::utils::observatory_utils::linear_algebra::calculate_pca;
use chrono_experiments::utils::observatory_utils::residual_filter::apply_single_sat_cleaning;
use chrono_experiments::utils::observatory_utils::windowing::{
    AnomalyTier, scan_temporal_anomalies,
};
use chrono_experiments::{IncidentRecord, OBSERVATORY_CONFIG, ObservatoryConfigParams};
use deep_causality_num::{Float106, RealField};
use deep_causality_physics::NEWTONIAN_CONSTANT_OF_GRAVITATION as G_REF;
use deep_causality_physics::{EARTH_GM, SPEED_OF_LIGHT};
// Added
use chrono_data_manager::{
    DataManager, ObservatoryEpoch, YEARS, get_data_output_path, get_gnss_data_input_path,
    get_year_datasets,
};
use chrono_experiments::print_utils::print_observatory_config;

use chrono::{TimeZone, Utc};
use rayon::prelude::*;
use std::env; // Added for cmd args
use std::fs;
use std::io;
use std::sync::Mutex;

/// Change this to `f64` for standard precision or `Float106` for quad-precision.
type FloatType = Float106;

// CONSTANTS
const WINDOW_SIZE_EPOCHS: usize = 2016; // 1 week
const WINDOW_STEP_EPOCHS: usize = 288; // 1 day
const NOMINAL_ORBIT_RADIUS_M: f64 = deep_causality_physics::GALILEO_NOMINAL_RADIUS_M;

// Config Struct

fn main() -> io::Result<()> {
    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!("║        Forensic Clock Analysis & Fleet Health Monitoring                 ║");
    println!(
        "║       Float Type: {:<56} ║",
        std::any::type_name::<FloatType>()
    );
    println!("╚══════════════════════════════════════════════════════════════════════════╝");
    println!();

    // Print configuration using parameter-based function
    let config_params = ObservatoryConfigParams {
        g_ref: G_REF,
        c: SPEED_OF_LIGHT,
        gm_earth: EARTH_GM,
        time_delay_s: TIME_DELAY_S,
    };
    print_observatory_config(&config_params);

    // Ensure output directory exists
    let out_path_buf = get_data_output_path("gqcf_space_gravitational_observatory");
    let out_path = out_path_buf.to_str().unwrap();

    // Get GNSS data path from data_manager
    let data_path_buf = get_gnss_data_input_path();
    let data_path = data_path_buf.to_str().unwrap();

    fs::create_dir_all(out_path)?;

    let mut all_incidents: Vec<IncidentRecord> = Vec::new();

    for year in YEARS.iter() {
        let year_incidents = analyze_year::<FloatType>(year, out_path, data_path)?;
        all_incidents.extend(year_incidents);
    }

    // Check for --incidents flag
    let args: Vec<String> = env::args().collect();
    let append_mode = args.contains(&"--incidents".to_string());

    if append_mode {
        // Append to incident_catalogue.md in the output dir.
        let report_path = format!("{}/incident_catalogue.md", out_path);
        generate_incident_catalogue(&all_incidents, Some(&report_path));
    } else {
        // Just print to stdout
        generate_incident_catalogue(&all_incidents, None);
    }

    Ok(())
}

fn analyze_year<R>(year: &str, out_path: &str, data_path: &str) -> io::Result<Vec<IncidentRecord>>
where
    R: RealField
        + Copy
        + From<f64>
        + Into<f64>
        + PartialOrd
        + std::fmt::Display
        + std::fmt::LowerExp,
{
    println!("╔══════════════════════════════════════════════════════════════════════════╗");
    println!(
        "║  ANALYZING YEAR: {}                                                    ║",
        year
    );
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    let year_out_path = format!("{}/{}", out_path, year);
    fs::create_dir_all(&year_out_path)?;

    // Use datasets defined in data_manager
    let datasets = get_year_datasets(year);

    if datasets.is_empty() {
        println!("  ⚠ No datasets configured for year {}", year);
        return Ok(Vec::new());
    }

    let sat_ids: Vec<&str> = OBSERVATORY_CONFIG.iter().map(|c| c.id).collect();
    println!("\n  Loading data for satellites: {:?} (parallel)", sat_ids);
    println!("  Processing {} datasets...", datasets.len());

    // Thread-safe epoch accumulator
    let all_epochs: Mutex<Vec<ObservatoryEpoch>> = Mutex::new(Vec::new());

    // Parallel processing of datasets
    datasets.into_par_iter().for_each(|dataset| {
        let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
        let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);

        if !std::path::Path::new(&clk_path).exists() {
            println!("    Missing data file: {}", clk_path);
            return;
        }

        if !std::path::Path::new(&sp3_path).exists() {
            println!("    Missing data file: {}", sp3_path);
            return;
        }

        let mut local_epochs: Vec<ObservatoryEpoch> = Vec::new();

        for sat_config in OBSERVATORY_CONFIG.iter() {
            let dm = DataManager::default();
            match dm.load_gnss_single_satellite(&clk_path, &sp3_path, sat_config.id) {
                Ok((clocks, orbits)) => {
                    merge_satellite_data::<R>(&mut local_epochs, sat_config.id, &clocks, &orbits);
                }
                Err(e) => {
                    eprintln!(
                        "    ⚠ Error loading {} from {}: {}",
                        sat_config.id, dataset, e
                    );
                }
            }
        }

        {
            let mut global = all_epochs.lock().unwrap();
            for local_epoch in local_epochs {
                if let Some(existing) = global
                    .iter_mut()
                    .find(|e| e.timestamp == local_epoch.timestamp)
                {
                    existing.satellites.extend(local_epoch.satellites);
                } else {
                    global.push(local_epoch);
                }
            }
        }
    });

    let all_epochs = all_epochs.into_inner().unwrap();
    println!("\n  Total epochs collected: {}", all_epochs.len());

    if all_epochs.is_empty() {
        println!("  ⚠ No data loaded. Check files and satellite IDs.");
        return Ok(Vec::new());
    }

    // 1. SWARM ANALYSIS
    let swarm_results = analyze_swarm_vectors::<R>(&all_epochs);
    println!(
        "  Swarm Analysis: Quad Power = {:.4e}, Mono Power = {:.4e}",
        swarm_results.quadrupole_power, swarm_results.monopole_power
    );

    // 2. 5T FIELD ANALYSIS
    let five_t_fields = build_5t_fields::<R>(&all_epochs);
    println!("  Valid 5T measurements: {}", five_t_fields.len());

    if five_t_fields.is_empty() {
        println!("  ⚠ Insufficient data for 5T analysis");
    } else {
        let mut results = analyze_5t_fields::<R>(&five_t_fields);
        results.swarm_quadrupole_power = swarm_results.quadrupole_power;
        results.swarm_monopole_power = swarm_results.monopole_power;
        results.print(year);
    }

    // 4. PCA ANOMALY DETECTION
    let incidents = perform_pca_anomaly_analysis::<R>(&all_epochs, year);

    Ok(incidents)
}

fn perform_pca_anomaly_analysis<R>(epochs: &[ObservatoryEpoch], year: &str) -> Vec<IncidentRecord>
where
    R: RealField + Copy + From<f64> + Into<f64> + std::fmt::Display + PartialOrd,
{
    let mut collected_incidents = Vec::new();

    // Dynamic satellite list from config
    let sats: Vec<&str> = OBSERVATORY_CONFIG.iter().map(|c| c.id).collect();
    let sat_names: Vec<String> = sats.iter().map(|&s| s.to_string()).collect();

    let mut aligned_timestamps: Vec<i64> = Vec::new();
    let mut sat_data: Vec<Vec<R>> = vec![Vec::new(); sats.len()];
    let mut sat_radii: Vec<Vec<R>> = vec![Vec::new(); sats.len()];

    // Align data
    let mut sorted: Vec<&ObservatoryEpoch> = epochs.iter().collect();
    sorted.sort_by_key(|e| e.timestamp);

    for epoch in sorted {
        let mut present = true;
        for id in &sats {
            if !epoch.satellites.contains_key(*id) {
                present = false;
                break;
            }
        }

        if present {
            aligned_timestamps.push(epoch.timestamp);
            for (i, id) in sats.iter().enumerate() {
                let s = &epoch.satellites[*id];
                sat_data[i].push(R::from(s.clock_bias_ns));
                sat_radii[i].push(R::from(s.radius_m));
            }
        }
    }

    // Explicit usize types for length comparison
    let len = aligned_timestamps.len();
    if len < 10 {
        return collected_incidents;
    }

    // Clean residuals
    let r_nominal = R::from(NOMINAL_ORBIT_RADIUS_M);
    let mut cleaned_vectors: Vec<Vec<R>> = Vec::new();

    for i in 0..sats.len() {
        let cleaned =
            apply_single_sat_cleaning(&aligned_timestamps, &sat_data[i], &sat_radii[i], r_nominal);
        cleaned_vectors.push(cleaned.gw_candidates);
    }

    println!("\n╔══════════════════════════════════════════════════════════════════════════╗");
    println!(
        "║  PCA ANOMALY DETECTION: {}                                                ║",
        year
    );
    println!(
        "  Computing Principal Component Analysis on {} epochs...",
        aligned_timestamps.len()
    );

    // PCA
    let pca = calculate_pca(&cleaned_vectors, 5);

    for (i, &var_exp) in pca.variance_explained.iter().enumerate() {
        let var_pct: f64 = R::into(var_exp) * 100.0;
        println!(
            "║  PC{}: Variance Explained = {:>6.2}%                                      ║",
            i + 1,
            var_pct
        );

        // Print Dominant Sat for this PC (from Loadings)
        let mut max_loading = 0.0;
        let mut dom_idx = 0;
        if i < pca.loadings.len() {
            for (j, val) in pca.loadings[i].iter().enumerate() {
                let v: f64 = (*val).into();
                if v.abs() > max_loading.abs() {
                    max_loading = v;
                    dom_idx = j;
                }
            }
            println!(
                "║      Dominant Sat:     {:<3} (Loading: {:>+6.3})                            ║",
                sat_names[dom_idx], max_loading
            );
        }
    }

    // Interpretation Check
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  INTERPRETATION:                                                         ║");

    // Robust heuristic: Check for coherence in PC1 loadings
    // Common Mode = All satellites have significant loadings with the SAME SIGN.
    let is_coherent = if !pca.loadings.is_empty() {
        let loadings = &pca.loadings[0];
        let threshold = 0.25; // 1/sqrt(7) is approx 0.37, so 0.25 is a reasonable inclusion floor
        let first_sign = R::into(loadings[0]).signum();

        loadings.iter().all(|&val| {
            let v = R::into(val);
            v.abs() > threshold && v.signum() == first_sign
        })
    } else {
        false
    };

    if is_coherent {
        println!("║  INTERPRETATION: Strong Common Mode (Monopole) detected in PC1.          ║");
        println!("║  Effect: All clocks drifting coherently (Systematic hardware noise).     ║");
    } else {
        println!("║  WARNING: Eccentric satellites (E14/E18) may dominate PCs.            ║");
        println!("║  Orbital eccentricity is introducing significant systematic variance.    ║");
    }
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    // TEMPORAL SCAN
    let scan_results = scan_temporal_anomalies(
        &cleaned_vectors,
        &sat_names,
        WINDOW_SIZE_EPOCHS,
        WINDOW_STEP_EPOCHS,
    );

    println!("\n╔══════════════════════════════════════════════════════════════════════════╗");
    println!(
        "║  TEMPORAL ANOMALY SCAN ({}) - Empirical Null Distribution             ║",
        year
    );
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!(
        "║  Baseline Windows (first 40%):   {:<5}                                 ║",
        scan_results.baseline_windows
    );
    println!(
        "║  PC2 Null Mean:               {:<.2}%                                    ║",
        R::into(scan_results.null_dist.mean) * 100.0
    );
    println!(
        "║  PC2 Null Std Dev:            {:<.2}%                                    ║",
        R::into(scan_results.null_dist.std_dev) * 100.0
    );
    println!("╠══════════════════════════════════════════════════════════════════════════╣");

    // Thresholds
    let p90 = R::into(scan_results.null_dist.percentile(0.90));
    let p99 = R::into(scan_results.null_dist.percentile(0.99));
    let p999 = R::into(scan_results.null_dist.percentile(0.999));

    // Dynamic Bonferroni correction
    let n_tests = scan_results.test_windows as f64;
    let alpha_base = 0.1 / n_tests;

    println!("║  THRESHOLDS (From Empirical Null Distribution, Bonferroni-Corrected)    ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!(
        "║  90th percentile:   {:<6.2}%  (α_eff = {:.2e})                       ║",
        p90 * 100.0,
        alpha_base
    );
    println!(
        "║  99th percentile:   {:<6.2}%  (α_eff = {:.2e})                       ║",
        p99 * 100.0,
        alpha_base / 10.0
    );
    println!(
        "║  99.9th percentile: {:<6.2}%  (α_eff = {:.2e})                       ║",
        p999 * 100.0,
        alpha_base / 100.0
    );
    println!(
        "║  Total Windows Tested:   {:<5}                                          ║",
        scan_results.test_windows
    );
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  DETECTED ANOMALIES BY TIER                                             ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");

    let count_extreme = scan_results
        .anomalies
        .iter()
        .filter(|a| a.tier == AnomalyTier::Extreme)
        .count();
    let count_significant = scan_results
        .anomalies
        .iter()
        .filter(|a| a.tier == AnomalyTier::Significant)
        .count();
    let count_elevated = scan_results
        .anomalies
        .iter()
        .filter(|a| a.tier == AnomalyTier::Elevated)
        .count();

    println!(
        "║  > 99.9% (extreme):       {:<5} events                                     ║",
        count_extreme
    );
    println!(
        "║  > 99% (significant):     {:<5} events                                     ║",
        count_significant
    );
    println!(
        "║  > 90% (elevated):        {:<5} events                                     ║",
        count_elevated
    );

    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  DETECTED ANOMALIES                                                      ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  Rank | Date/Time (UTC)      | Sat  | PC2    | Pctl   | Tier            ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");

    // Sort anomalies by PC2 var (descending)
    let mut anomalies = scan_results.anomalies.clone();
    anomalies.sort_by(|a, b| {
        b.pc2_variance
            .partial_cmp(&a.pc2_variance)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    for (rank, anomaly) in anomalies.iter().enumerate().take(15) {
        let tier_str = match anomaly.tier {
            AnomalyTier::Extreme => "EXTREME ***",
            AnomalyTier::Significant => "SIGNIF **",
            AnomalyTier::Elevated => "ELEVATED *",
            _ => "NORMAL",
        };

        // Convert index to timestamp
        let ts = aligned_timestamps[anomaly.start_index];
        // Create DateTime from timestamp (seconds)
        let dt = Utc
            .timestamp_opt(ts, 0)
            .single()
            .unwrap_or_else(|| Utc.timestamp_opt(0, 0).unwrap());
        let date_str = dt.format("%Y-%m-%d %H:%M:%S").to_string();

        let pctl = (1.0 - anomaly.p_value) * 100.0;

        println!(
            "║   {:<3} | {:<19} | {:<4} | {:<5.1}% | {:<5.1}% | {:<15} ║",
            rank + 1,
            date_str,
            anomaly.dominant_sat,
            R::into(anomaly.pc2_variance) * 100.0,
            pctl,
            tier_str
        );

        // Collect into global list
        if anomaly.tier == AnomalyTier::Significant || anomaly.tier == AnomalyTier::Extreme {
            collected_incidents.push(IncidentRecord {
                timestamp: ts,
                date_str,
                satellite: anomaly.dominant_sat.clone(),
                pc2_variance: R::into(anomaly.pc2_variance) * 100.0,
                tier: anomaly.tier,
            });
        }
    }
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  METHODOLOGY NOTE                                                        ║");
    println!("╠══════════════════════════════════════════════════════════════════════════╣");
    println!("║  This analysis is forensic (retrospective). Thresholds derived from     ║");
    println!("║  the first 40% of data as empirical null distribution.                  ║");
    println!("╚══════════════════════════════════════════════════════════════════════════╝");

    collected_incidents
}

fn generate_incident_catalogue(incidents: &[IncidentRecord], file_path: Option<&str>) {
    let mut output = String::new();

    output.push_str("\n\n");
    output
        .push_str("╔══════════════════════════════════════════════════════════════════════════╗\n");
    output
        .push_str("║  GNSS INCIDENT CATALOGUE (2016-2018)                                     ║\n");
    output
        .push_str("╠══════════════════════════════════════════════════════════════════════════╣\n");
    output
        .push_str("║  Full forensic timeline of detected anomalies.                           ║\n");
    output
        .push_str("╚══════════════════════════════════════════════════════════════════════════╝\n");
    output.push('\n');

    // Sort by Date
    let mut sorted: Vec<&IncidentRecord> = incidents.iter().collect();
    sorted.sort_by_key(|i| i.timestamp);

    // Header
    output.push_str("| Rank | Date/Time (UTC)      | Satellite | PC2 Var | Status       | Phase             |\n");
    output.push_str("|:-----|:---------------------|:----------|:--------|:-------------|:------------------|\n");

    let mut rank = 1;

    for incident in sorted {
        // Only show Elevated and above
        if incident.tier == AnomalyTier::Normal {
            continue;
        }

        let status = match incident.tier {
            AnomalyTier::Extreme => "CRITICAL",
            AnomalyTier::Significant => "SEVERE",
            AnomalyTier::Elevated => "MODERATE",
            _ => "NORMAL",
        };

        // Date-based Phase tagging (simple heuristic matching user's timeline)
        let phase = if incident.date_str.starts_with("2016-06") {
            "Initial Signs"
        } else if incident.date_str.starts_with("2016-09") {
            "Early Warning"
        } else if incident.date_str.starts_with("2016-11") {
            if incident.satellite == "E11" && incident.pc2_variance > 25.0 {
                "Crisis Origin"
            } else {
                "Pre-Crisis"
            }
        } else if incident.date_str.starts_with("2017-01") {
            "Public Ack"
        } else if incident.date_str.starts_with("2017-09") {
            if incident.pc2_variance > 22.0 {
                "Peak Crisis"
            } else {
                "Escalation"
            }
        } else {
            "Persistent"
        };

        use std::fmt::Write;
        writeln!(
            &mut output,
            "| **{}** | **{}** | **{}** | **{:.1}%** | **{}** | {} |",
            rank, incident.date_str, incident.satellite, incident.pc2_variance, status, phase
        )
        .unwrap();

        rank += 1;
    }

    // Print to stdout
    println!("{}", output);

    // If file_path is provided, write/append to file
    if let Some(path) = file_path {
        use std::io::Write;
        let mut file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .expect("Failed to open report file");

        writeln!(file, "{}", output).expect("Failed to write to report file");
        println!("  [Report appended to: {}]", path);
    }
}
