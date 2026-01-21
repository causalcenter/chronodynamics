/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Save 5T Observatory results to CSV files.

use crate::ObservatoryResults;
use std::fs::File;
use std::io::{self, Write};

/// Time series data point for 5T measurements
#[derive(Debug, Clone)]
pub struct TimeSeriesPoint {
    pub timestamp: i64,
    pub e12_clock_ns: f64,
    pub e12_radius_m: f64,
    pub e14_clock_ns: f64,
    pub e14_radius_m: f64,
    pub e18_clock_ns: f64,
    pub e18_radius_m: f64,
    pub clock_diff_ns: f64,
}

/// Save observatory results and time series to CSV files.
///
/// # Arguments
/// * `results` - The computed observatory results
/// * `time_series` - Vector of time series data points
/// * `year_out_path` - Output directory path for the year
/// * `year` - Year label for the output files
pub fn save_5t_results(
    results: &ObservatoryResults,
    time_series: &[TimeSeriesPoint],
    year_out_path: &str,
    year: &str,
) -> io::Result<()> {
    // Save summary
    let summary_path = format!("{}/{}_observatory_summary.csv", year_out_path, year);
    let mut summary_file = File::create(&summary_path)?;

    writeln!(summary_file, "Metric,Value")?;
    writeln!(summary_file, "Year,{}", year)?;
    writeln!(summary_file, "Total_Epochs,{}", results.total_epochs)?;
    writeln!(summary_file, "Valid_5T_Epochs,{}", results.valid_5t_epochs)?;
    writeln!(
        summary_file,
        "Mean_Clock_Diff_ns,{:.10e}",
        results.mean_clock_diff_e14_e12
    )?;
    writeln!(
        summary_file,
        "Std_Clock_Diff_ns,{:.10e}",
        results.std_clock_diff_e14_e12
    )?;
    writeln!(
        summary_file,
        "Mean_Altitude_Diff_m,{:.6}",
        results.mean_altitude_diff_m
    )?;
    writeln!(
        summary_file,
        "Predicted_GR_Shift_ns_s,{:.10e}",
        results.predicted_gr_shift_ns
    )?;
    writeln!(
        summary_file,
        "Observed_Shift_ns_s,{:.10e}",
        results.observed_shift_ns
    )?;
    writeln!(
        summary_file,
        "GR_Agreement_Pct,{:.6}",
        results.gr_agreement_pct
    )?;
    writeln!(
        summary_file,
        "Clock_Rate_Altitude_Correlation,{:.6}",
        results.clock_rate_altitude_correlation
    )?;
    writeln!(
        summary_file,
        "Measured_GR_Coefficient,{:.10e}",
        results.measured_gr_coefficient
    )?;
    writeln!(
        summary_file,
        "Predicted_GR_Coefficient,{:.10e}",
        results.predicted_gr_coefficient
    )?;
    writeln!(
        summary_file,
        "GR_Coefficient_Agreement_Pct,{:.6}",
        results.gr_coefficient_agreement_pct
    )?;
    writeln!(
        summary_file,
        "Detrended_Clock_Std_ns,{:.10e}",
        results.detrended_clock_std_ns
    )?;
    writeln!(
        summary_file,
        "E12_E18_Correlation,{:.6}",
        results.e12_e18_correlation
    )?;
    writeln!(
        summary_file,
        "E14_Temporal_Autocorr,{:.6}",
        results.temporal_autocorr_e14
    )?;
    // GW Detection metrics
    writeln!(
        summary_file,
        "GW_Strain_Upper_Limit,{:.10e}",
        results.gw_strain_upper_limit
    )?;
    writeln!(
        summary_file,
        "Quadrupole_Residual_ns,{:.10e}",
        results.quadrupole_residual
    )?;
    writeln!(
        summary_file,
        "Cross_Correlation_Coeff,{:.6}",
        results.cross_correlation_coeff
    )?;
    writeln!(
        summary_file,
        "Spectral_Power_nHz,{:.10e}",
        results.spectral_power_n_hz
    )?;
    // Hellings-Downs analysis
    writeln!(
        summary_file,
        "Geodesic_Strain,{:.10e}",
        results.geodesic_strain
    )?;
    writeln!(
        summary_file,
        "Angular_Separation_Deg,{:.2}",
        results.angular_separation_deg
    )?;
    writeln!(
        summary_file,
        "Expected_HD_Correlation,{:.6}",
        results.expected_hd_correlation
    )?;
    writeln!(
        summary_file,
        "HD_Significance_Pct,{:.2}",
        results.hd_significance_pct
    )?;

    println!("  → Saved: {}", summary_path);

    // Save time series (first 10000 epochs for analysis)
    let timeseries_path = format!("{}/{}_5t_timeseries.csv", year_out_path, year);
    let mut ts_file = File::create(&timeseries_path)?;

    writeln!(
        ts_file,
        "Timestamp,E12_Clock_ns,E12_Radius_m,E14_Clock_ns,E14_Radius_m,E18_Clock_ns,E18_Radius_m,Clock_Diff_ns"
    )?;

    for point in time_series.iter().take(10000) {
        writeln!(
            ts_file,
            "{},{:.6},{:.3},{:.6},{:.3},{:.6},{:.3},{:.6}",
            point.timestamp,
            point.e12_clock_ns,
            point.e12_radius_m,
            point.e14_clock_ns,
            point.e14_radius_m,
            point.e18_clock_ns,
            point.e18_radius_m,
            point.clock_diff_ns
        )?;
    }

    println!("  → Saved: {}", timeseries_path);

    Ok(())
}
