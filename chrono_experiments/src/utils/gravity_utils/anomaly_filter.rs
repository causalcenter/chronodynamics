/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Anomaly date filtering for chrono_mass experiment.
//!
//! Loads anomaly windows detected by the Gravitational Observatory
//! to filter epochs with compromised clock data.
//! Also filters to standard sample months (3, 6, 9, 12) for fair comparison.

use chrono::NaiveDate;

/// Standard sample months for fair year-to-year comparison.
/// These are March, June, September, December.
pub const STANDARD_SAMPLE_MONTHS: [u32; 4] = [3, 6, 9, 12];

/// An anomaly window representing a date range with elevated clock instability.
#[derive(Debug, Clone)]
pub struct AnomalyWindow {
    pub start: NaiveDate,
    pub end: NaiveDate,
    pub satellite: String,
    pub tier: String,
}

/// Returns known anomaly windows per year from Observatory analysis.
/// FILTERED TO E14 ONLY - chrono_mass uses only E14, so we only filter
/// dates where E14 specifically was flagged as anomalous.
pub fn get_anomaly_windows(year: &str) -> Vec<AnomalyWindow> {
    match year {
        "2016" => vec![
            // NO E14 anomalies detected in 2016
            // (E11, E12, E19 had anomalies, but not E14)
        ],
        "2017" => vec![
            // ONLY E14-specific anomaly detected by Observatory
            AnomalyWindow {
                start: NaiveDate::from_ymd_opt(2017, 9, 24).unwrap(),
                end: NaiveDate::from_ymd_opt(2017, 9, 24).unwrap(),
                satellite: "E14".to_string(),
                tier: "ELEVATED".to_string(),
            },
            // Note: E19, E11, E09 had many anomalies but E14 was generally healthy
        ],
        "2018" => vec![
            // NO E14 anomalies detected in 2018
            // (E09, E12 had massive anomaly cluster, but not E14)
        ],
        _ => vec![],
    }
}

/// Check if a timestamp is in a standard sample month (March, June, September, December).
pub fn is_standard_sample_month(timestamp: u64) -> bool {
    use chrono::{DateTime, Datelike, Utc};
    if let Some(dt) = DateTime::<Utc>::from_timestamp(timestamp as i64, 0) {
        let month = dt.month();
        STANDARD_SAMPLE_MONTHS.contains(&month)
    } else {
        false
    }
}

/// Check if a timestamp should be filtered (either in anomaly window OR non-standard month).
pub fn should_filter_epoch(timestamp: u64, windows: &[AnomalyWindow], filter_months: bool) -> bool {
    // Filter non-standard months if requested
    if filter_months && !is_standard_sample_month(timestamp) {
        return true;
    }

    // Check anomaly windows
    is_in_anomaly_window(timestamp, windows)
}

/// Check if a timestamp falls within any anomaly window.
pub fn is_in_anomaly_window(timestamp: u64, windows: &[AnomalyWindow]) -> bool {
    use chrono::{DateTime, Utc};
    let datetime = DateTime::<Utc>::from_timestamp(timestamp as i64, 0);
    if let Some(dt) = datetime {
        let date = dt.date_naive();
        for window in windows {
            if date >= window.start && date <= window.end {
                return true;
            }
        }
    }
    false
}

/// Count how many epochs fall within anomaly windows.
pub fn count_filtered_epochs(timestamps: &[u64], windows: &[AnomalyWindow]) -> (usize, usize) {
    let total = timestamps.len();
    let filtered = timestamps
        .iter()
        .filter(|&&ts| !is_in_anomaly_window(ts, windows))
        .count();
    (total, filtered)
}
