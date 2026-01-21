//! Data loaders for MHD-related datasets.

use crate::OmniRecord;
use chrono::NaiveDate;
use deep_causality_num::RealField;
use std::fs::File;
use std::io::{self, BufRead};

/// Loads OMNI data from the specified file, filtering for September 2017.
pub fn load_omni_data<R>(path: &str) -> io::Result<Vec<OmniRecord<R>>>
where
    R: RealField + From<f64>,
{
    let file = File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut records = Vec::new();

    println!("  Loading OMNI data from: {}", path);

    for line in reader.lines() {
        let line = line?;

        // Skip header lines (check if first char is numeric)
        let trimmed = line.trim();
        if trimmed.is_empty() || !trimmed.chars().next().unwrap().is_numeric() {
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();

        // Ensure we have enough columns (OMNI format depends on variables selected)
        // We expect: Year DOY Hour ...
        if parts.len() < 3 {
            continue;
        }

        // 1. Parse Time
        let year_res = parts[0].parse::<i32>();
        let doy_res = parts[1].parse::<u32>();
        let hour_res = parts[2].parse::<u32>();

        if year_res.is_err() || doy_res.is_err() || hour_res.is_err() {
            continue;
        }

        let year = year_res.unwrap();
        let doy = doy_res.unwrap();
        let hour = hour_res.unwrap();

        // **Filter: September 2017**
        // Sept 4 is Day 247, Sept 11 is Day 254 (roughly).
        // We'll load the whole month to be safe: Days 244-273.
        if year != 2017 || !(240..=280).contains(&doy) {
            continue;
        }

        // Convert DOY to Date
        if let Some(date) = NaiveDate::from_yo_opt(year, doy)
            && let Some(dt) = date.and_hms_opt(hour, 0, 0)
        {
            if parts.len() < 8 {
                continue;
            }

            let kp_f64 = parts[5].parse::<f64>().unwrap_or(0.0);
            let dst_f64 = parts[6].parse::<f64>().unwrap_or(0.0);
            let f10_7_f64 = parts[7].parse::<f64>().unwrap_or(0.0);

            let kp: R = R::from(kp_f64);
            let dst: R = R::from(dst_f64);
            let f10_7: R = R::from(f10_7_f64);

            records.push(OmniRecord::new(dt.and_utc().timestamp(), dst, kp, f10_7));
        }
    }

    println!(
        "  Loaded {} Magnetic Ground Truth records (Sept 2017).",
        records.len()
    );
    Ok(records)
}
