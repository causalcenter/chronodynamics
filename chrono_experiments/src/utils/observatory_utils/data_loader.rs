/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Fleet data loading - exact port of legacy gqcd/src/utils/chrono_gravity_utils/merge_satellite_data

use std::collections::HashMap;
use std::path::Path;
use std::sync::Mutex;

use chrono_data_manager::{
    ClockData, DataManager, FleetData, ObservatoryEpoch, OrbitData, SatelliteState,
    get_year_datasets,
};
use deep_causality_num::RealField;
use rayon::prelude::*;

/// Load fleet data for a specific year.
/// Direct port of legacy gqcd implementation.
pub fn load_fleet_data<R, P>(
    dm: &DataManager,
    year_path: P,
    satellite_ids: &[&str],
) -> Option<FleetData<R>>
where
    R: RealField + From<f64> + Into<f64> + Send + Sync + 'static + Copy + Clone,
    P: AsRef<Path> + Sync,
{
    let year_str = year_path.as_ref().file_name()?.to_str()?;
    let datasets = get_year_datasets(year_str);

    if datasets.is_empty() {
        return None;
    }

    // Thread-safe epoch accumulator
    let all_epochs: Mutex<Vec<ObservatoryEpoch>> = Mutex::new(Vec::new());

    // Parallel processing of datasets
    datasets.par_iter().for_each(|dataset| {
        let clk_path = year_path.as_ref().join(format!("{}.clk", dataset));
        let sp3_path = year_path.as_ref().join(format!("{}.sp3", dataset));

        if !clk_path.exists() {
            return;
        }

        // Collect epochs for this dataset
        let mut local_epochs: Vec<ObservatoryEpoch> = Vec::new();

        // Load data for each satellite
        for sat_id in satellite_ids {
            if let Ok((clocks, orbits)) =
                dm.load_gnss_single_satellite::<R, _>(&clk_path, &sp3_path, sat_id)
            {
                // EXACT port of legacy merge_satellite_data
                merge_satellite_data(&mut local_epochs, sat_id, &clocks, &orbits);
            }
        }

        // Merge local epochs into global
        let mut global = all_epochs.lock().unwrap();
        for local in local_epochs {
            if let Some(existing) = global.iter_mut().find(|e| e.timestamp == local.timestamp) {
                existing.satellites.extend(local.satellites);
            } else {
                global.push(local);
            }
        }
    });

    let mut epochs = all_epochs.into_inner().unwrap();
    if epochs.is_empty() {
        return None;
    }

    // Sort by timestamp - NO intersection filter (matches legacy behavior)
    epochs.sort_by_key(|e| e.timestamp);

    // Build Matrix (use 0.0 for missing satellites - legacy behavior)
    // NOTE: This is physically questionable but matches the legacy implementation.
    let mut matrix: Vec<Vec<R>> = Vec::with_capacity(epochs.len());
    let mut timestamp_to_index: HashMap<i64, usize> = HashMap::new();

    for (idx, epoch) in epochs.iter().enumerate() {
        timestamp_to_index.insert(epoch.timestamp, idx);
        let mut row: Vec<R> = Vec::with_capacity(satellite_ids.len());
        for id in satellite_ids {
            let bias_ns = epoch
                .satellites
                .get(*id)
                .map(|s| s.clock_bias_ns)
                .unwrap_or(0.0);
            row.push(R::from(bias_ns));
        }
        matrix.push(row);
    }

    Some(FleetData {
        satellites: satellite_ids.iter().map(|s| s.to_string()).collect(),
        epochs,
        timestamp_to_index,
        matrix,
    })
}

/// EXACT port of legacy gqcd/src/utils/chrono_gravity_utils/mod.rs:180-230
/// merge_satellite_data function
pub fn merge_satellite_data<R>(
    epochs: &mut Vec<ObservatoryEpoch>,
    sat_id: &str,
    clocks: &[ClockData<R>],
    orbits: &[OrbitData<R>],
) where
    R: RealField + From<f64> + Into<f64> + Copy,
{
    // Create orbit lookup by timestamp (unix timestamp as i64)
    let orbit_map: HashMap<i64, &OrbitData<R>> = orbits
        .iter()
        .map(|o| (o.timestamp().and_utc().timestamp(), o))
        .collect();

    for clock in clocks {
        let ts = clock.timestamp().and_utc().timestamp();

        // Find matching orbit
        let orbit = match orbit_map.get(&ts) {
            Some(o) => *o,
            None => continue,
        };

        let x_m: f64 = orbit.x_m().into();
        let y_m: f64 = orbit.y_m().into();
        let z_m: f64 = orbit.z_m().into();
        // Calculation of radius using explicit multiplication instead of powi
        let _radius_m = (x_m * x_m + y_m * y_m + z_m * z_m).sqrt();

        // Convert clock bias from seconds to nanoseconds
        let bias_s: f64 = clock.bias_s().into();
        let clock_bias_ns = bias_s * 1e9;

        let state = SatelliteState::new(sat_id.to_string(), clock_bias_ns, x_m, y_m, z_m);

        // Find or create epoch
        if let Some(epoch) = epochs.iter_mut().find(|e| e.timestamp == ts) {
            epoch.insert(state);
        } else {
            let mut epoch = ObservatoryEpoch::new(ts);
            epoch.insert(state);
            epochs.push(epoch);
        }
    }
}
