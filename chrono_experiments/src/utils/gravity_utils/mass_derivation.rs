/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::utils::proces_utils::{apply_mad_filter, interpolate_space_time};
use crate::{AnalysisConfig, SpaceTimeCoordinate};
use chrono_data_manager::DataManager;
use deep_causality_num::RealField;
use deep_causality_physics::{ChronoGauge, ChronoGaugeOps};
use deep_causality_physics::{
    EARTH_GM, EARTH_MASS_KG, EARTH_RADIUS, NEWTONIAN_CONSTANT_OF_GRAVITATION,
};
use std::io;

/// Extract GPS week from dataset name (e.g., "gbm19670" -> 1967)
pub fn extract_gps_week(dataset: &str) -> Option<u32> {
    if dataset.len() >= 7 && dataset.starts_with("gbm") {
        let week_str = &dataset[3..7];
        week_str.parse::<u32>().ok()
    } else {
        None
    }
}

/// Result of a mass derivation for a single dataset
#[derive(Debug, Clone)]
pub struct MassDerivationResult<R> {
    pub dataset_name: String,
    pub derived_gm: R,
    pub derived_mass: R,
    pub derived_gravity: R,
    pub gm_error_percent: R,
    pub mass_error_percent: R,
    pub filtered_values: Vec<R>,
}

/// Processes a single dataset to derive GM, Mass, and Gravity.
///
/// Uses `interpolate_space_time_two_pass` for high-precision mass retrieval.
/// Filters out anomalous weeks (1967, 1968).
pub fn derive_mass_from_dataset<R>(
    dataset_name: &str,
    clk_path: &str,
    sp3_path: &str,
    sat_id: &str,
    gauge_field: &ChronoGauge<R>,
    anomalous_weeks: &[u32],
) -> io::Result<Option<MassDerivationResult<R>>>
where
    R: RealField
        + Clone
        + From<f64>
        + Into<f64>
        + RealField
        + Default
        + std::fmt::Debug
        + deep_causality_num::FromPrimitive
        + deep_causality_num::ToPrimitive,
{
    // 1. Filter Anomalous Weeks
    if let Some(week) = extract_gps_week(dataset_name)
        && anomalous_weeks.contains(&week)
    {
        return Ok(None);
    }

    let config = AnalysisConfig::default();
    let dm = DataManager::default();

    // 2. Load Data
    let (clocks, orbits) = dm.load_gnss_single_satellite(clk_path, sp3_path, sat_id)?;
    if clocks.is_empty() || orbits.len() < 10 {
        return Ok(None);
    }

    // 3. Interpolate (Two Pass Filter)
    let data: Vec<SpaceTimeCoordinate<R>> = interpolate_space_time(&clocks, &orbits);

    if data.len() <= config.window_size_indices {
        return Ok(None);
    }

    // 4. Derive GM values using ChronoGaugeWitness
    let mut raw_gm_values: Vec<R> = Vec::new();

    let mut i = 0;
    while i < data.len() - config.window_size_indices {
        let idx_a = i;
        let idx_b = i + config.window_size_indices;

        // We need SpaceTimeCoordinate to implement ChronoGaugeWitness logic if not using the trait...
        // But assuming the trait is available as used in E01.
        if let Ok(gm) = gauge_field.solve_gm(&data[idx_a], &data[idx_b]) {
            raw_gm_values.push(gm);
        }

        i += config.step_size;
    }

    if raw_gm_values.is_empty() {
        return Ok(None);
    }

    // 5. Apply MAD Filter
    let filtered_gm =
        apply_mad_filter(&raw_gm_values, <R as From<f64>>::from(config.outlier_sigma));
    if filtered_gm.is_empty() {
        return Ok(None);
    }

    // 6. Calculate Means
    let n = <R as From<f64>>::from(filtered_gm.len() as f64);
    let derived_gm: R = filtered_gm.iter().fold(R::zero(), |acc, &x| acc + x) / n;

    // 7. Derive Mass and Gravity
    // Earth Mass ~ 5.9722e24 kg
    let mass_ref = <R as From<f64>>::from(EARTH_MASS_KG);

    // G constant 6.674_30e-11 m^3 kg^-1 s^-2
    let g_ref = <R as From<f64>>::from(NEWTONIAN_CONSTANT_OF_GRAVITATION);

    // Derived Mass: Mass = GM / G
    let derived_mass = derived_gm / g_ref;

    // Derived Gravity: g = GM / R_earth^2
    // R_earth approx 6371 km
    let r_earth = <R as From<f64>>::from(EARTH_RADIUS);
    let derived_gravity = derived_gm / (r_earth * r_earth);

    // Errors
    // GM Ref = 3.986004418e14
    let gm_ref = <R as From<f64>>::from(EARTH_GM);

    // Error calculation: |(val - ref) / ref| * 100
    let one_hundred = <R as From<f64>>::from(100.0);

    let gm_error_percent = RealField::abs(derived_gm - gm_ref) / gm_ref * one_hundred;
    let mass_error_percent = RealField::abs(derived_mass - mass_ref) / mass_ref * one_hundred;

    Ok(Some(MassDerivationResult {
        dataset_name: dataset_name.to_string(),
        derived_gm,
        derived_mass,
        derived_gravity,
        gm_error_percent,
        mass_error_percent,
        filtered_values: filtered_gm,
    }))
}
