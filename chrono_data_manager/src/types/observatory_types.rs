// SPDX-License-Identifier: MIT
// Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.

use deep_causality_num::RealField;
use std::collections::HashMap;

/// State of a single satellite at an epoch.
/// Contains clock bias and position information for relativistic analysis.
#[derive(Debug, Clone)]
pub struct SatelliteState {
    pub sat_id: String,
    pub clock_bias_ns: f64, // Clock bias in nanoseconds
    pub x_m: f64,           // ECEF X position in meters
    pub y_m: f64,           // ECEF Y position in meters
    pub z_m: f64,           // ECEF Z position in meters
    pub radius_m: f64,      // Distance from Earth center in meters
}

impl SatelliteState {
    /// Create a new SatelliteState from position and clock data.
    pub fn new(sat_id: String, clock_bias_ns: f64, x_m: f64, y_m: f64, z_m: f64) -> Self {
        let radius_m = (x_m * x_m + y_m * y_m + z_m * z_m).sqrt();
        Self {
            sat_id,
            clock_bias_ns,
            x_m,
            y_m,
            z_m,
            radius_m,
        }
    }
}

/// Multi-satellite observation at a single epoch.
/// Groups satellite states by satellite ID at a common timestamp.
#[derive(Debug, Clone)]
pub struct ObservatoryEpoch {
    pub timestamp: i64, // Unix timestamp
    pub satellites: HashMap<String, SatelliteState>,
}

impl ObservatoryEpoch {
    /// Create a new empty epoch at the given timestamp.
    pub fn new(timestamp: i64) -> Self {
        Self {
            timestamp,
            satellites: HashMap::new(),
        }
    }

    /// Get a satellite state by ID.
    pub fn get(&self, sat_id: &str) -> Option<&SatelliteState> {
        self.satellites.get(sat_id)
    }

    /// Insert a satellite state.
    pub fn insert(&mut self, state: SatelliteState) {
        self.satellites.insert(state.sat_id.clone(), state);
    }
}

/// Aligned fleet data ready for analysis.
pub struct FleetData<R>
where
    R: RealField,
{
    pub satellites: Vec<String>,
    pub epochs: Vec<ObservatoryEpoch>,
    pub timestamp_to_index: HashMap<i64, usize>,
    pub matrix: Vec<Vec<R>>, // [n_epochs][n_sats] clock biases
}
