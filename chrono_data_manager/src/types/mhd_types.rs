//! MHD-related data types for GQCD Chrono-MHD experiments.

use deep_causality_num::RealField;

/// OMNI magnetic field data record
#[derive(Debug, Clone)]
pub struct OmniRecord<R>
where
    R: RealField,
{
    timestamp: i64, // Unix Epoch Seconds
    dst_index: R,   // Magnetic Storm Intensity (The Signal)
    kp_index: R,    // Global Geomagnetic Activity
    f10_7_index: R, // Solar Radio Flux (The Control / Ionosphere)
}

impl<R> OmniRecord<R>
where
    R: RealField,
{
    pub fn new(timestamp: i64, dst_index: R, kp_index: R, f10_7_index: R) -> Self {
        Self {
            timestamp,
            dst_index,
            kp_index,
            f10_7_index,
        }
    }
}

impl<R> OmniRecord<R>
where
    R: RealField + Clone,
{
    pub fn timestamp(&self) -> i64 {
        self.timestamp
    }

    pub fn dst_index(&self) -> R {
        self.dst_index
    }

    pub fn kp_index(&self) -> R {
        self.kp_index
    }

    pub fn f10_7_index(&self) -> R {
        self.f10_7_index
    }
}

/// Hourly swarm statistics for MHD analysis
#[derive(Debug, Clone)]
pub struct HourlyStats<R>
where
    R: RealField,
{
    timestamp: i64,   // Unix timestamp of the hour start
    variance_ns2: R,  // Median variance of the swarm for this hour
    sat_count: usize, // Number of satellites contributing
}

impl<R> HourlyStats<R>
where
    R: RealField,
{
    pub fn new(timestamp: i64, variance_ns2: R, sat_count: usize) -> Self {
        Self {
            timestamp,
            variance_ns2,
            sat_count,
        }
    }
}

impl<R> HourlyStats<R>
where
    R: RealField + Clone,
{
    pub fn timestamp(&self) -> i64 {
        self.timestamp
    }

    pub fn variance_ns2(&self) -> R {
        self.variance_ns2
    }

    pub fn sat_count(&self) -> usize {
        self.sat_count
    }
}

/// Scenario configuration for MHD simulation
#[derive(Debug, Clone)]
pub struct MhdScenario<R>
where
    R: RealField,
{
    name: String,
    b_field_tesla: R,
    volume_m3: R,  // Characteristic volume of the field
    distance_m: R, // Distance from field center to clock
}

impl<R> MhdScenario<R>
where
    R: RealField,
{
    pub fn new(name: impl Into<String>, b_field_tesla: R, volume_m3: R, distance_m: R) -> Self {
        Self {
            name: name.into(),
            b_field_tesla,
            volume_m3,
            distance_m,
        }
    }
}

impl<R> MhdScenario<R>
where
    R: RealField + Clone,
{
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn b_field_tesla(&self) -> R {
        self.b_field_tesla
    }

    pub fn volume_m3(&self) -> R {
        self.volume_m3
    }

    pub fn distance_m(&self) -> R {
        self.distance_m
    }
}
