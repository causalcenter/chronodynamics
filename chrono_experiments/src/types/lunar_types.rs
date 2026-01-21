/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

/// Statistics for a single day
#[derive(Debug, Clone)]
pub struct DayStats {
    day: usize,
    count: usize,
    sum: f64,
    sum_sq: f64,
    min: f64,
    max: f64,
}

impl DayStats {
    pub fn new(day: usize) -> Self {
        Self {
            day,
            count: 0,
            sum: 0.0,
            sum_sq: 0.0,
            min: f64::MAX,
            max: f64::MIN,
        }
    }

    pub fn add(&mut self, g: f64) {
        self.count += 1;
        self.sum += g;
        self.sum_sq += g * g;
        if g < self.min {
            self.min = g;
        }
        if g > self.max {
            self.max = g;
        }
    }

    pub fn mean(&self) -> f64 {
        if self.count > 0 {
            self.sum / self.count as f64
        } else {
            0.0
        }
    }

    pub fn std_dev(&self) -> f64 {
        if self.count > 1 {
            let mean = self.mean();
            let variance =
                (self.sum_sq - self.count as f64 * mean * mean) / (self.count - 1) as f64;
            variance.sqrt()
        } else {
            0.0
        }
    }

    pub fn deviation_pct(&self, g_ref: f64) -> f64 {
        (self.mean() - g_ref) / g_ref * 100.0
    }

    pub fn day(&self) -> usize {
        self.day
    }

    pub fn count(&self) -> usize {
        self.count
    }

    pub fn sum(&self) -> f64 {
        self.sum
    }

    pub fn sum_sq(&self) -> f64 {
        self.sum_sq
    }

    pub fn min(&self) -> f64 {
        self.min
    }

    pub fn max(&self) -> f64 {
        self.max
    }
}

/// Monthly analysis results
#[derive(Debug, Clone)]
pub struct MonthAnalysis {
    name: String,
    quarter: String,
    total_points: usize,
    days: Vec<DayStats>,
    overall_mean: f64,
    overall_std_dev: f64,
    peak_deviation: f64,
    trough_deviation: f64,
    amplitude: f64,
    peak_day: usize,
    trough_day: usize,
    // Tier A: Derived Lunar Distance (Kepler)
    detected_period_days: f64,
    derived_lunar_distance: f64,
    theoretical_lunar_distance: f64,
}

impl MonthAnalysis {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: String,
        quarter: String,
        total_points: usize,
        days: Vec<DayStats>,
        overall_mean: f64,
        overall_std_dev: f64,
        peak_deviation: f64,
        trough_deviation: f64,
        amplitude: f64,
        peak_day: usize,
        trough_day: usize,
        detected_period_days: f64,
        derived_lunar_distance: f64,
        theoretical_lunar_distance: f64,
    ) -> Self {
        Self {
            name,
            quarter,
            total_points,
            days,
            overall_mean,
            overall_std_dev,
            peak_deviation,
            trough_deviation,
            amplitude,
            peak_day,
            trough_day,
            detected_period_days,
            derived_lunar_distance,
            theoretical_lunar_distance,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn quarter(&self) -> &str {
        &self.quarter
    }

    pub fn total_points(&self) -> usize {
        self.total_points
    }

    pub fn days(&self) -> &Vec<DayStats> {
        &self.days
    }

    pub fn overall_mean(&self) -> f64 {
        self.overall_mean
    }

    pub fn overall_std_dev(&self) -> f64 {
        self.overall_std_dev
    }

    pub fn peak_deviation(&self) -> f64 {
        self.peak_deviation
    }

    pub fn trough_deviation(&self) -> f64 {
        self.trough_deviation
    }

    pub fn amplitude(&self) -> f64 {
        self.amplitude
    }

    pub fn peak_day(&self) -> usize {
        self.peak_day
    }

    pub fn trough_day(&self) -> usize {
        self.trough_day
    }

    pub fn detected_period_days(&self) -> f64 {
        self.detected_period_days
    }

    /// Derived lunar distance from clock periodicity (m)
    pub fn derived_lunar_distance(&self) -> f64 {
        self.derived_lunar_distance
    }

    /// Theoretical lunar distance (m)
    pub fn theoretical_lunar_distance(&self) -> f64 {
        self.theoretical_lunar_distance
    }

    /// Error of derived vs theoretical (%)
    pub fn distance_error_pct(&self) -> f64 {
        if self.theoretical_lunar_distance != 0.0 {
            ((self.derived_lunar_distance - self.theoretical_lunar_distance)
                / self.theoretical_lunar_distance)
                .abs()
                * 100.0
        } else {
            0.0
        }
    }
}
