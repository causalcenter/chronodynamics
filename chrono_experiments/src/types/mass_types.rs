//! Mass experiment types.
//!
//! Generic over `R: RealField` to support both `f64` and `DoubleFloat`.

use crate::types::gravity_types::GravitySolverResult;

/// Detailed record of mass experiment results including full vectors.
/// Tuple structure: (Timestamp, R_A, R_B, Delta_H, Result, Position, Velocity)
pub type DetailedMassRecord<R> = (
    u64,                    // Timestamp
    R,                      // R_A
    R,                      // R_B
    R,                      // Delta H
    GravitySolverResult<R>, // Result
    [R; 3],                 // Position [x,y,z]
    [R; 3],                 // Velocity [vx,vy,vz]
);

/// Configuration for mass analysis.
///
/// Note: Config values remain f64 as they are thresholds/parameters, not data.
pub struct AnalysisConfig {
    pub window_size_indices: usize, // e.g., 240 (2 hours)
    pub min_height_diff_m: f64,     // e.g., 2,000,000.0
    pub outlier_sigma: f64,         // e.g., 3.0 (3-sigma rejection)
    pub step_size: usize,           // Step size to reduce autocorrelation
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            window_size_indices: 240,       // 2 hours (240 * 30s)
            min_height_diff_m: 2_000_000.0, // 2000 km
            outlier_sigma: 3.0,
            step_size: 1, // Keep 1 for high-res plots, use higher for stats
        }
    }
}
