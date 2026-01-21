//! Save detailed analysis results to CSV.

use crate::types::gravity_types::GravitySolverResult;
use deep_causality_physics::EARTH_GM;
use std::fs::{self, File};
use std::io::{self, Write};

/// Saves detailed analysis results to a CSV file.
/// Returns the vector of GM values for further processing.
pub fn save_detailed_analysis(
    year_out_path: &str,
    dataset_name: &str,
    results: &[(u64, f64, f64, f64, GravitySolverResult<f64>)],
) -> io::Result<()> {
    fs::create_dir_all(year_out_path)?;

    let mut file = File::create(format!(
        "{}/{}_detailed_analysis.csv",
        year_out_path, dataset_name
    ))?;

    writeln!(
        file,
        "Timestamp,R_A,R_B,Delta_H_km,GM,Error_Pct,Term_Time,Term_Kinetic,Term_Potential"
    )?;

    for (timestamp, r_a, r_b, d_h_km, res) in results {
        let error = (res.gm - EARTH_GM).abs() / EARTH_GM * 100.0;
        writeln!(
            file,
            "{},{},{},{:.3},{:.5e},{:.4},{:.5e},{:.5e},{:.5e}",
            timestamp,
            r_a,
            r_b,
            d_h_km,
            res.gm,
            error,
            res.term_time,
            res.term_kinetic,
            res.term_potential
        )?;
    }

    Ok(())
}
