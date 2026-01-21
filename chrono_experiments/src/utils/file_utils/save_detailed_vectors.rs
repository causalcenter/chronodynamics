//! Save detailed vectors analysis to CSV.

use crate::types::mass_types::DetailedMassRecord;
use deep_causality_physics::EARTH_GM;
use std::fs::{self, File};
use std::io::{self, Write};

/// Saves detailed vector analysis results to a CSV file.
/// Includes full 3D position and velocity vectors.
pub fn save_detailed_vectors(
    year_out_path: &str,
    dataset_name: &str,
    results: &[DetailedMassRecord<f64>],
) -> io::Result<()> {
    let final_path = year_out_path;

    fs::create_dir_all(final_path)?;

    let mut file = File::create(format!(
        "{}/{}_detailed_vectors.csv",
        final_path, dataset_name
    ))?;

    // Header with all vector components
    writeln!(
        file,
        "Timestamp,R_A,R_B,Delta_H_km,GM,Error_Pct,Term_Time,Term_Kinetic,Term_Potential,X,Y,Z,VX,VY,VZ"
    )?;

    for (timestamp, r_a, r_b, d_h_km, res, pos, vel) in results {
        let error = (res.gm - EARTH_GM).abs() / EARTH_GM * 100.0;
        writeln!(
            file,
            "{},{},{},{:.3},{:.5e},{:.4},{:.5e},{:.5e},{:.5e},{:.4},{:.4},{:.4},{:.5},{:.5},{:.5}",
            timestamp,
            r_a,
            r_b,
            d_h_km,
            res.gm,
            error,
            res.term_time,
            res.term_kinetic,
            res.term_potential,
            pos[0],
            pos[1],
            pos[2],
            vel[0],
            vel[1],
            vel[2]
        )?;
    }

    Ok(())
}
