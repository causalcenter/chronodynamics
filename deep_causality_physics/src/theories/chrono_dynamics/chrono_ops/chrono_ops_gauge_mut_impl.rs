// ============================================================================
// ChronoGaugeMutOps Implementation
// ============================================================================

use crate::theories::chrono_dynamics::chrono_ops::chrono_utils;
use crate::{ChronoGauge, ChronoOpsGaugeMut, SpaceTimeCoordinate};
use deep_causality_num::{Complex, RealField};
use deep_causality_topology::{CWComplex, LinkVariable, SU2_U1, TopologyError};
use std::fmt::Debug;
use std::iter::Sum;

impl<R> ChronoOpsGaugeMut<R> for ChronoGauge<R>
where
    R: RealField
        + Clone
        + From<f64>
        + Into<f64>
        + Default
        + Debug
        + Send
        + Sync
        + Sum
        + deep_causality_num::FromPrimitive
        + deep_causality_num::ToPrimitive,
{
    fn populate_links_from_source(&mut self) -> Result<(), TopologyError> {
        let data: &Vec<SpaceTimeCoordinate<R>> = self.source();
        let shape = self.lattice().shape();
        let n_radial = shape[1];
        let n_temporal = shape[0];

        if data.is_empty() {
            return Err(TopologyError::LatticeGaugeError(
                "Source data is empty, cannot populate links".to_string(),
            ));
        }

        // Bin data by radial index
        let mut radial_bins: Vec<Vec<R>> = vec![Vec::new(); n_radial];
        for coord in data {
            let r_idx = chrono_utils::radius_to_lattice_index(coord.r_m, n_radial);
            radial_bins[r_idx].push(coord.clock_drift_rate);
        }

        // Compute average clock drift per radial shell
        let mut avg_drifts: Vec<Option<R>> = Vec::with_capacity(n_radial);
        for bin in &radial_bins {
            if bin.is_empty() {
                avg_drifts.push(None);
            } else {
                let sum: R = bin.iter().cloned().fold(R::zero(), |a, b| a + b);
                let avg = sum / <R as From<f64>>::from(bin.len() as f64);
                avg_drifts.push(Some(avg));
            }
        }

        // Encoding scale: phase = (1 - drift) * N_t
        // This ensures Polyakov loop |P| = exp(-sum_phase) encodes the integrated potential
        let n_t_scale = <R as From<f64>>::from(n_temporal as f64);

        // Update links at each radial shell
        // Collect all cells (Temporal and Spatial)
        let all_cells: Vec<_> = self.lattice().cells(1).collect();

        // Prepare velocity bins if available, otherwise assume 0
        let mut radial_vel_bins: Vec<Vec<R>> = vec![Vec::new(); n_radial];
        for coord in data {
            let r_idx = chrono_utils::radius_to_lattice_index(coord.r_m, n_radial);
            radial_vel_bins[r_idx].push(coord.v_ms);
        }
        
        let c = <R as From<f64>>::from(3e8);
        let v_scale = n_t_scale / c;

        // Now iterate over collected cells and set links
        for cell in all_cells {
            let pos = cell.position();
            let r_idx = pos[1];
            let dir = cell.orientation().trailing_zeros() as usize;

            if let Some(avg_drift) = avg_drifts.get(r_idx).and_then(|opt| opt.as_ref()) {
                // Temporal Lines
                if dir == 0 {
                    let phase = -(*avg_drift) * n_t_scale;
                    let link = LinkVariable::<SU2_U1, Complex<R>, R>::from_phase(phase);
                    self.set_link(cell, link);
                } 
                // Spatial Links (Isotropic approximation)
                else {
                    // Compute average velocity for this shell
                    // (In a real implementation, we'd cache this like avg_drifts)
                    let bin = &radial_vel_bins[r_idx];
                    if !bin.is_empty() {
                         let sum: R = bin.iter().cloned().fold(R::zero(), |a, b| a + b);
                         let avg_vel = sum / <R as From<f64>>::from(bin.len() as f64);
                         
                         let phase = avg_vel * v_scale;
                         let link = LinkVariable::<SU2_U1, Complex<R>, R>::from_phase(phase);
                         self.set_link(cell, link);
                    }
                }
            }
        }

        Ok(())
    }

    fn populate_smooth_links_from_source(&mut self) -> Result<(), TopologyError> {
        let data: &Vec<SpaceTimeCoordinate<R>> = self.source();
        let shape = self.lattice().shape();
        let n_radial = shape[1];
        let n_temporal = shape[0];

        if data.is_empty() {
            return Err(TopologyError::LatticeGaugeError(
                "Source data is empty, cannot populate links".to_string(),
            ));
        }

        // 1. Prepare Data: Collect (r, drift, velocity) and sort by radius
        // We use the magnitude of velocity v_ms as a proxy, or components if available.
        // Assuming spherical symmetry for now, we encode v_ms into the spatial links uniformly.
        // Better would be to project onto lattice directions, but v_ms is a scalar in the struct.
        let mut points: Vec<(R, R, R)> = data
            .iter()
            .map(|coord| (coord.r_m, coord.clock_drift_rate, coord.v_ms))
            .collect();

        // Sort by radius
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // 2. Encoding scales
        let n_t_scale = <R as From<f64>>::from(n_temporal as f64);
        // Velocity scale: Normalize relative to c for phase encoding?
        // Or simple linear scaling? v is small compared to c.
        // Let's use v_scale = v / c * N_t to be consistent with time dilation scaling?
        // Actually, Virial relation is v^2 ~ phi.
        // Let's encode v directly into phase angle.
        // For consistency such that S_B ~ v^2, we need U_spatial = exp(i * v * scale).
        // Let's use scale = 1.0 / c (dimensionless v/c) * N_t.
        let c = <R as From<f64>>::from(3e8); // Approx C
        let v_scale = n_t_scale / c;

        // 3. Pre-calculate smoothed values
        let mut smooth_drifts = Vec::with_capacity(n_radial);
        let mut smooth_velocities = Vec::with_capacity(n_radial);

        // Helper vectors for interpolation
        let drift_points: Vec<(R, R)> = points.iter().map(|(r, d, _)| (*r, *d)).collect();
        let velocity_points: Vec<(R, R)> = points.iter().map(|(r, _, v)| (*r, *v)).collect();

        for r_idx in 0..n_radial {
            let r = chrono_utils::lattice_index_to_radius(r_idx, n_radial);

            let drift = chrono_utils::linear_interpolate_sorted(&drift_points, r);
            let vel = chrono_utils::linear_interpolate_sorted(&velocity_points, r);

            smooth_drifts.push(drift);
            smooth_velocities.push(vel);
        }

        // 4. Update Links
        // Collect all cells (Temporal and Spatial)
        let all_cells: Vec<_> = self.lattice().cells(1).collect();

        for cell in all_cells {
            let pos = cell.position();
            let r_idx = pos[1];
            let dir = cell.orientation().trailing_zeros() as usize;

            if r_idx >= n_radial {
                continue;
            } // Safety

            // Temporal Link (dir=0) -> Potential
            if dir == 0 {
                let drift = smooth_drifts[r_idx];
                let phase = -drift * n_t_scale;
                let link = LinkVariable::<SU2_U1, Complex<R>, R>::from_phase(phase);
                self.set_link(cell, link);
            }
            // Spatial Links (dir=1,2,3) -> Kinetic
            else {
                let v = smooth_velocities[r_idx];
                // Distribute v across spatial dimensions?
                // Alternatively, isotropic encoding: each spatial link gets v/sqrt(3)?
                // Or just encode v magnitude into all spatial links?
                // Let's try isotropic: phase ~ v.
                let phase = v * v_scale;
                let link = LinkVariable::<SU2_U1, Complex<R>, R>::from_phase(phase);
                self.set_link(cell, link);
            }
        }

        Ok(())
    }
}
