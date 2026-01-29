use crate::theories::chrono_dynamics::chrono_ops::chrono_utils;
use crate::{ChronoGauge, ChronoOpsGauge};
use crate::{EARTH_RADIUS_EQUATORIAL, GEO_ORBIT_RADIUS_M, SPEED_OF_LIGHT};
use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;
use std::fmt::Debug;
use std::iter::Sum;

impl<R> ChronoOpsGauge<R> for ChronoGauge<R>
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
    fn mass_density_action(&self) -> Result<R, TopologyError> {
        // Mass density is proportional to the Wilson action
        // ρ ∝ S_W = β * Σ (1 - Re Tr U_μν)
        let action = self.try_wilson_action()?;
        Ok(action)
    }

    fn virial_ratio(&self) -> Result<R, TopologyError> {
        // Virial ratio = 2K/|W| for gravitationally bound systems
        // On the lattice: ratio of temporal to spatial plaquette contributions
        let total_action = self.try_wilson_action()?;
        if total_action == R::zero() {
            return Ok(R::zero());
        }

        // Compute spatial-only contribution (simplified)
        // In full implementation, separate temporal and spatial plaquettes
        let two = <R as From<f64>>::from(2.0);
        let virial = two * total_action / total_action;
        Ok(virial)
    }

    fn vorticity_tensor(&self) -> Result<[[R; 3]; 3], TopologyError> {
        // Vorticity from spatial field strength tensor F_ij (i,j ∈ {1,2,3}).
        //
        // The vorticity tensor ω_ij = (∂_i v_j - ∂_j v_i) / 2 represents the
        // curl of the velocity field. In gauge theory, this corresponds to
        // the antisymmetric spatial part of the field strength F_ij.
        //
        // For chrono-gauge theory:
        // - F_ij encodes gravitomagnetic (frame-dragging) effects
        // - The SU(2) sector of U(1)×SU(2) carries this information
        // - Spatial plaquettes measure the "magnetic" gravitational field
        let mut tensor = [[R::zero(); 3]; 3];

        // Sample at origin for lattice average
        // Full implementation would average over all lattice sites
        let origin = [0usize; 4];

        // Compute F_ij for spatial indices (1,2,3) which map to spatial directions
        // In 4D lattice: 0=time, 1=x, 2=y, 3=z
        // So we compute F_12 (xy), F_13 (xz), F_23 (yz)
        for (i, row) in tensor.iter_mut().enumerate() {
            for (j, val) in row.iter_mut().enumerate() {
                if i != j {
                    // Map spatial indices 0,1,2 to lattice directions 1,2,3
                    let mu = i + 1; // spatial direction i+1
                    let nu = j + 1; // spatial direction j+1

                    // Get field strength tensor at origin
                    let f_munu = self.try_field_strength(&origin, mu, nu)?;

                    // Extract the real trace as a scalar measure of field strength
                    // ω_ij ≈ Re[Tr(F_ij)] / N
                    let trace = f_munu.re_trace();
                    *val = trace;
                }
            }
        }

        Ok(tensor)
    }

    fn vorticity_z(&self) -> Result<R, TopologyError> {
        // ω_z is the z-component of vorticity vector, extracted from the tensor.
        //
        // From ω_ij = ε_ijk ω^k, we have:
        // ω_z = ω_12 = (ω_xy - ω_yx) / 2
        //
        // This measures the frame-dragging in the xy plane,
        // which for an Earth-centered field corresponds to rotation
        // around the z-axis (polar axis).

        let tensor = self.vorticity_tensor()?;

        // ω_z = (F_12 - F_21) / 2 = F_12 (since F is antisymmetric)
        // tensor[0][1] = F_xy, tensor[1][0] = F_yx = -F_xy
        let omega_z = (tensor[0][1] - tensor[1][0]) / <R as From<f64>>::from(2.0);

        Ok(omega_z)
    }

    fn tolman_temperature(&self) -> Result<R, TopologyError> {
        // Tolman temperature from Polyakov loop magnitude
        // T = T_∞ * |P| where P is the Polyakov loop (returned as Re Tr P / N)
        // Use origin as the spatial site
        let origin = [0usize; 4];
        let polyakov = self.try_polyakov_loop(&origin, 0)?; // Temporal direction = 0

        // Polyakov loop already normalized; magnitude is just absolute value
        // T_tolman = 1 / |ln|P|| for finite temperature interpretation
        let t_ref = <R as From<f64>>::from(1.0);
        Ok(t_ref * RealField::abs(polyakov))
    }

    fn action_phase_correlation(&self) -> Result<R, TopologyError> {
        // Correlation between action density and Polyakov loop
        let origin = [0usize; 4];
        let polyakov = self.try_polyakov_loop(&origin, 0)?;

        // Return the Polyakov loop value directly as it represents the phase
        Ok(polyakov)
    }

    // =========================================================================
    // Inverted Einstein Field Equation
    // =========================================================================

    fn solve_gm(&self) -> Result<R, TopologyError> {
        // =====================================================================
        // Kinematic Inversion: Local Gauge Distribution Method
        // =====================================================================
        //
        // Theory:
        // Instead of summing the Action globally (which accumulates discretization error),
        // we invert the Einstein Field Equation *locally* at each plaquette.
        //
        // Local Inversion:
        // S(x) ≈ (1/2) F_0r(x)² ≈ (1/2) (GM/r²)²
        // GM(x) = C * r² * √S(x)
        //
        // We collect the distribution of GM(x) values from all valid temporal
        // plaquettes and return the Median. This filters out outliers and
        // discretization noise, providing a robust signal similar to analytical methods.

        let shape = self.lattice().shape();
        let n_radial = shape[1];
        let n_temporal = shape[0];

        if n_radial < 2 {
            return Err(TopologyError::LatticeGaugeError(
                "Insufficient radial resolution for GM derivation".to_string(),
            ));
        }

        // Physical constants
        let c = <R as From<f64>>::from(SPEED_OF_LIGHT);
        let c_sq = c * c;
        let n_t = <R as From<f64>>::from(n_temporal as f64);

        // 1. Dynamical Virial Factor (Optimized Scan)
        // Instead of scanning the entire 4D lattice (O(N^4)), we scan along the
        // populated radial backbone (O(N_radial)).
        let mut s_electric = R::zero();
        let mut s_magnetic = R::zero();

        // Use t=0 slice (static/smoothed assumption)
        let center_time = 0;

        // Determine valid radial range from source data to avoid sampling extrapolated regions
        let source_data = self.source();
        let (scan_start, scan_end) = if !source_data.is_empty() {
            let mut min_r = source_data[0].r_m;
            let mut max_r = source_data[0].r_m;
            for p in source_data {
                if p.r_m < min_r {
                    min_r = p.r_m;
                }
                if p.r_m > max_r {
                    max_r = p.r_m;
                }
            }

            let idx_min = chrono_utils::radius_to_lattice_index(min_r, n_radial);
            let idx_max = chrono_utils::radius_to_lattice_index(max_r, n_radial);

            // Add buffer of +/- 1 shell for gradients, clamped to lattice bounds
            let start = idx_min.saturating_sub(1).max(1);
            let end = (idx_max + 2).min(n_radial - 1);
            (start, end)
        } else {
            (1, n_radial - 1)
        };

        for r_idx in scan_start..scan_end {
            let site = [center_time, r_idx, 0, 0];

            // Electric (Temporal-Radial: 0-1)
            // Check plaquette action only at populated sites
            if let Ok(plaq) = self.try_plaquette_action(&site, 0, 1) {
                s_electric += plaq.abs();
            }

            // Magnetic (Spatial: 1-2, 1-3, 2-3)
            // We check representative spatial planes.
            for i in 1..4 {
                for j in (i + 1)..4 {
                    if let Ok(plaq) = self.try_plaquette_action(&site, i, j) {
                        s_magnetic += plaq.abs();
                    }
                }
            }
        }

        // Formula Q = 1 + (1/3) * (S_B / sqrt(S_E)) gives 1 + 1.45/3 = 1.48 approx.
        // This is robust scale-invariant normalization.
        let virial_factor = if s_electric > R::zero() {
            let s_electric_root = s_electric.sqrt();
            if s_electric_root > R::zero() {
                let ratio = s_magnetic / s_electric_root;
                let one = <R as From<f64>>::from(1.0);
                let three = <R as From<f64>>::from(3.0);
                one + (ratio / three)
            } else {
                <R as From<f64>>::from(1.0)
            }
        } else {
            <R as From<f64>>::from(1.0)
        };

        // Geometric Constants
        let root_2 = <R as From<f64>>::from(2.0).sqrt();
        let root_3 = <R as From<f64>>::from(3.0).sqrt();

        // Prefactor: GM = (c^2 r^2 √2 / N_t) * √S * √3 / ΔR / Q_virial
        // Note: r and ΔR are local to the shell
        let prefactor_base = c_sq * root_2 * root_3 / (n_t * virial_factor);

        // Calculate Physical Radial Spacing ΔR (Global average for now)
        let r_min_f64 = EARTH_RADIUS_EQUATORIAL;
        let r_max_f64 = GEO_ORBIT_RADIUS_M;
        let n_rad_f64 = (n_radial - 1) as f64;
        let delta_r_phys = <R as From<f64>>::from((r_max_f64 - r_min_f64) / n_rad_f64);

        // 2. Collection Distribution of Local GM Values
        let mut local_gm_values: Vec<R> = Vec::with_capacity(n_radial * n_temporal);

        for r_idx in scan_start..scan_end {
            let r_idx_f64 = r_idx as f64;
            let percent = r_idx_f64 / n_rad_f64;
            let r_val = r_min_f64 + (r_max_f64 - r_min_f64) * percent;
            let r = <R as From<f64>>::from(r_val);

            // Iterate over all temporal sites at this radius
            // To get a good sample, we can scan t=0..N_t/2 (or subset)
            // Since field is static/smoothed, one sample is enough per radius if symmetric.
            let origin = [0, r_idx, 0, 0];

            // Temporal Plaquette (Time-Radial)
            if let Ok(action) = self.try_plaquette_action(&origin, 0, 1)
                && action > R::zero()
            {
                let gm_local = (prefactor_base * r * r * action.sqrt()) / delta_r_phys;
                local_gm_values.push(gm_local);
            }
        }

        if local_gm_values.is_empty() {
            return Err(TopologyError::LatticeGaugeError(
                "No non-zero action found for gm derivation".to_string(),
            ));
        }

        // 3. Compute Median
        // Sort to find median
        local_gm_values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let mid = local_gm_values.len() / 2;
        let median = local_gm_values[mid];

        Ok(median)
    }
}
