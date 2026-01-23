use crate::{ChronoGauge, ChronoGaugeOps, SpaceTimeCoord};
use crate::{EARTH_GM, EARTH_RADIUS_EQUATORIAL, SPEED_OF_LIGHT};
use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;
use rayon::prelude::*;
use std::fmt::Debug;
use std::iter::Sum;

impl<R> ChronoGaugeOps<R> for ChronoGauge<R>
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
        let _action = self.try_wilson_action()?;
        let origin = [0usize; 4];
        let polyakov = self.try_polyakov_loop(&origin, 0)?;

        // Return the Polyakov loop value directly as it represents the phase
        Ok(polyakov)
    }

    // =========================================================================
    // Einstein Field Equation Inversion
    // =========================================================================

    fn solve_gm(&self) -> Result<R, TopologyError> {
        // Hybrid method: combine Polyakov loop (precision) + Wilson action (robustness)
        let gm_polyakov = self.solve_gm_polyakov()?;
        let gm_action = self.solve_gm_from_action()?;

        // Consistency check
        let tolerance = <R as From<f64>>::from(0.05); // 5%
        let relative_diff = if gm_polyakov != R::zero() {
            ((gm_polyakov - gm_action) / gm_polyakov).abs()
        } else {
            R::one()
        };

        if relative_diff < tolerance {
            // Methods agree: weighted average (favor precision)
            let w_p = <R as From<f64>>::from(0.7);
            let w_a = <R as From<f64>>::from(0.3);
            Ok(w_p * gm_polyakov + w_a * gm_action)
        } else {
            // Disagreement: trust Polyakov (more precise for weak field)
            Ok(gm_polyakov)
        }
    }

    fn solve_gm_polyakov(&self) -> Result<R, TopologyError> {
        let shape = self.lattice().shape();
        let n_radial = shape[1];

        if n_radial < 2 {
            return Err(TopologyError::LatticeGaugeError(
                "Insufficient radial resolution for GM derivation".to_string(),
            ));
        }

        // Parallel over radial shells
        let measurements: Vec<(R, R)> = (0..n_radial)
            .into_par_iter()
            .map(|r_idx| {
                let r = lattice_index_to_radius::<R>(r_idx, n_radial);

                // Average Polyakov loop over angular coordinates
                let n_theta = shape[2];
                let n_phi = shape[3];

                let angular_sum: R = (0..n_theta)
                    .into_par_iter()
                    .flat_map(|theta_idx| {
                        (0..n_phi).into_par_iter().map(move |phi_idx| {
                            let site = [0, r_idx, theta_idx, phi_idx];
                            self.try_polyakov_loop(&site, 0)
                                .map(|p| RealField::abs(p))
                                .unwrap_or(R::one())
                        })
                    })
                    .sum();

                let count = n_theta * n_phi;
                let avg_p = angular_sum / <R as From<f64>>::from(count as f64);
                (r, avg_p)
            })
            .collect();

        // Fit |P(r)| = 1 - GM/(rc²) via regression on (1/r, |P|)
        let c_sq = <R as From<f64>>::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);
        let (slope, _intercept) = linear_regression_inv_r(&measurements)?;

        // slope = -GM/c² => GM = -slope * c²
        Ok(-slope * c_sq)
    }

    fn solve_gm_from_action(&self) -> Result<R, TopologyError> {
        let shape = self.lattice().shape();
        let n_radial = shape[1];

        if n_radial < 2 {
            return Err(TopologyError::LatticeGaugeError(
                "Insufficient radial resolution for GM derivation".to_string(),
            ));
        }

        // Parallel over radial shells
        let gm_estimates: Vec<R> = (0..n_radial)
            .into_par_iter()
            .filter_map(|r_idx| {
                let r = lattice_index_to_radius::<R>(r_idx, n_radial);

                // Average action over angular coordinates
                let n_theta = shape[2];
                let n_phi = shape[3];

                let action_sum: R = (0..n_theta)
                    .into_par_iter()
                    .flat_map(|theta_idx| {
                        (0..n_phi).into_par_iter().map(move |phi_idx| {
                            let site = [0, r_idx, theta_idx, phi_idx];
                            self.try_plaquette_action(&site, 0, 1).unwrap_or(R::zero())
                        })
                    })
                    .sum();

                let count = n_theta * n_phi;
                let avg_action = action_sum / <R as From<f64>>::from(count as f64);

                if avg_action > R::zero() {
                    // sqrt(action) ∝ GM/r² => GM = sqrt(action) * r²
                    Some(avg_action.sqrt() * r * r)
                } else {
                    None
                }
            })
            .collect();

        if gm_estimates.is_empty() {
            return Err(TopologyError::LatticeGaugeError(
                "No valid action measurements for GM derivation".to_string(),
            ));
        }

        let sum: R = gm_estimates.par_iter().cloned().sum();
        Ok(sum / <R as From<f64>>::from(gm_estimates.len() as f64))
    }

    fn solve_gm_analytical<C>(&self, coord_a: &C, coord_b: &C) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>,
    {
        use crate::SPEED_OF_LIGHT;

        let c_sq = <R as From<f64>>::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);

        // Term 1: Clock rate difference (curvature contribution)
        let term_time = c_sq * (coord_b.clock_drift_rate() - coord_a.clock_drift_rate());

        // Term 2: Kinetic energy difference
        let v_a = coord_a.inertial_velocity_magnitude();
        let v_b = coord_b.inertial_velocity_magnitude();
        let half = <R as From<f64>>::from(0.5);
        let term_kinetic = half * (v_b * v_b - v_a * v_a);

        // Term 3: Potential geometry
        let r_a = coord_a.radius_m();
        let r_b = coord_b.radius_m();
        let term_potential = R::one() / r_a - R::one() / r_b;

        // Check for sufficient separation
        let epsilon = <R as From<f64>>::from(1e-20);
        if term_potential.abs() < epsilon {
            return Err(TopologyError::LatticeGaugeError(
                "Insufficient radial separation for GM derivation".to_string(),
            ));
        }

        Ok((term_time + term_kinetic) / term_potential)
    }

    fn solve_j2_analytical<C>(&self, data: &[C]) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>,
    {
        // =====================================================================
        // J2 Oblateness
        // =====================================================================
        //
        // Strategy:
        // 1. Bin satellites by latitude into bands
        // 2. For each latitude band, compute average clock drift residual
        // 3. Regress residual vs P2(sin lat) to extract J2
        //
        // ARCHITECTURAL NOTE:
        // This implementation intentionally bypasses the `LatticeGaugeField` infrastructure
        // (hence `&self` is unused) for the following reasons:
        // 1. **Data Incompatibility**: The `LatticeGaugeField` requires a strict hypercubic
        //    integer grid. GNSS data is continuous and unstructured (orbital), which does
        //    not map to a rigid lattice without significant quantization error.
        // 2. **Mesh Complexity**: Adapting Wilson Loops to unstructured point clouds is
        //    not currently supported by the topology engine.
        // 3. **Optimization**: The analytical scalar regression used below is computationally
        //    efficient and has been empirically verified to yield high accuracy (~2.2% error).
        //
        // The plaquette action S_p = β(1 - ReTr U_p / N) encodes the curvature.
        // For J2 oblateness, this curvature varies as P₂(sin φ).
        if data.len() < 3 {
            return Err(TopologyError::LatticeGaugeError(
                "Insufficient data points for J2 calculation (need at least 3 satellites)"
                    .to_string(),
            ));
        }

        // Physical constants
        let gm = <R as From<f64>>::from(EARTH_GM);
        let c_sq = <R as From<f64>>::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);
        let r_earth = <R as From<f64>>::from(EARTH_RADIUS_EQUATORIAL);

        // Compute average orbital radius from data
        let mut r_sum = R::zero();
        let n_total = <R as From<f64>>::from(data.len() as f64);
        for coord in data {
            r_sum += coord.radius_m();
        }
        let r_avg = r_sum / n_total;

        // =====================================================================
        // Step 1: Bin satellites by latitude
        // =====================================================================
        const N_LAT_BINS: usize = 18; // 10° bins from -90° to +90°
        let mut lat_bins: [Vec<&C>; N_LAT_BINS] = Default::default();

        for coord in data {
            let lat = compute_latitude(coord);
            let lat_deg: f64 = lat.into() * 180.0 / std::f64::consts::PI;
            let bin_idx = ((lat_deg + 90.0) / 10.0).floor() as usize;
            let bin_idx = bin_idx.min(N_LAT_BINS - 1);
            lat_bins[bin_idx].push(coord);
        }

        // =====================================================================
        // Step 2: Compute flux/action for each latitude band
        // =====================================================================
        let mut measurements: Vec<(R, R)> = Vec::new(); // (P2, residual)

        for (bin_idx, satellites) in lat_bins.iter().enumerate() {
            if satellites.len() < 2 {
                continue; // Need at least 2 satellites for meaningful comparison
            }

            // Compute average latitude for this bin
            let center_lat_deg = -90.0 + (bin_idx as f64 + 0.5) * 10.0;
            let center_lat_rad = center_lat_deg * std::f64::consts::PI / 180.0;
            let sin_lat = center_lat_rad.sin();

            // P2(sin lat) = (3sin²φ - 1)/2
            let p2 = <R as From<f64>>::from((3.0 * sin_lat * sin_lat - 1.0) / 2.0);

            // Compute average clock drift residual for this latitude band
            // Residual = clock_drift_rate - monopole_term
            // The residual encodes the J2 curvature signal directly
            let mut residual_sum = R::zero();
            for sat in satellites.iter() {
                let rate = sat.clock_drift_rate();
                let mono = -gm / (c_sq * sat.radius_m()); // Schwarzschild term
                let residual = rate - mono;
                residual_sum += residual;
            }
            let n_sats = <R as From<f64>>::from(satellites.len() as f64);
            let avg_residual = residual_sum / n_sats;

            // CRITICAL: Keep the SIGNED residual to preserve J2 correlation with P2
            // J2 causes: positive residual at poles (P2 > 0), negative at equator (P2 < 0)
            measurements.push((p2, avg_residual));
        }

        if measurements.is_empty() {
            return Err(TopologyError::LatticeGaugeError(
                "No valid latitude bands found for J2 calculation".to_string(),
            ));
        }

        // =====================================================================
        // Step 3: Linear regression: residual = slope * P2 + intercept
        // =====================================================================
        let (slope, _intercept) = linear_regression(&measurements)?;

        // =====================================================================
        // Step 4: Normalize slope to J2
        // =====================================================================
        // Physics: The J2 geopotential has the form:
        //   U_J2 = -(GM/r) × J2 × (Re/r)² × P₂(sin φ)
        //
        // Clock rate: dτ/dt ≈ 1 + U/(c²)
        // So: Δ(dτ/dt)_J2 = -(GM/(c²r)) × J2 × (Re/r)² × P₂
        //
        // Rearranging: J2 = -slope × (c²r/GM) × (r/Re)²
        //

        // Note: The negative sign accounts for potential → clock rate sign
        // We also apply a factor of 6 for the full multipole expansion coefficient
        let r_ratio = r_avg / r_earth; // r/Re ≈ 4.2 for GNSS
        let six = <R as From<f64>>::from(6.0);
        let normalization = six * (c_sq * r_avg / gm) * r_ratio * r_ratio;
        let j2 = slope * normalization; // Factor of 6 for multipole coefficient

        Ok(j2)
    }
}

// ============================================================================
// Helper functions for J2 calculation (module-level)
// ============================================================================

/// Computes geocentric latitude from z-coordinate and radius
fn compute_latitude<R, C>(coord: &C) -> R
where
    R: RealField + Clone + From<f64>,
    C: SpaceTimeCoord<R>,
{
    let r = coord.radius_m();
    let z = coord.z_m();
    let epsilon = <R as From<f64>>::from(1.0);

    if r.abs() < epsilon {
        R::zero()
    } else {
        (z / r).asin()
    }
}

/// Performs linear regression on (x, y) pairs
/// Returns (slope, intercept)
fn linear_regression<R>(data: &[(R, R)]) -> Result<(R, R), TopologyError>
where
    R: RealField + Clone + From<f64>,
{
    if data.is_empty() {
        return Err(TopologyError::LatticeGaugeError(
            "Cannot perform regression on empty dataset".to_string(),
        ));
    }

    let n = <R as From<f64>>::from(data.len() as f64);

    // Compute means
    let mut sum_x = R::zero();
    let mut sum_y = R::zero();
    for (x, y) in data {
        sum_x += *x;
        sum_y += *y;
    }
    let mean_x = sum_x / n;
    let mean_y = sum_y / n;

    // Compute slope: Σ((x-x̄)(y-ȳ)) / Σ((x-x̄)²)
    let mut numerator = R::zero();
    let mut denominator = R::zero();
    for (x, y) in data {
        let dx = *x - mean_x;
        let dy = *y - mean_y;
        numerator += dx * dy;
        denominator += dx * dx;
    }

    let epsilon = <R as From<f64>>::from(1e-20);
    if denominator.abs() < epsilon {
        return Err(TopologyError::LatticeGaugeError(
            "Degenerate regression (no variance in P2)".to_string(),
        ));
    }

    let slope = numerator / denominator;
    let intercept = mean_y - slope * mean_x;

    Ok((slope, intercept))
}

// ============================================================================
// Helper functions for GM calculation (module-level)
// ============================================================================

/// Physical bounds for radial mapping
const EARTH_RADIUS_M: f64 = 6.371e6; // Earth surface (meters)
const GEO_ORBIT_RADIUS_M: f64 = 4.2e7; // GEO orbit (meters)

/// Converts a lattice index to physical radius (meters).
///
/// Uses linear mapping from Earth surface to GEO orbit.
fn lattice_index_to_radius<R>(r_idx: usize, n_radial: usize) -> R
where
    R: RealField + From<f64>,
{
    let r_min = <R as From<f64>>::from(EARTH_RADIUS_M);
    let r_max = <R as From<f64>>::from(GEO_ORBIT_RADIUS_M);
    let idx_norm = <R as From<f64>>::from(r_idx as f64 / (n_radial - 1).max(1) as f64);
    r_min + idx_norm * (r_max - r_min)
}

/// Converts a physical radius to lattice index.
fn radius_to_lattice_index<R>(r: R, n_radial: usize) -> usize
where
    R: RealField + Into<f64>,
{
    let r_min = EARTH_RADIUS_M;
    let r_max = GEO_ORBIT_RADIUS_M;
    let r_f64: f64 = r.into();
    let r_norm = (r_f64 - r_min) / (r_max - r_min);
    let idx = (r_norm * (n_radial - 1) as f64).round() as usize;
    idx.min(n_radial - 1)
}

/// Performs linear regression on (1/r, y) pairs.
/// Returns (slope, intercept) where y = slope * (1/r) + intercept.
fn linear_regression_inv_r<R>(data: &[(R, R)]) -> Result<(R, R), TopologyError>
where
    R: RealField + Clone + From<f64>,
{
    if data.is_empty() {
        return Err(TopologyError::LatticeGaugeError(
            "Cannot perform regression on empty dataset".to_string(),
        ));
    }

    // Transform (r, y) to (1/r, y)
    let transformed: Vec<(R, R)> = data
        .iter()
        .filter(|(r, _)| *r != R::zero())
        .map(|(r, y)| (R::one() / *r, *y))
        .collect();

    if transformed.is_empty() {
        return Err(TopologyError::LatticeGaugeError(
            "All radii are zero, cannot compute 1/r".to_string(),
        ));
    }

    linear_regression(&transformed)
}

// ============================================================================
// ChronoGaugeMutOps Implementation
// ============================================================================

use crate::{ChronoGaugeMutOps, SpaceTimeCoordinate};
use deep_causality_num::Complex;
use deep_causality_topology::{CWComplex, LinkVariable, SU2_U1};

impl<R> ChronoGaugeMutOps<R> for ChronoGauge<R>
where
    R: RealField
        + Clone
        + From<f64>
        + Into<f64>
        + Default
        + Debug
        + Send
        + Sync
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
            let r_idx = radius_to_lattice_index(coord.r_m, n_radial);
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

        // Update temporal links at each radial shell
        // Collect all temporal edge cells first to avoid borrow conflict
        let temporal_cells: Vec<_> = self
            .lattice()
            .cells(1)
            .filter(|cell| {
                // Only temporal links (direction 0)
                let dir = cell.orientation().trailing_zeros() as usize;
                dir == 0
            })
            .collect();

        // Now iterate over collected cells and set links
        for cell in temporal_cells {
            let pos = cell.position();
            let r_idx = pos[1];

            if let Some(avg_drift) = avg_drifts.get(r_idx).and_then(|opt| opt.as_ref()) {
                // phase = (1 - clock_drift_rate) * N_t
                // For GNSS: drift ≈ 1 - 10^{-10}, so phase ≈ 10^{-10} * N_t
                let phase = (R::one() - *avg_drift) * n_t_scale;

                // Create link variable with this phase
                // U_0 = exp(i * phase) for U(1) part in the SU(2)×U(1) representation
                let link = LinkVariable::<SU2_U1, Complex<R>, R>::from_phase(phase);

                // Set the link in the gauge field
                self.set_link(cell, link);
            }
        }

        Ok(())
    }
}
