//! Observable computation utilities for chrono-gauge experiments.
//!
//! These formulas match the legacy gqcd_chrono_manifold DEC calculations.

use crate::SpaceTimeCoordinate;
use crate::types::gauge_types::EpochMetrics;
use crate::utils::proces_utils::apply_mad_filter;
use deep_causality_num::{Float, FromPrimitive, RealField, ToPrimitive};
use deep_causality_physics::{EARTH_GM, EARTH_RADIUS, EARTH_ROTATION_RATE, SPEED_OF_LIGHT};

pub const C_SQ: f64 = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

/// Computes time-velocity correlation across ALL coordinates at dataset level.
///
/// This is the correct way to compute the correlation - accumulating E14/E18
/// data across all epochs before correlation (matching legacy approach).
///
/// Reference value: +0.9991
pub fn compute_dataset_time_velocity_correlation<R>(all_coordinates: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    compute_time_velocity_correlation(all_coordinates)
}

/// Computes all epoch metrics DIRECTLY from satellite data.
pub fn compute_epoch_metrics_fast<R>(satellites: &[SpaceTimeCoordinate<R>]) -> EpochMetrics<R>
where
    R: RealField + Clone + From<f64> + Into<f64> + Default + Float + FromPrimitive + ToPrimitive,
{
    // ----------------------------------
    // 1. Laplacian (mass density) - 0 in vacuum
    // ----------------------------------
    let mass_density = compute_laplacian(satellites);

    // ----------------------------------
    // 2. Curl magnitude - exactly 0 for scalar fields
    // ----------------------------------
    let curl_magnitude = R::zero();

    // ----------------------------------
    // 3. Gradient magnitude from clock drift rate
    // ----------------------------------
    let gradient_magnitude = compute_gradient_magnitude(satellites);

    // ----------------------------------
    // 4. Vorticity z-component - 0 for irrotational field
    // ----------------------------------
    let vorticity_z = R::zero();

    // ----------------------------------
    // 5. Time-velocity correlation (Momentum)
    // ----------------------------------
    let time_velocity_corr = compute_time_velocity_correlation(satellites);

    // ----------------------------------
    // 6. Virial ratio (Energy)
    // ----------------------------------
    // ----------------------------------
    // 6. Virial ratio (Energy)
    // ----------------------------------
    let virial_ratio = compute_virial_ratio(satellites);

    // ----------------------------------
    // 8. Derived Earth Rotation (from Momentum optimization)
    // ----------------------------------
    let derived_rotation = estimate_earth_rotation(satellites);

    // ----------------------------------
    // 7. Tolman temperature consistency
    // ----------------------------------
    let tolman_temp = compute_tolman_consistency(satellites);

    // Fleet coherence from clock rate statistics
    let (coherence_mean, coherence_var) = compute_fleet_coherence(satellites);

    EpochMetrics {
        mass_density,
        curl_magnitude,
        gradient_magnitude,
        vorticity_z,
        time_velocity_corr,
        virial_ratio,
        tolman_temp,
        coherence_mean,
        coherence_var,
        derived_rotation,
        derived_gm: R::default(),
        derived_mass: R::default(),
        derived_gravity: R::default(),
    }
}

/// Computes discrete Laplacian approximation of clock drift field.
fn compute_laplacian<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.len() < 4 {
        return R::zero();
    }

    // Use finite differences: Δf ≈ Σ(f_neighbor - f_center) / h²
    let mut laplacian_sum = 0.0f64;
    let n = satellites.len();

    for i in 0..n {
        let rate_i: f64 = satellites[i].clock_drift_rate.into();
        let pos_i = [
            satellites[i].position[0].into(),
            satellites[i].position[1].into(),
            satellites[i].position[2].into(),
        ];

        let mut neighbor_contrib = 0.0f64;
        let mut neighbor_count = 0;

        for (j, sat_j) in satellites.iter().enumerate().take(n) {
            if i == j {
                continue;
            }

            let pos_j: [f64; 3] = [
                sat_j.position[0].into(),
                sat_j.position[1].into(),
                sat_j.position[2].into(),
            ];

            let dx = pos_j[0] - pos_i[0];
            let dy = pos_j[1] - pos_i[1];
            let dz = pos_j[2] - pos_i[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;

            if dist_sq > 1e4 {
                let rate_j: f64 = sat_j.clock_drift_rate.into();
                neighbor_contrib += (rate_j - rate_i) / dist_sq;
                neighbor_count += 1;
            }
        }

        if neighbor_count > 0 {
            laplacian_sum += neighbor_contrib / neighbor_count as f64;
        }
    }

    <R as From<f64>>::from(laplacian_sum / n as f64)
}

/// Computes gradient magnitude approximating DEC exterior derivative.
///
/// Legacy: Uses manifold data = 1.0 + clock_drift_rate
/// Gradient is the L2 norm of differences (not divided by distance)
/// Reference: ~1e-10 to ~1e-14
fn compute_gradient_magnitude<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.len() < 2 {
        return R::zero();
    }

    // Legacy manifold data: chrono_value = 1.0 + clock_drift_rate
    let chrono_values: Vec<f64> = satellites
        .iter()
        .map(|s| 1.0 + s.clock_drift_rate.into())
        .collect();

    // Compute gradient as pairwise differences (simulating DEC exterior derivative)
    // Legacy builds tetrahedron from 4 satellites, gradient is on edges
    let mut grad_components: Vec<f64> = Vec::new();

    for i in 0..chrono_values.len() {
        for j in (i + 1)..chrono_values.len() {
            // Edge gradient = difference between vertex values
            let d_chrono = chrono_values[j] - chrono_values[i];
            grad_components.push(d_chrono);
        }
    }

    if grad_components.is_empty() {
        return R::zero();
    }

    // L2 norm of all gradient components (matches legacy grad_mag calculation)
    let grad_sq_sum: f64 = grad_components.iter().map(|g| g * g).sum();
    <R as From<f64>>::from(grad_sq_sum.sqrt())
}

/// Computes time-velocity correlation matching legacy formula.
///
/// Legacy: Uses E14/E18 satellites, correlates kinetic_rate vs v²_eci
/// Reference value: +0.999
fn compute_time_velocity_correlation<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.len() < 2 {
        return R::one();
    }

    let mut kinetic_rates: Vec<f64> = Vec::new();
    let mut eci_velocity_sq: Vec<f64> = Vec::new();

    for sat in satellites {
        // CRITICAL: Only E14 and E18 have elliptical orbits (launch failure).
        // They have significant velocity variation enabling strong correlation.
        // Regular circular-orbit satellites have constant velocity = weak correlation.
        if sat.sat_id != 14 && sat.sat_id != 18 {
            continue;
        }

        // Compute r from position (matching legacy exactly)
        let pos_x: f64 = sat.position[0].into();
        let pos_y: f64 = sat.position[1].into();
        let pos_z: f64 = sat.position[2].into();
        let r = (pos_x * pos_x + pos_y * pos_y + pos_z * pos_z).sqrt();

        if r < 1e6 {
            continue; // Skip invalid positions
        }

        // Gravitational time dilation: Δτ/τ = -GM/(r·c²)
        let grav_dilation = -EARTH_GM / (r * C_SQ);

        // Kinetic rate = total_rate - gravitational_contribution
        let total_rate: f64 = sat.clock_drift_rate.into();
        let rate_kinetic = total_rate - grav_dilation;

        // Convert ECEF velocity to ECI (matching legacy exactly)
        let vel_x: f64 = sat.velocity[0].into();
        let vel_y: f64 = sat.velocity[1].into();
        let vel_z: f64 = sat.velocity[2].into();
        let vx_eci = vel_x + (-EARTH_ROTATION_RATE * pos_y);
        let vy_eci = vel_y + (EARTH_ROTATION_RATE * pos_x);
        let vz_eci = vel_z;
        let v_sq_eci = vx_eci * vx_eci + vy_eci * vy_eci + vz_eci * vz_eci;

        kinetic_rates.push(rate_kinetic);
        eci_velocity_sq.push(v_sq_eci);
    }

    if kinetic_rates.is_empty() {
        return R::one();
    }

    // Apply MAD Filter to remove clock jumps/outliers (Legacy robustness matching)
    // Sigma = 3.5 captures 99.9% of Gaussian data, rejecting only extreme outliers
    let rates_f64: Vec<f64> = kinetic_rates.clone();
    let rates_clean = apply_mad_filter(&rates_f64, 3.5);

    // Filter index-matched velocity
    let mut clean_rates = Vec::new();
    let mut clean_vel = Vec::new();

    // Create a set of "good" rates
    let valid_rate_set: std::collections::HashSet<String> = rates_clean
        .iter()
        .map(|r| format!("{:.20}", r)) // Use string representation for safe float comparison
        .collect();

    for (i, rate) in kinetic_rates.iter().enumerate() {
        if valid_rate_set.contains(&format!("{:.20}", rate)) {
            clean_rates.push(*rate);
            clean_vel.push(eci_velocity_sq[i]);
        }
    }

    if clean_rates.len() < 2 {
        return R::one();
    }

    // Pearson correlation on CLEANED data
    let n = clean_rates.len() as f64;
    let mean_rate: f64 = clean_rates.iter().sum::<f64>() / n;
    let mean_vel: f64 = clean_vel.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_rate = 0.0;
    let mut var_vel = 0.0;

    for i in 0..clean_rates.len() {
        let d_rate = clean_rates[i] - mean_rate;
        let d_vel = clean_vel[i] - mean_vel;
        cov += d_rate * d_vel;
        var_rate += d_rate * d_rate;
        var_vel += d_vel * d_vel;
    }

    let denom = (var_rate * var_vel).sqrt();
    if denom > 1e-30 {
        <R as From<f64>>::from(cov / denom)
    } else {
        R::one()
    }
}

/// Computes virial ratio matching legacy formula.
///
/// Legacy: mean_kinetic / (GM * mean_potential) where potential = 1/r
/// Reference: ~0.5 for bound orbits
fn compute_virial_ratio<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.is_empty() {
        return R::zero();
    }

    let mut term_kinetic_sum = 0.0f64;
    let mut term_potential_sum = 0.0f64;
    let mut count = 0usize;

    for sat in satellites {
        let r: f64 = sat.r_m.into();
        if r < 1e6 {
            continue;
        }

        // Term_Kinetic = 0.5 × v²
        let vel_x: f64 = sat.velocity[0].into();
        let vel_y: f64 = sat.velocity[1].into();
        let vel_z: f64 = sat.velocity[2].into();
        let v_sq = vel_x * vel_x + vel_y * vel_y + vel_z * vel_z;
        let term_kinetic = 0.5 * v_sq;

        // Term_Potential = 1/r
        let term_potential = 1.0 / r;

        term_kinetic_sum += term_kinetic;
        term_potential_sum += term_potential;
        count += 1;
    }

    if count > 0 {
        let mean_kinetic = term_kinetic_sum / count as f64;
        let mean_potential = term_potential_sum / count as f64;
        let potential_energy = EARTH_GM * mean_potential;

        if potential_energy.abs() > 1e-10 {
            <R as From<f64>>::from(mean_kinetic.abs() / potential_energy.abs())
        } else {
            R::zero()
        }
    } else {
        R::zero()
    }
}

/// Computes Tolman consistency matching legacy formula.
///
/// Legacy: 1 - (std_dev/mean).abs() of tau_ratios
/// where tau_ratio = sqrt(dilation_sat / dilation_ground)
/// Reference: 1.00
fn compute_tolman_consistency<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.is_empty() {
        return R::one();
    }

    // Collect tau ratios
    let mut tau_ratios: Vec<f64> = Vec::new();

    for sat in satellites {
        let r: f64 = sat.r_m.into();
        if r < 1e6 {
            continue;
        }

        // Legacy formula: tau_ratio = sqrt(dilation_sat / dilation_ground)
        // where dilation = 1 - GM/(c² × r)
        let dilation_ground = 1.0 - EARTH_GM / (C_SQ * EARTH_RADIUS);
        let dilation_sat = 1.0 - EARTH_GM / (C_SQ * r);

        // Avoid division by zero
        if dilation_ground.abs() > 1e-20 {
            let tau_ratio = (dilation_sat / dilation_ground).sqrt();
            tau_ratios.push(tau_ratio);
        }
    }

    if tau_ratios.len() < 2 {
        return R::one();
    }

    // Calculate mean and std dev
    let n = tau_ratios.len() as f64;
    let mean: f64 = tau_ratios.iter().sum::<f64>() / n;
    let variance: f64 = tau_ratios.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();

    // Legacy: 1 - (std_dev/mean).abs()
    if mean.abs() > 1e-20 {
        <R as From<f64>>::from(1.0 - (std_dev / mean).abs())
    } else {
        R::one()
    }
}

/// Computes fleet coherence (mean, variance) of clock drift rates.
fn compute_fleet_coherence<R>(satellites: &[SpaceTimeCoordinate<R>]) -> (R, R)
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.is_empty() {
        return (R::zero(), R::zero());
    }

    let rates: Vec<f64> = satellites
        .iter()
        .map(|s| s.clock_drift_rate.into())
        .collect();

    let n = rates.len() as f64;
    let mean: f64 = rates.iter().sum::<f64>() / n;
    let variance: f64 = rates.iter().map(|r| (r - mean).powi(2)).sum::<f64>() / n;

    (
        <R as From<f64>>::from(mean),
        <R as From<f64>>::from(variance),
    )
}

/// Estimates Earth's rotation rate by maximizing time-velocity correlation.
///
/// Iterates over candidate Omega values around the known Earth rotation rate.
/// The Omega that maximizes the correlation is the "derived" rotation.
///
/// Returns the derived rate (rad/s).
pub fn estimate_earth_rotation<R>(satellites: &[SpaceTimeCoordinate<R>]) -> R
where
    R: RealField + Clone + From<f64> + Into<f64> + Default,
{
    if satellites.len() < 2 {
        return R::zero();
    }

    // Filter for E14/E18
    let valid_indices: Vec<usize> = satellites
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            let id = s.sat_id;
            let r: f64 = s.r_m.into();
            (id == 14 || id == 18) && r > 1e6
        })
        .map(|(i, _)| i)
        .collect();

    if valid_indices.len() < 10 {
        return R::zero(); // Insufficient data for estimation
    }

    // Pre-calculate invariant terms
    let rates_kinetic: Vec<f64> = valid_indices
        .iter()
        .map(|&i| {
            let sat = &satellites[i];
            let r: f64 = sat.r_m.into();
            let grav_d = -EARTH_GM / (r * C_SQ);
            let total: f64 = sat.clock_drift_rate.into();
            total - grav_d
        })
        .collect();

    let ecef_vels: Vec<([f64; 3], [f64; 3])> = valid_indices
        .iter()
        .map(|&i| {
            let sat = &satellites[i];
            let p = &sat.position;
            let v = &sat.velocity;
            (
                [p[0].into(), p[1].into(), p[2].into()],
                [v[0].into(), v[1].into(), v[2].into()],
            )
        })
        .collect();

    // Clean rates for outliers (using same MAD logic as main correlation)
    // We do this ONCE on the kinetic rates assuming outliers are rate-based
    let rates_clean = apply_mad_filter(&rates_kinetic, 3.5);
    let valid_rate_set: std::collections::HashSet<String> =
        rates_clean.iter().map(|r| format!("{:.20}", r)).collect();

    // Filter indices again based on clean rates
    let final_indices: Vec<usize> = rates_kinetic
        .iter()
        .enumerate()
        .filter(|(_, r)| valid_rate_set.contains(&format!("{:.20}", r)))
        .map(|(i, _)| i)
        .collect();

    if final_indices.len() < 5 {
        return R::zero();
    }

    let clean_rates: Vec<f64> = final_indices.iter().map(|&i| rates_kinetic[i]).collect();
    let clean_data_refs: Vec<&([f64; 3], [f64; 3])> =
        final_indices.iter().map(|&i| &ecef_vels[i]).collect();

    // Optimization: Scan Omega from 0.5x to 1.5x Earth Rate
    // Coarse scan then fine tune? No, simple grid scan is enough for demonstration.
    // Range: 0.0 to 1.5e-4. (Earth is ~7.29e-5).
    // Let's use 20 steps.

    let mut best_omega = 0.0;
    let mut max_corr = -1.0;

    let steps = 20;
    let start_omega = 0.0;
    let end_omega = 1.5e-4;
    let step_size = (end_omega - start_omega) / steps as f64;

    for i in 0..=steps {
        let omega = start_omega + i as f64 * step_size;

        // Calculate v^2 list for this omega
        let eci_sq: Vec<f64> = clean_data_refs
            .iter()
            .map(|(p, v)| {
                let vx = v[0] + (-omega * p[1]);
                let vy = v[1] + (omega * p[0]);
                let vz = v[2];
                vx * vx + vy * vy + vz * vz
            })
            .collect();

        // Calc correlation
        let corr = calculate_pearson_f64(&clean_rates, &eci_sq);

        if corr > max_corr {
            max_corr = corr;
            best_omega = omega;
        }
    }

    // Fine tune around best (zoom in 10x)
    let fine_start = (best_omega - step_size).max(0.0);
    let fine_end = best_omega + step_size;
    let fine_step = (fine_end - fine_start) / 20.0;

    for i in 0..=20 {
        let omega = fine_start + i as f64 * fine_step;
        let eci_sq: Vec<f64> = clean_data_refs
            .iter()
            .map(|(p, v)| {
                let vx = v[0] + (-omega * p[1]);
                let vy = v[1] + (omega * p[0]);
                let vz = v[2];
                vx * vx + vy * vy + vz * vz
            })
            .collect();

        let corr = calculate_pearson_f64(&clean_rates, &eci_sq);
        if corr > max_corr {
            max_corr = corr;
            best_omega = omega;
        }
    }

    <R as From<f64>>::from(best_omega)
}

fn calculate_pearson_f64(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;

    for (xi, yi) in x.iter().zip(y.iter()) {
        let dx = xi - mx;
        let dy = yi - my;
        cov += dx * dy;
        vx += dx * dx;
        vy += dy * dy;
    }

    if vx * vy == 0.0 {
        0.0
    } else {
        cov / (vx * vy).sqrt()
    }
}
