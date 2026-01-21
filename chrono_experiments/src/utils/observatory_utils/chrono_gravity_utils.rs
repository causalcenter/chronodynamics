/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Core analysis utilities for the observatory.
//! Ported with R: RealField generics.
use chrono_data_manager::{ObservatoryEpoch, SatelliteState};

use crate::utils::observatory_utils::residual_filter::{
    apply_gqcd_cleaning, apply_single_sat_cleaning, check_for_red_noise, detrend_polynomial,
};
use crate::utils::observatory_utils::swarm_analysis::{
    DynamicSwarmResult, analyze_swarm_covariance_dynamic,
};
use deep_causality_num::RealField;
use deep_causality_physics::{EARTH_GM, GALILEO_NOMINAL_RADIUS_M, SPEED_OF_LIGHT};
use std::collections::HashMap;

/// Time delay for pseudo-2T measurement in seconds (2 hours)
pub const TIME_DELAY_S: u64 = 7200;

/// Generic SatelliteState for R calculations
#[derive(Debug, Clone)]
pub struct SatelliteStateR<R> {
    pub sat_id: String,
    pub clock_bias_ns: R,
    pub x_m: R,
    pub y_m: R,
    pub z_m: R,
    pub radius_m: R,
}

#[derive(Debug, Clone)]
pub struct FiveT<R> {
    pub timestamp: i64,
    pub t1_e12: Option<SatelliteStateR<R>>,
    pub t2_e18: Option<SatelliteStateR<R>>,
    pub t3_e14_t0: Option<SatelliteStateR<R>>,
    pub t4_e14_t1h: Option<SatelliteStateR<R>>,
    pub t5_e14_t2h: Option<SatelliteStateR<R>>,
}

#[derive(Debug, Clone)]
pub struct ObservatoryResults<R> {
    pub total_epochs: usize,
    pub valid_5t_epochs: usize,
    pub mean_clock_diff_e14_e12: R,
    pub std_clock_diff_e14_e12: R,
    pub mean_altitude_diff_m: R,
    pub predicted_gr_shift_ns: R,
    pub observed_shift_ns: R,
    pub gr_agreement_pct: R,
    pub clock_rate_altitude_correlation: R,
    pub measured_gr_coefficient: R,
    pub predicted_gr_coefficient: R,
    pub gr_coefficient_agreement_pct: R,
    pub detrended_clock_std_ns: R,
    pub e12_e18_correlation: R,
    pub temporal_autocorr_e14: R,
    pub gw_strain_upper_limit: R,
    pub quadrupole_residual: R,
    pub cross_correlation_coeff: R,
    pub spectral_power_n_hz: R,
    pub geodesic_strain: R,
    pub angular_separation_deg: R,
    pub expected_hd_correlation: R,
    pub hd_significance_pct: R,
    // GQCD Cleaning Pipeline
    pub raw_rms_ns: R,
    pub physics_rms_ns: R,
    pub cleaned_rms_ns: R,
    pub cleaning_ratio: R,
    pub red_noise_coefficient: R,
    pub drift_rate_ns_hr: R,
    // Frequency Domain Analysis (PSD)
    pub peak_12h_power: R,
    pub peak_14h_power: R,
    pub peak_27d_power: R,
    pub nhz_background_power: R,
    // Swarm Analysis
    pub swarm_quadrupole_power: R,
    pub swarm_monopole_power: R,
    pub anomaly_start: String,
    pub anomaly_end: String,
}

pub fn build_5t_fields<R>(epochs: &[ObservatoryEpoch]) -> Vec<FiveT<R>>
where
    R: RealField + Copy + From<f64>,
{
    // Sort epochs by timestamp
    let mut sorted_epochs: Vec<_> = epochs.iter().collect();
    sorted_epochs.sort_by_key(|e| e.timestamp);

    // Build lookup map
    let epoch_map: HashMap<i64, &ObservatoryEpoch> =
        sorted_epochs.iter().map(|e| (e.timestamp, *e)).collect();

    let mut five_t_fields = Vec::new();

    let to_generic = |s: &SatelliteState| -> SatelliteStateR<R> {
        SatelliteStateR {
            sat_id: s.sat_id.clone(),
            clock_bias_ns: R::from(s.clock_bias_ns),
            x_m: R::from(s.x_m),
            y_m: R::from(s.y_m),
            z_m: R::from(s.z_m),
            radius_m: R::from(s.radius_m),
        }
    };

    for epoch in sorted_epochs.iter() {
        let ts = epoch.timestamp;

        let t1_e12 = epoch.satellites.get("E12").map(to_generic);
        let t2_e18 = epoch.satellites.get("E18").map(to_generic);
        let t3_e14_t0 = epoch.satellites.get("E14").map(to_generic);

        let t4_e14_t1h = epoch_map
            .get(&(ts + 3600))
            .and_then(|e| e.satellites.get("E14").map(to_generic));

        let t5_e14_t2h = epoch_map
            .get(&(ts + TIME_DELAY_S as i64))
            .and_then(|e| e.satellites.get("E14").map(to_generic));

        let mut count = 0;
        if t1_e12.is_some() {
            count += 1;
        }
        if t2_e18.is_some() {
            count += 1;
        }
        if t3_e14_t0.is_some() {
            count += 1;
        }
        if t4_e14_t1h.is_some() {
            count += 1;
        }
        if t5_e14_t2h.is_some() {
            count += 1;
        }

        if count >= 3 {
            five_t_fields.push(FiveT {
                timestamp: ts,
                t1_e12,
                t2_e18,
                t3_e14_t0,
                t4_e14_t1h,
                t5_e14_t2h,
            });
        }
    }

    five_t_fields
}

pub fn analyze_5t_fields<R>(fields: &[FiveT<R>]) -> ObservatoryResults<R>
where
    R: RealField
        + Copy
        + From<f64>
        + Into<f64>
        + PartialOrd
        + std::fmt::Display
        + std::fmt::LowerExp,
{
    let mut clock_diffs: Vec<R> = Vec::new();
    let mut altitude_diffs: Vec<R> = Vec::new();
    let mut timestamps: Vec<i64> = Vec::new();

    let mut e12_e18_paired: Vec<(i64, R, R)> = Vec::new();
    let mut e14_t0_t2h_paired: Vec<(R, R)> = Vec::new();
    let mut e18_e14_diffs: Vec<R> = Vec::new();

    for field in fields {
        // Clock differential E14 - E12
        if let (Some(e14), Some(e12)) = (&field.t3_e14_t0, &field.t1_e12) {
            let diff = e14.clock_bias_ns - e12.clock_bias_ns;
            clock_diffs.push(diff);
            altitude_diffs.push(e14.radius_m - e12.radius_m);
            timestamps.push(field.timestamp);
        }

        if let (Some(e12), Some(e18)) = (&field.t1_e12, &field.t2_e18) {
            e12_e18_paired.push((field.timestamp, e12.clock_bias_ns, e18.clock_bias_ns));
        }

        if let (Some(e14_t0), Some(e14_t2h)) = (&field.t3_e14_t0, &field.t5_e14_t2h) {
            e14_t0_t2h_paired.push((e14_t0.clock_bias_ns, e14_t2h.clock_bias_ns));
        }

        if let (Some(e14), Some(e18)) = (&field.t3_e14_t0, &field.t2_e18) {
            e18_e14_diffs.push(e14.clock_bias_ns - e18.clock_bias_ns);
        }
    }

    let n = R::from(clock_diffs.len() as f64);
    let zero = R::zero();

    if clock_diffs.is_empty() {
        return ObservatoryResults {
            total_epochs: fields.len(),
            valid_5t_epochs: 0,
            mean_clock_diff_e14_e12: zero,
            std_clock_diff_e14_e12: zero,
            mean_altitude_diff_m: zero,
            predicted_gr_shift_ns: zero,
            observed_shift_ns: zero,
            gr_agreement_pct: zero,
            clock_rate_altitude_correlation: zero,
            measured_gr_coefficient: zero,
            predicted_gr_coefficient: zero,
            gr_coefficient_agreement_pct: zero,
            detrended_clock_std_ns: zero,
            e12_e18_correlation: zero,
            temporal_autocorr_e14: zero,
            gw_strain_upper_limit: zero,
            quadrupole_residual: zero,
            cross_correlation_coeff: zero,
            spectral_power_n_hz: zero,
            geodesic_strain: zero,
            angular_separation_deg: zero,
            expected_hd_correlation: zero,
            hd_significance_pct: zero,
            raw_rms_ns: zero,
            physics_rms_ns: zero,
            cleaned_rms_ns: zero,
            cleaning_ratio: R::one(),
            red_noise_coefficient: zero,
            drift_rate_ns_hr: zero,
            peak_12h_power: zero,
            peak_14h_power: zero,
            peak_27d_power: zero,
            nhz_background_power: zero,
            swarm_quadrupole_power: zero,
            swarm_monopole_power: zero,
            anomaly_start: "None".into(),
            anomaly_end: "None".into(),
        };
    }

    let sum_diff: R = clock_diffs.iter().fold(zero, |acc, &x| acc + x);
    let mean_clock_diff = sum_diff / n;

    // Detrend
    let detrended_diffs = detrend_polynomial(&timestamps, &clock_diffs);
    let var_sum: R = detrended_diffs.iter().fold(zero, |acc, &x| acc + x * x);
    let std_clock_diff = (var_sum / (n - R::one())).sqrt();

    let sum_alt: R = altitude_diffs.iter().fold(zero, |acc, &x| acc + x);
    let mean_altitude_diff = sum_alt / n;

    // GR Shift - Use actual measured radii, not assumed values
    // E12 is circular at nominal Galileo radius
    let r_e12 = R::from(GALILEO_NOMINAL_RADIUS_M);
    let r_e14 = r_e12 + mean_altitude_diff;
    let c = R::from(SPEED_OF_LIGHT);
    let gm = R::from(EARTH_GM);

    let gr_factor = gm / (c * c);
    let predicted_shift_per_s = gr_factor * (R::one() / r_e14 - R::one() / r_e12);
    let predicted_gr_shift_ns = predicted_shift_per_s * R::from(1e9);

    let obs_span = n * R::from(30.0);
    let observed_shift_ns = mean_clock_diff / obs_span;

    let r_e14_avg = r_e12 + mean_altitude_diff;
    let r_e12_avg = r_e12;
    let cleaned = apply_gqcd_cleaning(&timestamps, &clock_diffs, r_e14_avg, r_e12_avg, None, None);

    let raw_rms_ns = cleaned.raw_rms_ns;
    let physics_rms_ns = cleaned.physics_rms_ns;
    let cleaned_rms_ns = cleaned.gw_rms_ns;

    let red_noise_coeff = check_for_red_noise(&cleaned.gw_candidates);
    let (_, b_drift, _) = cleaned.drift_coefficients;

    // --- New Calculations ---
    let gr_agreement_pct = if predicted_gr_shift_ns.abs() > R::from(1e-12) {
        let diff = (observed_shift_ns - predicted_gr_shift_ns).abs();
        if diff < predicted_gr_shift_ns {
            (R::one() - diff / predicted_gr_shift_ns) * R::from(100.0)
        } else {
            zero
        }
    } else {
        zero
    };

    let v1: Vec<R> = e12_e18_paired.iter().map(|t| t.1).collect();
    let v2: Vec<R> = e12_e18_paired.iter().map(|t| t.2).collect();
    let e12_e18_corr = compute_correlation(&v1, &v2);

    // Calculate clock rates from DETRENDED diffs (not raw diffs)
    let mut clock_rates_ns_s: Vec<R> = Vec::new();
    let mut rate_altitudes: Vec<R> = Vec::new();
    for i in 1..detrended_diffs.len() {
        let dt_s = R::from((timestamps[i] - timestamps[i - 1]) as f64);
        if dt_s > R::zero() && dt_s < R::from(3600.0) {
            let d_clock = detrended_diffs[i] - detrended_diffs[i - 1];
            let clock_rate = d_clock / dt_s;
            clock_rates_ns_s.push(clock_rate);
            rate_altitudes.push((altitude_diffs[i] + altitude_diffs[i - 1]) / R::from(2.0));
        }
    }

    let n_rates = clock_rates_ns_s.len();
    let clock_alt_corr = if n_rates > 10 {
        compute_correlation(&clock_rates_ns_s, &rate_altitudes)
    } else {
        compute_correlation(&clock_diffs, &altitude_diffs)
    };

    // GR Coefficient Calculation: Linear regression of rate vs altitude
    let (measured_gr_coeff, detrend_std) = if n_rates > 10 {
        let mean_rate =
            clock_rates_ns_s.iter().fold(zero, |acc, &x| acc + x) / R::from(n_rates as f64);
        let mean_alt =
            rate_altitudes.iter().fold(zero, |acc, &x| acc + x) / R::from(n_rates as f64);

        let mut numerator = zero;
        let mut denominator = zero;
        for i in 0..n_rates {
            let alt_dev = rate_altitudes[i] - mean_alt;
            let rate_dev = clock_rates_ns_s[i] - mean_rate;
            numerator += alt_dev * rate_dev;
            denominator += alt_dev * alt_dev;
        }

        let slope_ns_s_per_m = if denominator.abs() > R::from(1e-20) {
            numerator / denominator
        } else {
            zero
        };
        let slope_ns_s_per_km = slope_ns_s_per_m * R::from(1000.0);

        // Compute detrended residuals
        let detrended_rates: Vec<R> = clock_rates_ns_s
            .iter()
            .zip(rate_altitudes.iter())
            .map(|(&rate, &alt)| rate - slope_ns_s_per_m * (alt - mean_alt) - mean_rate)
            .collect();
        let detrended_var = detrended_rates.iter().fold(zero, |acc, &x| acc + x * x)
            / R::from((n_rates - 1).max(1) as f64);
        (slope_ns_s_per_km, detrended_var.sqrt())
    } else {
        (zero, std_clock_diff)
    };

    // Predicted GR coefficient (from theory)
    let r_galileo = R::from(GALILEO_NOMINAL_RADIUS_M);
    let predicted_gr_coeff =
        (gm / (c * c)) / (r_galileo * r_galileo) * R::from(1e9) * R::from(1000.0);

    // Agreement calculation (log-ratio based)
    let gr_coeff_agree =
        if predicted_gr_coeff.abs() > R::from(1e-15) && measured_gr_coeff.abs() > R::from(1e-15) {
            let sign_match = (measured_gr_coeff > zero) == (predicted_gr_coeff > zero);
            if sign_match {
                let ratio: f64 = (measured_gr_coeff.abs() / predicted_gr_coeff.abs()).into();
                let log_ratio = ratio.log10().abs();
                if log_ratio < 1.0 {
                    R::from((1.0 - log_ratio) * 100.0)
                } else {
                    zero
                }
            } else {
                zero
            }
        } else {
            zero
        };

    let corr_auto = if clock_diffs.len() > 1 {
        let c1 = clock_diffs[0..clock_diffs.len() - 1].to_vec();
        let c2 = clock_diffs[1..clock_diffs.len()].to_vec();
        compute_correlation(&c1, &c2)
    } else {
        zero
    };
    // ------------------------

    ObservatoryResults {
        total_epochs: fields.len(),
        valid_5t_epochs: clock_diffs.len(),
        mean_clock_diff_e14_e12: mean_clock_diff,
        std_clock_diff_e14_e12: std_clock_diff,
        mean_altitude_diff_m: mean_altitude_diff,
        predicted_gr_shift_ns,
        observed_shift_ns,
        gr_agreement_pct,
        clock_rate_altitude_correlation: clock_alt_corr,
        measured_gr_coefficient: measured_gr_coeff,
        predicted_gr_coefficient: predicted_gr_coeff,
        gr_coefficient_agreement_pct: gr_coeff_agree,
        detrended_clock_std_ns: detrend_std,
        e12_e18_correlation: e12_e18_corr,
        temporal_autocorr_e14: corr_auto,
        gw_strain_upper_limit: zero,
        quadrupole_residual: zero,
        cross_correlation_coeff: zero,
        spectral_power_n_hz: zero,
        geodesic_strain: zero,
        angular_separation_deg: zero,
        expected_hd_correlation: zero,
        hd_significance_pct: zero,
        raw_rms_ns,
        physics_rms_ns,
        cleaned_rms_ns,
        cleaning_ratio: if cleaned_rms_ns > zero {
            raw_rms_ns / cleaned_rms_ns
        } else {
            R::one()
        },
        red_noise_coefficient: red_noise_coeff,
        drift_rate_ns_hr: b_drift,
        peak_12h_power: zero,
        peak_14h_power: zero,
        peak_27d_power: zero,
        nhz_background_power: zero,
        swarm_quadrupole_power: zero,
        swarm_monopole_power: zero,
        anomaly_start: "None".into(),
        anomaly_end: "None".into(),
    }
}

pub fn analyze_swarm_vectors<R>(epochs: &[ObservatoryEpoch]) -> DynamicSwarmResult<R>
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    // Use all 7 satellites for Swarm Analysis (matching "7-Sat Fleet Health")
    let sats = ["E11", "E12", "E14", "E18", "E19", "E08", "E09"];
    let mut aligned_timestamps: Vec<i64> = Vec::new();
    // Dynamic sizing based on satellite list
    let mut sat_data_biases: Vec<Vec<R>> = vec![Vec::new(); sats.len()];
    let mut sat_data_radii: Vec<Vec<R>> = vec![Vec::new(); sats.len()];

    let mut sorted: Vec<&ObservatoryEpoch> = epochs.iter().collect();
    sorted.sort_by_key(|e| e.timestamp);

    for epoch in sorted {
        let mut present = true;
        for id in &sats {
            if !epoch.satellites.contains_key(*id) {
                present = false;
                break;
            }
        }

        if present {
            aligned_timestamps.push(epoch.timestamp);
            for (i, id) in sats.iter().enumerate() {
                let s = &epoch.satellites[*id];
                sat_data_biases[i].push(R::from(s.clock_bias_ns));
                sat_data_radii[i].push(R::from(s.radius_m));
            }
        }
    }

    if aligned_timestamps.len() < 10 {
        return DynamicSwarmResult {
            eigenvalues: vec![],
            eigenvectors: vec![],
            quadrupole_power: R::zero(),
            monopole_power: R::zero(),
            num_satellites: 0,
        };
    }

    let mut cleaned_vectors: Vec<Vec<R>> = Vec::new();
    let r_nominal = R::from(GALILEO_NOMINAL_RADIUS_M);

    for i in 0..sats.len() {
        let biases = &sat_data_biases[i];
        let radii = &sat_data_radii[i];

        let cleaned = apply_single_sat_cleaning(&aligned_timestamps, biases, radii, r_nominal);
        cleaned_vectors.push(cleaned.gw_candidates);
    }

    let n = cleaned_vectors[0].len();
    let mut residual_rows = Vec::with_capacity(n);

    for t in 0..n {
        let mut row_vec = Vec::with_capacity(sats.len());
        // Dynamic iteration
        for (_s, input_vec) in cleaned_vectors.iter().enumerate().take(sats.len()) {
            row_vec.push(input_vec[t]);
        }
        residual_rows.push(row_vec);
    }

    analyze_swarm_covariance_dynamic(&residual_rows)
}

pub fn compute_correlation<R>(x: &[R], y: &[R]) -> R
where
    R: RealField + Copy + From<f64> + PartialOrd,
{
    let n = x.len().min(y.len());
    if n == 0 {
        return R::zero();
    }

    let n_gen = R::from(n as f64);
    let zero = R::zero();

    let sum_x = x.iter().take(n).fold(zero, |acc, &v| acc + v);
    let sum_y = y.iter().take(n).fold(zero, |acc, &v| acc + v);
    let mean_x = sum_x / n_gen;
    let mean_y = sum_y / n_gen;

    let mut cov = zero;
    let mut var_x = zero;
    let mut var_y = zero;

    for i in 0..n {
        let dx = x[i] - mean_x;
        let dy = y[i] - mean_y;
        cov += dx * dy;
        var_x += dx * dx;
        var_y += dy * dy;
    }

    if var_x > zero && var_y > zero {
        cov / (var_x.sqrt() * var_y.sqrt())
    } else {
        zero
    }
}

impl<R> ObservatoryResults<R>
where
    R: RealField
        + Copy
        + From<f64>
        + Into<f64>
        + PartialOrd
        + std::fmt::Display
        + std::fmt::LowerExp,
{
    pub fn print(&self, year: &str) {
        let mean_alt: f64 = self.mean_altitude_diff_m.into();
        let gr_agree: f64 = self.gr_agreement_pct.into();
        let corr_alt: f64 = self.clock_rate_altitude_correlation.into();
        let gr_coeff_agree: f64 = self.gr_coefficient_agreement_pct.into();
        let corr_e12_e18: f64 = self.e12_e18_correlation.into();
        let corr_auto: f64 = self.temporal_autocorr_e14.into();

        println!("\n╔══════════════════════════════════════════════════════════════════════════╗");
        println!(
            "║  OBSERVATORY RESULTS - {}                                                ║",
            year
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  DATA QUALITY                                                            ║");
        println!(
            "║    Total epochs analyzed:         {:<15}                        ║",
            self.total_epochs
        );
        println!(
            "║    Valid 5T measurements:         {:<15}                        ║",
            self.valid_5t_epochs
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  CLOCK DIFFERENTIAL (E14 - E12)                                          ║");
        println!(
            "║    Mean clock difference:      {:<+10.4e} ns                             ║",
            self.mean_clock_diff_e14_e12
        );
        println!(
            "║    Std deviation:              {:<10.4e} ns                              ║",
            self.std_clock_diff_e14_e12
        );
        println!(
            "║    Mean altitude difference:   {:<10.2} m                                ║",
            mean_alt
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  GENERAL RELATIVITY VERIFICATION                                         ║");
        println!(
            "║    Predicted GR shift:         {:<10.4e} ns/s                            ║",
            self.predicted_gr_shift_ns
        );
        println!(
            "║    Observed shift:             {:<10.4e} ns/s                            ║",
            self.observed_shift_ns
        );
        println!(
            "║    GR agreement:               {:<.2}%                                     ║",
            gr_agree
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  ISOLATED GR ANALYSIS (Systematics Removed)                              ║");
        println!(
            "║    Clock-Altitude correlation: {:<+7.4}                                      ║",
            corr_alt
        );
        println!(
            "║    Measured GR coeff:          {:<10.4e} ns/km                            ║",
            self.measured_gr_coefficient
        );
        println!(
            "║    Predicted GR coeff:         {:<10.4e} ns/km                            ║",
            self.predicted_gr_coefficient
        );
        println!(
            "║    GR coefficient agreement:   {:<.2}%                                     ║",
            gr_coeff_agree
        );
        println!(
            "║    Detrended residual:         {:<10.4e} ns (GR removed)                 ║",
            self.detrended_clock_std_ns
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  CORRELATION ANALYSIS                                                    ║");
        println!(
            "║    E12-E18 correlation:        {:<+7.4}                                      ║",
            corr_e12_e18
        );
        println!(
            "║    E14 temporal autocorr:      {:<+7.4}                                      ║",
            corr_auto
        );
        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  SWARM COVARIANCE ANALYSIS (7-Sat Fleet Health)                          ║");
        println!(
            "║    Differential Power:         {:<10.4e}                                 ║",
            self.swarm_quadrupole_power
        );
        println!(
            "║    Common Mode Power:          {:<10.4e}                                 ║",
            self.swarm_monopole_power
        );

        let mono_abs: f64 = self.swarm_monopole_power.into();
        let quad_abs: f64 = self.swarm_quadrupole_power.into();

        if mono_abs.abs() > quad_abs.abs() {
            println!("║    ✓ Common mode dominant - Fleet operating coherently                  ║");
        } else {
            println!("║    ⚠ DIFFERENTIAL DOMINANT - Fleet coherence issue detected             ║");
        }

        println!("╠══════════════════════════════════════════════════════════════════════════╣");
        println!("║  INTERPRETATION                                                          ║");
        println!("├──────────────────────────────────────────────────────────────────────────┤");

        if corr_e12_e18 > 0.5 {
            println!("║  ~ Moderate E12-E18 correlation - partial common-mode rejection         ║");
        } else if corr_e12_e18 < -0.8 {
            println!("║  ! Low E12-E18 correlation - satellites experiencing different fields   ║");
        } else {
            println!("║  - Neutral E12-E18 correlation                                          ║");
        }

        if corr_auto > 0.9 {
            println!("║  ✓ High E14 temporal autocorr indicates stable orbital dynamics         ║");
        } else {
            println!("║  ⚠ Low E14 temporal autocorr indicates orbital instability              ║");
        }
        println!("╚══════════════════════════════════════════════════════════════════════════╝");
    }
}
