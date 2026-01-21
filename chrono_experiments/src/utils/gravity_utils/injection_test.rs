/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use chrono_data_manager::ObservatoryEpoch;
use std::f64::consts::PI;

/// Inject a synthetic Gravitational Wave signal into the clock data.
///
/// Simulates a GW passing through the observatory.
///
/// properties:
/// - Strain h = 1e-9 (Massive signal for validation)
/// - Frequency = 5.4 mHz (To match the detected "hum" region or higher)
/// - Source: Z-axis propagation (+ polarization)
///
/// Physics:
/// A GW stretches/squeezes space transverse to propagation.
/// For + polarization propagating in Z:
/// Metric perturbation h_ab = diag(h(t), -h(t), 0, 0).
///
/// Effect on a clock at position (x, y, z):
/// The "proper time" integral is complex, but for a 1-way link (Sat to Earth),
/// the dominant effect is the "Doppler" shift of the link frequency:
/// y(t) = dnu/nu = 0.5 * n_i * n_j * h_ij(t - L/c)
///
/// Here, we simplify:
/// We assume the "clock bias" is the integral of the rate y(t).
/// Bias perturbation dB(t) = Integral(y(t) dt) ~ (0.5 * h_int / omega) * pattern.
///
/// Pattern for Z-propagation:
/// h_xx = h, h_yy = -h.
/// n = (nx, ny, nz) unit vector from Sat to Earth.
/// n_x = x/r, n_y = y/r.
/// Response F = n_x^2 - n_y^2 = (x^2 - y^2) / r^2 = cos(2*phi) * sin^2(theta).
///
/// So Rate Shift y(t) = 0.5 * h(t) * cos(2phi) * sin^2(theta).
/// Bias Shift dB(t) = Integrate[ 0.5 * h0 * sin(wt) * F ]
///                  = -0.5 * (h0/w) * cos(wt) * F.
///
/// We will inject this dB(t) into `sat.clock_bias_ns`.
///
pub fn inject_gw_signal(epochs: &mut [ObservatoryEpoch], strain_h: f64, freq_hz: f64) {
    println!(">>> INJECTION: Synthesizing GW Signal...");
    println!("    Strain:     {:.2e}", strain_h);
    println!("    Frequency:  {:.2e} Hz", freq_hz);
    println!("    Source:     Z-axis, + Polarization");

    let omega = 2.0 * PI * freq_hz;
    let start_time = epochs.first().map(|e| e.timestamp).unwrap_or(0);

    for epoch in epochs.iter_mut() {
        let t_sec = (epoch.timestamp - start_time) as f64;

        // Signal phase (simple sine wave)
        // h(t) = h0 * sin(omega * t)
        // Integral h(t) = -(h0/omega) * cos(omega * t)
        let gw_integral = -(strain_h / omega) * (omega * t_sec).cos();

        // Assuming bias is in ns.
        // dB (sec) = 0.5 * gw_integral * Pattern
        // dB (ns) = dB (sec) * 1e9

        for (_, sat) in epoch.satellites.iter_mut() {
            let x = sat.x_m;
            let y = sat.y_m;
            let z = sat.z_m;
            let r2 = x * x + y * y + z * z;

            if r2 > 0.0 {
                // Geometric Factor F = (x^2 - y^2) / r^2
                // This is the quadrupole pattern for + pol along Z
                let pattern = (x * x - y * y) / r2;

                // Effect
                let db_ns = 0.5 * gw_integral * pattern * 1e9;

                // Inject
                sat.clock_bias_ns += db_ns;
            }
        }
    }

    println!(
        ">>> INJECTION: Complete. Synthetic signal added to {} epochs.",
        epochs.len()
    );
}
