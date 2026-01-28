use crate::theories::chrono_dynamics::chrono_ops::chrono_ops_gauge::ChronoOpsGauge;
use crate::theories::chrono_dynamics::chrono_ops::chrono_ops_gauge_mut::ChronoOpsGaugeMut;
use crate::{
    ChronoGauge, EARTH_GM, EARTH_RADIUS_EQUATORIAL, GEO_ORBIT_RADIUS_M, SPEED_OF_LIGHT,
    SpaceTimeCoordinate,
};
use deep_causality_topology::Lattice;
use std::sync::Arc;

#[test]
fn test_solve_gm_smooth_synthetic() {
    // 1. Setup ideal parameters
    // Use Real Earth GM for physical consistency with calibration
    let gm = 3.986004418e14; // Earth GM (m^3/s^2)
    let c = SPEED_OF_LIGHT;
    let _c_sq = c * c;

    // 2. Generate synthetic data for circular orbits (Virial assumption)
    // Range from Earth surface to GEO.
    // Generate 100 points.
    let _n_points = 100;
    let _r_min = EARTH_RADIUS_EQUATORIAL;
    let _r_max = GEO_ORBIT_RADIUS_M;
    // 2. Generate Synthetic Data (Schwarzschild Metric)
    let c = 299_792_458.0;
    let gm = 3.986004418e14 * 100000.0; // Scaled up for numerical stability
    let r_min = 6_378_137.0;
    let r_max = 42_164_000.0;

    let mut data = Vec::new();
    let n_points = 1000; // Increased density for better integral

    for i in 0..n_points {
        let _t = i as f64;
        let r = r_min + (r_max - r_min) * (i as f64 / n_points as f64);

        // Schwarzschild Metric for Circular Orbit
        // dtau = sqrt(1 - 3GM/rc^2) dt
        // 1 - 2GM/rc^2 (Gravitational) - v^2/c^2 (Kinetic)
        // v^2 = GM/r => Total = 1 - 3GM/rc^2
        let term = 1.0 - 3.0 * gm / (r * c * c);
        let dtau_dt = term.sqrt();
        let drift = dtau_dt - 1.0;

        // Circular Velocity: v = (GM/r)
        let v = (gm / r).sqrt();

        data.push(SpaceTimeCoordinate {
            timestamp: i as u64,
            sat_id: 0,
            r_m: r,
            v_ms: v,
            clock_bias_s: 0.0,
            position: [r, 0.0, 0.0],
            velocity: [0.0, v, 0.0],
            clock_drift_rate: drift,
        });
    }

    // 3. Create Gauge Field
    // 32^4 resolution
    let lattice = Arc::new(Lattice::new([32, 32, 32, 32], [true, true, true, true]));
    // Use smoothed population which encodes velocity
    let mut field = ChronoGauge::<f64>::identity(lattice, 1.0).with_source(data);
    field.populate_smooth_links_from_source().unwrap();

    // 4. Solve for GM
    let derived_gm = field.solve_gm().expect("Failed to solve GM");

    println!("Input GM:   {}", gm);
    println!("Derived GM: {}", derived_gm);

    let error = (derived_gm - gm).abs() / gm;
    println!("Error:      {:.2}%", error * 100.0);

    // We used Geometric Wilson Action with Dynamical Virial Factor.
    // Error comes from finite difference discretization (~4%) and virial estimation noise.
    // Expect < 5% error.
    assert!(
        error < 0.05,
        "GM error too high (expected < 5%): {:.2}%",
        error * 100.0
    );
}
