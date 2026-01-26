use crate::{EARTH_RADIUS_EQUATORIAL, GEO_ORBIT_RADIUS_M, SpaceTimeCoord};
use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;

// ============================================================================
// Helper functions for GM calculation
// ============================================================================
/// Converts a physical radius to lattice index.
pub(crate) fn radius_to_lattice_index<R>(r: R, n_radial: usize) -> usize
where
    R: RealField + Into<f64>,
{
    let r_min = EARTH_RADIUS_EQUATORIAL;
    let r_max = GEO_ORBIT_RADIUS_M;
    let r_f64: f64 = r.into();
    let r_norm = (r_f64 - r_min) / (r_max - r_min);
    let idx = (r_norm * (n_radial - 1) as f64).round() as usize;
    idx.min(n_radial - 1)
}

/// Converts a lattice radial index back to physical radius.
pub(crate) fn lattice_index_to_radius<R>(idx: usize, n_radial: usize) -> R
where
    R: RealField + From<f64>,
{
    let idx_f64 = idx as f64;
    let n_rad_f64 = (n_radial - 1) as f64;
    let percent = idx_f64 / n_rad_f64;
    let r_min = EARTH_RADIUS_EQUATORIAL;
    let r_max = GEO_ORBIT_RADIUS_M;
    let r_val = r_min + (r_max - r_min) * percent;
    <R as From<f64>>::from(r_val)
}

/// Performs linear interpolation on sorted (x, y) data.
/// Assumes points are sorted by x.
/// Extrapolates using the nearest boundary value (clamping).
pub(crate) fn linear_interpolate_sorted<R>(points: &[(R, R)], x_query: R) -> R
where
    R: RealField + Copy,
{
    if points.is_empty() {
        return R::zero();
    }
    if points.len() == 1 {
        return points[0].1;
    }

    // Handle extrapolation (Clamp)
    if x_query <= points[0].0 {
        return points[0].1;
    }
    if x_query >= points.last().unwrap().0 {
        return points.last().unwrap().1;
    }

    // Interpolate
    for i in 0..points.len() - 1 {
        let (x0, y0) = points[i];
        let (x1, y1) = points[i + 1];

        if x_query >= x0 && x_query <= x1 {
            // Found segment
            if x1 == x0 {
                return y0; // Avoid div by zero
            }
            let t = (x_query - x0) / (x1 - x0);
            return y0 + (y1 - y0) * t;
        }
    }

    // Should be unreachable due to boundary checks
    points.last().unwrap().1
}

// ============================================================================
// Helper functions for J2 calculation (module-level)
// ============================================================================

/// Computes geocentric latitude from z-coordinate and radius
pub(crate) fn compute_latitude<R, C>(coord: &C) -> R
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
pub(crate) fn linear_regression<R>(data: &[(R, R)]) -> Result<(R, R), TopologyError>
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
