# Physics Specification: Chrono-Source Integration

**Status:** Draft
**Context:** Unifying `ChronoGauge` Physics with `LatticeGaugeField` Source Architecture.

## 1. Objective

Refactor `deep_causality_physics` to fully leverage the new `Source` ($S$) capability of `LatticeGaugeField`.
Eliminate the method signature `solve_j2(&self, data: &[C])` entirely in favor of `solve_j2(&self)`, where the data comes directly from `self.source()`.

## 2. Type Architecture

### 2.1 ChronoGauge Alias

The `ChronoGauge` alias will now default to carrying data, or be explicit about it.
To support the "Inverse Problem" (Metrology) directly:

```rust
// Proposed Alias in alias/mod.rs
pub type ChronoGauge<FloatType, FieldSource = Vec<SpaceTimeCoordinate<FloatType>>> =
    LatticeGaugeField<SU2_U1, DIM_4D_SPACE_TIME, Complex<FloatType>, FloatType, FieldSource>;
```

This means:
*   `ChronoGauge<f64>` now implies a data-carrying field of type Vec<SpaceTimeCoordinate<FloatType>>.
*   Users can still opt-out with `ChronoGauge<f64, ()>` for pure vacuum if needed, but the default favors the Chrono-Gauge use case (Data -> Geometry).

SpaceTimeCoordinate is defeined 

```rust
// File  types/space_time_coordinate.rs
/// - `R`: Real field type (e.g., `f64`, `DoubleFloat`)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpaceTimeCoordinate<R: RealField> {
    /// The UTC timestamp in seconds (Unix Epoch).
    pub timestamp: u64,
    /// Satellite ID (e.g., E14)
    pub sat_id: u32,
    /// Radius from Earth's center of mass (meters) $|r|$
    pub r_m: R,
    /// Velocity magnitude relative to Earth (m/s) $|v|$
    pub v_ms: R,
    /// Raw clock bias (seconds), uncorrected for relativistic effects.
    pub clock_bias_s: R,
    /// Full 3D position vector [x, y, z] (ITRF frame)
    pub position: [R; 3],
    /// Full 3D velocity vector [vx, vy, vz] (ITRF frame)
    pub velocity: [R; 3],
    /// Rate of change of clock bias (seconds/second), i.e., frequency offset.
    /// This corresponds to $\mathcal{T} - 1$.
    pub clock_drift_rate: R,
}

impl<R: RealField + From<f64>> SpaceTimeCoordinate<R> {
    /// Helper to restore relativistic effects removed by IGS.
    /// Calculates $\Delta t_{periodic} = -2(\vec{r} \cdot \vec{v}) / c^2$
    pub fn get_total_bias(&self) -> R {
        let dot_rv = self.position[0] * self.velocity[0]
            + self.position[1] * self.velocity[1]
            + self.position[2] * self.velocity[2];
        let c_sq = R::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);
        let two = R::from(2.0);
        let rel_correction = -two * dot_rv / c_sq;
        self.clock_bias_s + rel_correction
    }
}

/// Implements the `SpaceTimeCoord` trait for direct use with
/// `solve_gm` and other physics operations.
impl<R: RealField + From<f64>> SpaceTimeCoord<R> for SpaceTimeCoordinate<R> {
    fn clock_drift_rate(&self) -> R {
        self.clock_drift_rate
    }

    fn radius_m(&self) -> R {
        self.r_m
    }

    fn inertial_velocity_magnitude(&self) -> R {
        self.v_ms
    }

    fn z_m(&self) -> R {
        self.position[2]
    }
}
```