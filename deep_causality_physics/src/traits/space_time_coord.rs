use deep_causality_num::RealField;

/// Trait for types that can provide space-time coordinate data.
pub trait SpaceTimeCoord<T: RealField> {
    /// Returns the clock drift rate (proper time / coordinate time - 1).
    fn clock_drift_rate(&self) -> T;

    /// Returns the radial distance from Earth's center in meters.
    fn radius_m(&self) -> T;

    /// Returns the inertial velocity magnitude in m/s.
    fn inertial_velocity_magnitude(&self) -> T;

    /// Returns the Z coordinate (ECEF) in meters.
    fn z_m(&self) -> T;
}
