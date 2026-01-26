use crate::SpaceTimeCoord;
use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;

pub trait ChronoOpsAnalytical<R: RealField> {
    /// Derives GM from source chrono data (satellite clocks) using kinetic energy corrections
    /// without the lattice gauge field.
    ///
    /// This method performs a regression analysis on the clock drift rates of satellites,
    /// correcting for their kinematic state (velocity) to isolate the gravitational potential.
    ///
    /// # Physics
    ///
    /// The clock drift rate $d\tau/dt$ is related to the potential $\Phi$ and velocity $v$ by:
    /// $$ \frac{d\tau}{dt} \approx 1 + \frac{\Phi}{c^2} - \frac{v^2}{2c^2} $$
    ///
    /// By isolating $\Phi = -GM/r$, we can solve for GM:
    /// $$ \Phi \approx c^2 \left( \frac{d\tau}{dt} - 1 + \frac{v^2}{2c^2} \right) $$
    ///
    /// # Returns
    ///
    /// The gravitational parameter GM in m³/s².
    fn solve_gm_analytical<C>(&self, coord_a: &C, coord_b: &C) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>;

    /// Computes J2 oblateness coefficient analytically (without the Gauge Field)
    /// from observed satellite time data.
    fn solve_j2_analytical<C>(&self, data: &[C]) -> Result<R, TopologyError>
    where
        C: SpaceTimeCoord<R>;
}
