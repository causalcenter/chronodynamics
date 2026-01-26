use deep_causality_num::RealField;
use deep_causality_topology::TopologyError;

/// Mutable operations for populating the gauge field from source data.
pub trait ChronoOpsGaugeMut<R: RealField> {
    /// Populates temporal link variables from the internal source data.
    ///
    /// This method reads clock drift rates from `self.source()`, bins them
    /// by radius, and encodes the average drift into temporal link phases:
    ///
    /// $$U_0(x) = \exp(i \cdot (1 - \dot{\tau}) \cdot N_t)$$
    ///
    /// # Requirements
    ///
    /// The field must have source data attached (S = Vec<SpaceTimeCoordinate>).
    fn populate_links_from_source(&mut self) -> Result<(), TopologyError>;

    /// Populates temporal link variables using interpolation for smoothness.
    ///
    /// Unlike `populate_links_from_source` which uses binning (creating potential steps),
    /// this method generates a smooth observable function Drift(r) from the data points
    /// and assigns phases based on this continuous function.
    ///
    /// This is essential for gradient-based observables like the Wilson Action
    /// to avoid artificial spikes at bin boundaries.
    fn populate_smooth_links_from_source(&mut self) -> Result<(), TopologyError>;
}
