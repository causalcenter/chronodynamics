//! SPARC Galaxy Data Types
//!
//! Data structures for the SPARC (Spitzer Photometry and Accurate Rotation Curves) dataset.
//! Contains 175 disk galaxies with rotation curves and mass models.
//!
//! Reference: Lelli, McGaugh, Schombert (2016). AJ 152, 157.

use deep_causality_num::RealField;

/// A single point on a galaxy rotation curve.
#[derive(Debug, Clone)]
pub struct RotationCurvePoint<R>
where
    R: RealField,
{
    /// Galactocentric radius [kpc]
    pub radius_kpc: R,
    /// Observed rotation velocity [km/s]
    pub v_obs: R,
    /// Error on observed velocity [km/s]
    pub v_err: R,
    /// Gas contribution to rotation velocity [km/s]
    pub v_gas: R,
    /// Disk contribution to rotation velocity [km/s]
    pub v_disk: R,
    /// Bulge contribution to rotation velocity [km/s]
    pub v_bul: R,
    /// Disk surface brightness [L☉/pc²]
    pub sb_disk: R,
    /// Bulge surface brightness [L☉/pc²]
    pub sb_bul: R,
}

impl<R> RotationCurvePoint<R>
where
    R: RealField + Clone + From<f64>,
{
    /// Compute the total baryonic velocity: sqrt(v_gas² + v_disk² + v_bul²)
    pub fn v_baryonic(&self) -> R {
        let v_gas2 = self.v_gas * self.v_gas;
        let v_disk2 = self.v_disk * self.v_disk;
        let v_bul2 = self.v_bul * self.v_bul;
        (v_gas2 + v_disk2 + v_bul2).sqrt()
    }

    /// Compute the velocity anomaly (observed - baryonic)
    pub fn v_anomaly(&self) -> R {
        self.v_obs - self.v_baryonic()
    }

    /// Compute the Newtonian acceleration from baryonic matter [m/s²]
    pub fn a_baryonic(&self) -> R {
        let kpc_to_m = R::from(3.086e19_f64);
        let km_s_to_m_s = R::from(1000.0_f64);
        let v_bar = self.v_baryonic() * km_s_to_m_s;
        let r = self.radius_kpc * kpc_to_m;
        if r > R::zero() {
            v_bar * v_bar / r
        } else {
            R::zero()
        }
    }

    /// Compute the observed centripetal acceleration [m/s²]
    pub fn a_observed(&self) -> R {
        let kpc_to_m = R::from(3.086e19_f64);
        let km_s_to_m_s = R::from(1000.0_f64);
        let v = self.v_obs * km_s_to_m_s;
        let r = self.radius_kpc * kpc_to_m;
        if r > R::zero() { v * v / r } else { R::zero() }
    }
}

/// Complete rotation curve for a galaxy.
#[derive(Debug, Clone)]
pub struct GalaxyRotationCurve<R>
where
    R: RealField,
{
    /// Galaxy name (e.g., "NGC6503")
    pub name: String,
    /// Distance to galaxy [Mpc]
    pub distance_mpc: R,
    /// Rotation curve data points
    pub points: Vec<RotationCurvePoint<R>>,
}

impl<R> GalaxyRotationCurve<R>
where
    R: RealField + Clone + From<f64>,
{
    /// Get radii as vector [kpc]
    pub fn radii(&self) -> Vec<R> {
        self.points.iter().map(|p| p.radius_kpc).collect()
    }

    /// Get observed velocities as vector [km/s]
    pub fn v_observed(&self) -> Vec<R> {
        self.points.iter().map(|p| p.v_obs).collect()
    }

    /// Get baryonic velocities as vector [km/s]
    pub fn v_baryonic(&self) -> Vec<R> {
        self.points.iter().map(|p| p.v_baryonic()).collect()
    }

    /// Get velocity anomalies as vector [km/s]
    pub fn v_anomaly(&self) -> Vec<R> {
        self.points.iter().map(|p| p.v_anomaly()).collect()
    }

    /// Get observed accelerations [m/s²]
    pub fn a_observed(&self) -> Vec<R> {
        self.points.iter().map(|p| p.a_observed()).collect()
    }

    /// Get baryonic accelerations [m/s²]
    pub fn a_baryonic(&self) -> Vec<R> {
        self.points.iter().map(|p| p.a_baryonic()).collect()
    }

    /// Get the asymptotically flat velocity (average of outer 3 points)
    pub fn v_flat(&self) -> R {
        if self.points.len() < 3 {
            return self.points.last().map(|p| p.v_obs).unwrap_or(R::zero());
        }
        let n = self.points.len();
        let outer: R = self.points[n - 3..]
            .iter()
            .map(|p| p.v_obs)
            .fold(R::zero(), |acc, x| acc + x);
        outer / R::from(3.0)
    }

    /// Get maximum radius [kpc]
    pub fn r_max(&self) -> R {
        self.points
            .last()
            .map(|p| p.radius_kpc)
            .unwrap_or(R::zero())
    }

    /// Number of data points
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Galaxy properties from the SPARC catalog.
#[derive(Debug, Clone)]
pub struct GalaxyProperties<R>
where
    R: RealField,
{
    /// Galaxy name
    pub name: String,
    /// Hubble type (numeric)
    pub hubble_type: i32,
    /// Distance [Mpc]
    pub distance_mpc: R,
    /// Distance error [Mpc]
    pub distance_err: R,
    /// Inclination [degrees]
    pub inclination: R,
    /// Inclination error [degrees]
    pub inclination_err: R,
    /// Total luminosity at 3.6μm [10⁹ L☉]
    pub luminosity_3p6: R,
    /// Effective radius [kpc]
    pub r_eff: R,
    /// Disk scale length [kpc]
    pub r_disk: R,
    /// Total HI mass [10⁹ M☉]
    pub m_hi: R,
    /// Flat rotation velocity [km/s]
    pub v_flat: R,
    /// Flat velocity error [km/s]
    pub v_flat_err: R,
    /// Quality flag
    pub quality: i32,
}

/// Complete SPARC galaxy with properties and rotation curve.
#[derive(Debug, Clone)]
pub struct SparcGalaxy<R>
where
    R: RealField,
{
    /// Galaxy properties from catalog
    pub properties: Option<GalaxyProperties<R>>,
    /// Rotation curve data
    pub rotation_curve: GalaxyRotationCurve<R>,
}

impl<R> SparcGalaxy<R>
where
    R: RealField + Clone,
{
    /// Get galaxy name
    pub fn name(&self) -> &str {
        &self.rotation_curve.name
    }

    /// Get distance [Mpc]
    pub fn distance(&self) -> R {
        self.rotation_curve.distance_mpc
    }

    /// Get rotation curve reference
    pub fn curve(&self) -> &GalaxyRotationCurve<R> {
        &self.rotation_curve
    }
}
