//! Unified Data Manager for GQCD Experiments
//!
//! This module provides a centralized API for loading all types of data used across
//! GQCD experiments. The DataManager wraps individual data loaders and provides a
//! consistent, discoverable interface.

use std::io;
use std::path::{Path, PathBuf};

use deep_causality_num::RealField;

use crate::{
    BinaryPulsarParams, ClockData, EinsteinElevatorData, GnssDataResult, OmniRecord, OrbitData,
    PulsarObservation, RamanShot, SparcGalaxy,
};

use crate::data_loader::load_pulsar::load_raw_toas_data;
// Import internal loaders
use crate::data_loader::load_pulsar::RawTOA;
use crate::{
    discover_binary_systems, get_ipta_data_path, get_nanograv_data_path, get_sparc_data_path,
    get_zenodo_data_path, list_galaxies, load_all_pta, load_all_rotation_curves,
    load_all_satellites, load_arb_rot_data, load_clock_data, load_data,
    load_einstein_elevator_data, load_galaxy_by_name, load_gm_values_from_csv, load_omni_data,
    load_orbit_data,
};

/// Unified data manager for all GQCD experiment data loading.
///
/// Provides a centralized API for loading:
/// - GNSS satellite data (clock, orbit, single/multi-satellite)
/// - Atom interferometry data
/// - MHD/magnetic storm data
/// - Pulsar timing array data
/// - Binary pulsar data
/// - Derived analysis results
#[derive(Debug, Clone)]
pub struct DataManager {
    /// Base path for GNSS data (optional, uses default if None)
    #[allow(dead_code)] // Reserved for future custom path configuration
    gnss_data_path: Option<PathBuf>,
    /// Base path for Zenodo data (optional, uses default if None)
    zenodo_data_path: Option<PathBuf>,
    /// Base path for IPTA data (optional, uses default if None)
    ipta_data_path: Option<PathBuf>,
    /// Base path for Pulsar NANOGrav data (optional, uses default if None)
    nanograv_data_path: Option<PathBuf>,
    /// Base path for SPARC galaxy data (optional, uses default if None)
    sparc_data_path: Option<PathBuf>,
}

impl Default for DataManager {
    fn default() -> Self {
        Self::new()
    }
}

impl DataManager {
    /// Create a new DataManager with default paths.
    pub fn new() -> Self {
        Self {
            gnss_data_path: None,
            zenodo_data_path: None,
            ipta_data_path: None,
            nanograv_data_path: None,
            sparc_data_path: None,
        }
    }

    /// Create a DataManager with custom paths.
    pub fn with_paths(
        gnss_data_path: Option<PathBuf>,
        zenodo_data_path: Option<PathBuf>,
        ipta_data_path: Option<PathBuf>,
        nanograv_data_path: Option<PathBuf>,
        sparc_data_path: Option<PathBuf>,
        _pulsar_data_path: Option<PathBuf>,
    ) -> Self {
        Self {
            gnss_data_path,
            zenodo_data_path,
            ipta_data_path,
            nanograv_data_path,
            sparc_data_path,
        }
    }

    // ========================================================================
    // GNSS Data Loading Methods
    // ========================================================================

    /// Load GNSS data for a single satellite.
    ///
    /// # Arguments
    /// * `clk_path` - Path to the .clk (clock) file
    /// * `sp3_path` - Path to the .sp3 (orbit) file
    /// * `sat_id` - Satellite ID (e.g., "E14", "E18")
    ///
    /// # Returns
    /// Tuple of (clock_data, orbit_data) vectors
    pub fn load_gnss_single_satellite<R, P>(
        &self,
        clk_path: P,
        sp3_path: P,
        sat_id: &str,
    ) -> GnssDataResult<R>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_data(clk_path, sp3_path, sat_id)
    }

    /// Load GNSS data for all Galileo satellites in the files.
    ///
    /// # Arguments
    /// * `clk_path` - Path to the .clk (clock) file
    /// * `sp3_path` - Path to the .sp3 (orbit) file
    ///
    /// # Returns
    /// Tuple of (clock_data, orbit_data) vectors containing all satellites
    pub fn load_gnss_all_satellites<R, P>(&self, clk_path: P, sp3_path: P) -> GnssDataResult<R>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_all_satellites(clk_path, sp3_path)
    }

    /// Load only GNSS clock data for a single satellite.
    ///
    /// # Arguments
    /// * `clk_path` - Path to the .clk (clock) file
    /// * `sat_id` - Satellite ID (e.g., "E14", "E18")
    ///
    /// # Returns
    /// Vector of clock data
    pub fn load_gnss_clock_data<R, P>(
        &self,
        clk_path: P,
        sat_id: &str,
    ) -> io::Result<Vec<ClockData<R>>>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_clock_data(clk_path, sat_id)
    }

    /// Load only GNSS orbit data for a single satellite.
    ///
    /// # Arguments
    /// * `sp3_path` - Path to the .sp3 (orbit) file
    /// * `sat_id` - Satellite ID (e.g., "E14", "E18")
    ///
    /// # Returns
    /// Vector of orbit data
    pub fn load_gnss_orbit_data<R, P>(
        &self,
        sp3_path: P,
        sat_id: &str,
    ) -> io::Result<Vec<OrbitData<R>>>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_orbit_data(sp3_path, sat_id)
    }

    // ========================================================================
    // Interferometry Data Loading Methods
    // ========================================================================

    /// Load Raman interferometry data from arbitrary orientations dataset.
    ///
    /// # Arguments
    /// * `base_path` - Base path to the arb_rot dataset
    ///
    /// # Returns
    /// Vector of Raman shot measurements
    pub fn load_interferometry_arb_rot<R, P>(&self, base_path: P) -> io::Result<Vec<RamanShot<R>>>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_arb_rot_data(base_path.as_ref())
    }

    /// Load Einstein Elevator interferometry data from MATLAB file.
    ///
    /// # Arguments
    /// * `path` - Path to the .m MATLAB file
    ///
    /// # Returns
    /// Einstein Elevator dataset
    pub fn load_interferometry_einstein_elevator<R, P>(
        &self,
        path: P,
    ) -> io::Result<EinsteinElevatorData<R>>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_einstein_elevator_data(path.as_ref())
    }

    /// Get the default path to Zenodo interferometry data.
    ///
    /// # Returns
    /// PathBuf to the Zenodo data directory
    pub fn get_interferometry_data_path(&self) -> PathBuf {
        self.zenodo_data_path
            .clone()
            .unwrap_or_else(get_zenodo_data_path)
    }

    // ========================================================================
    // MHD Data Loading Methods
    // ========================================================================

    /// Load OMNI magnetic storm data.
    ///
    /// # Arguments
    /// * `path` - Path to the OMNI data file
    ///
    /// # Returns
    /// Vector of OMNI records (Dst, Kp, F10.7 indices)
    pub fn load_mhd_omni_data<R>(&self, path: &str) -> io::Result<Vec<OmniRecord<R>>>
    where
        R: RealField + From<f64>,
    {
        load_omni_data(path)
    }

    // ========================================================================
    // PTA (Pulsar Timing Array) Data Loading Methods
    // ========================================================================

    /// Load all pulsar timing array data from IPTA directory.
    ///
    /// # Arguments
    /// * `dir_path` - Path to the IPTA data directory
    ///
    /// # Returns
    /// Vector of pulsar observations
    pub fn load_pta_all<R>(&self, dir_path: &str) -> io::Result<Vec<PulsarObservation<R>>>
    where
        R: RealField + From<f64>,
    {
        load_all_pta(dir_path)
    }

    /// Get the default path to IPTA data.
    ///
    /// # Returns
    /// PathBuf to the IPTA data directory
    pub fn get_pta_ipta_path(&self) -> PathBuf {
        self.ipta_data_path
            .clone()
            .unwrap_or_else(get_ipta_data_path)
    }

    /// Get the default path to NANOGrav data.
    ///
    /// # Returns
    /// PathBuf to the NANOGrav data directory
    pub fn get_pta_nanograv_path(&self) -> PathBuf {
        self.nanograv_data_path
            .clone()
            .unwrap_or_else(get_nanograv_data_path)
    }

    // ========================================================================
    // Pulsar Data Loading Methods
    // ========================================================================

    /// Discover and load binary pulsar systems from NANOGrav dataset.
    ///
    /// # Arguments
    /// * `dir_path` - Path to the NANOGrav data directory
    ///
    /// # Returns
    /// Vector of binary pulsar parameters
    pub fn load_pulsar_binary_systems<R, P>(
        &self,
        dir_path: P,
    ) -> io::Result<Vec<BinaryPulsarParams<R>>>
    where
        R: RealField + From<f64> + Into<f64>,
        P: AsRef<Path>,
    {
        discover_binary_systems(dir_path.as_ref())
    }

    /// Load raw TOAs for a pulsar from NANOGrav .tim file
    pub fn load_raw_toas<R>(nanograv_dir: &Path, pulsar_name: &str) -> io::Result<Vec<RawTOA<R>>>
    where
        R: RealField + From<f64>,
    {
        load_raw_toas_data(nanograv_dir, pulsar_name)
    }

    // ========================================================================
    // Utility Data Loading Methods
    // ========================================================================

    /// Load derived GM values from CSV file (from gqcd_chrono_mass output).
    ///
    /// # Arguments
    /// * `path` - Path to the CSV file containing GM values
    ///
    /// # Returns
    /// Vector of GM (geocentric gravitational constant) values
    pub fn load_derived_gm_values<R, P>(&self, path: P) -> io::Result<Vec<R>>
    where
        R: RealField + From<f64>,
        P: AsRef<Path>,
    {
        load_gm_values_from_csv(path.as_ref())
    }

    // ========================================================================
    // SPARC Galaxy Data Loading Methods
    // ========================================================================

    /// Load all SPARC galaxy rotation curves.
    ///
    /// # Returns
    /// Vector of all 175 galaxies with rotation curves
    pub fn load_sparc_all_galaxies<R>(&self) -> io::Result<Vec<SparcGalaxy<R>>>
    where
        R: RealField + From<f64> + Clone,
    {
        let path = self
            .sparc_data_path
            .clone()
            .unwrap_or_else(get_sparc_data_path);
        load_all_rotation_curves(&path)
    }

    /// Load a specific SPARC galaxy by name.
    ///
    /// # Arguments
    /// * `galaxy_name` - Galaxy name (e.g., "NGC6503", "DDO154")
    ///
    /// # Returns
    /// The galaxy rotation curve data
    pub fn load_sparc_galaxy<R>(&self, galaxy_name: &str) -> io::Result<SparcGalaxy<R>>
    where
        R: RealField + From<f64>,
    {
        let path = self
            .sparc_data_path
            .clone()
            .unwrap_or_else(get_sparc_data_path);
        load_galaxy_by_name(&path, galaxy_name)
    }

    /// List all available SPARC galaxy names.
    ///
    /// # Returns
    /// Vector of galaxy names
    pub fn list_sparc_galaxies(&self) -> io::Result<Vec<String>> {
        let path = self
            .sparc_data_path
            .clone()
            .unwrap_or_else(get_sparc_data_path);
        list_galaxies(&path)
    }

    /// Get the default path to SPARC data.
    ///
    /// # Returns
    /// PathBuf to the SPARC data directory
    pub fn get_sparc_data_path(&self) -> PathBuf {
        self.sparc_data_path
            .clone()
            .unwrap_or_else(get_sparc_data_path)
    }
}
