use chrono_data_manager::{
    ANOMALOUS_WEEKS, DataManager, YEARS, extract_gps_dataset_id, get_gnss_data_input_path,
    get_year_datasets,
};
use chrono_experiments::print_utils::{print_chrono_mass_header, print_mass_summary};
use chrono_experiments::proces_utils::interpolate_space_time;
use chrono_experiments::*;
use deep_causality_num::Float106;
use deep_causality_physics::{ChronoGauge, ChronoGaugeMutOps, ChronoGaugeOps, SpaceTimeCoordinate};
use deep_causality_topology::Lattice;
use rayon::prelude::*;
use std::io;
use std::io::Error;
use std::sync::Arc;
use std::sync::Mutex;

/// Satellite ID for this experiment (Galileo E14).
const SAT_ID: &str = "E14";
/// Enable verbose logging
const DBG: bool = false;

/// Active mode for this experiment run
const GM_MODE: GmDeriveMode = GmDeriveMode::Gauge;

/// Change this to `f64` for standard precision or `Float106` for high precision.
pub type FloatType = Float106;

/// Macro to convert f64 literals to FloatType.
macro_rules! flt {
    ($val:expr) => {
        FloatType::from($val)
    };
}

fn main() -> io::Result<()> {
    print_chrono_mass_header(GM_MODE == GmDeriveMode::Analytical);

    // Verify data path
    let data_path_buf = get_gnss_data_input_path();
    let data_path = data_path_buf.to_str().unwrap();

    // 4D lattice, size 32^4, periodic boundaries
    let lattice = Arc::new(Lattice::new([32, 32, 32, 32], [true, true, true, true]));

    let mut yearly_results: Vec<(String, Vec<FloatType>)> = Vec::new();
    for year in YEARS {
        let results = run_year_analysis(GM_MODE, year, lattice.clone(), data_path)?;
        yearly_results.push((year.to_string(), results));
    }

    let mut master_results: Vec<FloatType> = Vec::new();
    for (year, results) in &yearly_results {
        print_mass_summary(results, year);
        master_results.extend(results.iter().cloned());
    }

    print_mass_summary(&master_results, "TOTAL (All Years)");

    Ok(())
}

fn run_year_analysis(
    mode: GmDeriveMode,
    year: &str,
    lattice: Arc<Lattice<4>>,
    data_path: &str,
) -> io::Result<Vec<FloatType>> {
    match mode {
        GmDeriveMode::Analytical => run_year_analysis_analytical(year, lattice, data_path),
        GmDeriveMode::Gauge => run_year_analysis_gauge(year, lattice, data_path),
    }
}

/// Analyze a single year of data.
fn run_year_analysis_analytical(
    year: &str,
    lattice: Arc<Lattice<4>>,
    data_path: &str,
) -> io::Result<Vec<FloatType>> {
    let datasets = get_year_datasets(year);
    let global_results = Mutex::new(Vec::new());

    // Process datasets in parallel
    datasets.par_iter().for_each(|dataset| {
        // Extract GPS dataset ID from filename (e.g., "gbm19670" -> 19670)
        if let Some(dataset_id) = extract_gps_dataset_id(dataset)
            // Filter out "Broken" GPS weeks e.g. Clocks crisis 2017, IGS08 to IGS14 transition, or Sept. anomaly 2018
            // See E00 Chrono Experiment for how these were found via chrono forensic
            // and see chrono_data_manager/src/lib.rs for a complete list of filtered out GPS weeks
            && ANOMALOUS_WEEKS.contains(&dataset_id)
        {
            if DBG {
                println!(
                    "[{}] Skipping anomalous dataset {} (2017/2018 Data Crisis)",
                    dataset, dataset_id
                );
            }
            return;
        }

        let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
        let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);

        let process_result = (|| -> Result<Vec<FloatType>, io::Error> {
            // Load GNSS data
            let dm = DataManager::default();
            let (clocks, orbits) = dm.load_gnss_single_satellite(&clk_path, &sp3_path, SAT_ID)?;

            // Interpolate orbits to clock timestamps using 10th-order Lagrange polynomial
            let data: Vec<SpaceTimeCoordinate<FloatType>> =
                interpolate_space_time(&clocks, &orbits);

            // 1. Create field and attach source data
            // Use beta = 1.0 (standard stiffness)
            let mut gauge_field =
                ChronoGauge::<FloatType>::identity(lattice.clone(), flt!(1.0)).with_source(data);

            // 2. Populate link variables from source (clock drift → link phase)
            // This maps the clock drift rates into the U(1) phase of temporal links
            if let Err(e) = gauge_field.populate_links_from_source() {
                return Err(Error::other(e));
            }

            // 3. Extract GM via Gauge field
            match gauge_field.solve_gm_analytical::<SpaceTimeCoordinate<FloatType>>() {
                Ok(gm) => Ok(vec![gm]),
                Err(e) => Err(Error::other(e)),
            }
        })();

        match process_result {
            Ok(results) => {
                let mut acc = global_results.lock().unwrap();
                acc.extend(results);
            }
            Err(e) => {
                eprintln!("Failed to process dataset {}: {}", dataset, e);
            }
        }
    });

    let final_results = global_results.into_inner().unwrap();
    Ok(final_results)
}

/// Analyze a single year of data using Gauge Action (Kinematic Inversion).
fn run_year_analysis_gauge(
    year: &str,
    lattice: Arc<Lattice<4>>,
    data_path: &str,
) -> io::Result<Vec<FloatType>> {
    let datasets = get_year_datasets(year);
    let global_results = Mutex::new(Vec::new());

    // Process datasets in parallel
    datasets.par_iter().for_each(|dataset| {
        // Extract GPS dataset ID from filename (e.g., "gbm19670" -> 19670)
        if let Some(dataset_id) = extract_gps_dataset_id(dataset)
            // Filter out "Broken" GPS weeks
            && ANOMALOUS_WEEKS.contains(&dataset_id)
        {
            if DBG {
                println!(
                    "[{}] Skipping anomalous dataset {} (2017/2018 Data Crisis)",
                    dataset, dataset_id
                );
            }
            return;
        }

        let clk_path = format!("{}/{}/{}.clk", data_path, year, dataset);
        let sp3_path = format!("{}/{}/{}.sp3", data_path, year, dataset);

        let process_result = (|| -> Result<Vec<FloatType>, io::Error> {
            // Load GNSS data
            let dm = DataManager::default();
            let (clocks, orbits) = dm.load_gnss_single_satellite(&clk_path, &sp3_path, SAT_ID)?;

            // Interpolate orbits to clock timestamps using 10th-order Lagrange polynomial
            let data: Vec<SpaceTimeCoordinate<FloatType>> =
                interpolate_space_time(&clocks, &orbits);

            // 1. Create field and attach source data
            let mut gauge_field =
                ChronoGauge::<FloatType>::identity(lattice.clone(), flt!(1.0)).with_source(data);

            // 2. Populate link variables from source
            if let Err(e) = gauge_field.populate_links_from_source() {
                return Err(Error::other(e));
            }

            // 3. Extract GM via Gauge Field
            match gauge_field.solve_gm() {
                Ok(gm) => Ok(vec![gm]),
                Err(e) => Err(Error::other(e)),
            }
        })();

        match process_result {
            Ok(results) => {
                let mut acc = global_results.lock().unwrap();
                acc.extend(results);
            }
            Err(e) => {
                eprintln!("Failed to process dataset {}: {}", dataset, e);
            }
        }
    });

    let final_results = global_results.into_inner().unwrap();
    Ok(final_results)
}
