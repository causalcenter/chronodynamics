use chrono_data_manager::{
    ANOMALOUS_WEEKS, DataManager, YEARS, extract_gps_dataset_id, get_gnss_data_input_path,
    get_year_datasets,
};
use chrono_experiments::print_utils::{print_chrono_mass_header, print_mass_summary};
use chrono_experiments::proces_utils::{apply_mad_filter_points, interpolate_space_time};
use chrono_experiments::statistics_utils::{
    analyze_variance_normalization, calculate_pearson_correlation, calculate_spearman_correlation,
};
use chrono_experiments::*;
use deep_causality_num::{Float106, RealField, ToPrimitive};
use deep_causality_physics::{
    ChronoGauge, ChronoOpsAnalytical, ChronoOpsGauge, ChronoOpsGaugeMut, EARTH_GM,
    SpaceTimeCoordinate,
};
use deep_causality_topology::Lattice;
use rayon::prelude::*;
use std::io;
use std::io::Error;
use std::sync::Arc;
use std::sync::Mutex;

/// Satellite ID for this experiment (Galileo E14).
/// Alterinatively, use E18 for validaton albeit its known to be less precise.
const SAT_ID: &str = "E14";
/// Enable verbose logging
const DBG: bool = false;

/// Active mode for this experiment run. Chose between Analytical and Gauge.
const GM_MODE: GmDeriveMode = GmDeriveMode::Gauge;

/// Change this to `f64` for standard precision or `Float106` for high precision.
pub type FloatType = Float106;

/// Helper Macro to convert f64 literals to FloatType.
macro_rules! flt {
    ($val:expr) => {
        FloatType::from($val)
    };
}

fn main() -> io::Result<()> {
    print_chrono_mass_header(SAT_ID, GM_MODE == GmDeriveMode::Analytical);

    // Verify data path
    let data_path_buf = get_gnss_data_input_path();
    let data_path = data_path_buf.to_str().unwrap();

    // 4D lattice, size 32^4, periodic boundaries
    let lattice = Arc::new(Lattice::new([32, 32, 32, 32], [true, true, true, true]));

    let mut yearly_results: Vec<(String, Vec<GmDataPoint<FloatType>>)> = Vec::new();
    for year in YEARS {
        let results = run_year_analysis(GM_MODE, year, lattice.clone(), data_path)?;
        yearly_results.push((year.to_string(), results));
    }

    let mut master_points: Vec<GmDataPoint<FloatType>> = Vec::new();

    for (year, results) in &yearly_results {
        // Extract GMs for standard summary
        let gms: Vec<FloatType> = results.iter().map(|p| p.gm).collect();
        print_mass_summary(&gms, year);

        master_points.extend(results.iter().cloned());

        // Per year correlation analysis
        analyze_error(results, year);
    }

    let master_gms: Vec<FloatType> = master_points.iter().map(|p| p.gm).collect();
    print_mass_summary(&master_gms, "TOTAL (All Years)");

    // Total correlation analysis
    analyze_error(&master_points, "TOTAL (All Years)");

    Ok(())
}

fn run_year_analysis(
    mode: GmDeriveMode,
    year: &str,
    lattice: Arc<Lattice<4>>,
    data_path: &str,
) -> io::Result<Vec<GmDataPoint<FloatType>>> {
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
) -> io::Result<Vec<GmDataPoint<FloatType>>> {
    let datasets = get_year_datasets(year);
    let global_results = Mutex::new(Vec::new());

    let gauge_field = ChronoGauge::<FloatType>::identity(lattice, flt!(1.0));

    // Process datasets in parallel
    datasets.par_iter().for_each(|dataset| {
        // Extract GPS dataset ID from filename (e.g., "gbm19670" -> 19670)
        if let Some(dataset_id) = extract_gps_dataset_id(dataset)
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

        let process_result = (|| -> Result<Vec<GmDataPoint<FloatType>>, io::Error> {
            let config = AnalysisConfig::default();

            // Load GNSS data
            let dm = DataManager::default();
            let (clocks, orbits) = dm.load_gnss_single_satellite(&clk_path, &sp3_path, SAT_ID)?;

            // Interpolate orbits to clock timestamps using 10th-order Lagrange polynomial
            let data: Vec<SpaceTimeCoordinate<FloatType>> =
                interpolate_space_time(&clocks, &orbits);

            // Skip datasets with insufficient data
            if data.len() <= config.window_size_indices {
                if DBG {
                    println!("  → Insufficient data points for window size");
                }
                return Ok(Vec::new());
            }

            let mut raw_points: Vec<GmDataPoint<FloatType>> = Vec::new();

            let mut i = 0;
            while i < data.len() - config.window_size_indices {
                let idx_a = i;
                let idx_b = i + config.window_size_indices;
                let r_a: FloatType = data[idx_a].r_m;
                let r_b: FloatType = data[idx_b].r_m;
                let d_h: FloatType = (r_a - r_b).abs();

                // Skip if height difference is too small
                if d_h < config.min_height_diff_m {
                    i += 1;
                    continue;
                }

                // Use solve_gm_analytical() to invert the Einstein field equation
                if let Ok(gm) = gauge_field.solve_gm_analytical(&data[idx_a], &data[idx_b]) {
                    // Calculate Radial Velocity v_r = (r . v) / |r|
                    let r_vec_a = data[idx_a].position;
                    let v_vec_a = data[idx_a].velocity;
                    let dot_rv_a =
                        r_vec_a[0] * v_vec_a[0] + r_vec_a[1] * v_vec_a[1] + r_vec_a[2] * v_vec_a[2];
                    let vr_a = dot_rv_a / r_a;

                    let r_vec_b = data[idx_b].position;
                    let v_vec_b = data[idx_b].velocity;
                    let dot_rv_b =
                        r_vec_b[0] * v_vec_b[0] + r_vec_b[1] * v_vec_b[1] + r_vec_b[2] * v_vec_b[2];
                    let vr_b = dot_rv_b / r_b;

                    let vr_avg = (vr_a + vr_b) / flt!(2.0);

                    // Calculate Position Angles (Lat/Long)
                    // Use midpoint position for estimating angles
                    let x = (r_vec_a[0] + r_vec_b[0]) / flt!(2.0);
                    let y = (r_vec_a[1] + r_vec_b[1]) / flt!(2.0);
                    let z = (r_vec_a[2] + r_vec_b[2]) / flt!(2.0);
                    let r = (r_a + r_b) / flt!(2.0);

                    // Latitude phi = asin(z/r)
                    let lat = (z / r).asin();

                    // Longitude lambda = atan2(y, x)
                    let long = y.atan2(x);

                    raw_points.push(GmDataPoint::new(r_a, gm, vr_avg, lat, long));
                }

                i += config.step_size;
            }

            if raw_points.is_empty() {
                println!("  → No valid GM derivations");
                return Ok(Vec::new());
            }

            // Apply MAD filter for outlier rejection (generic over GmDataPoint)
            let filtered = apply_mad_filter_points(&raw_points, flt!(config.outlier_sigma));

            Ok(filtered)
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
) -> io::Result<Vec<GmDataPoint<FloatType>>> {
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

        let process_result = (|| -> Result<Vec<GmDataPoint<FloatType>>, io::Error> {
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
            if let Err(e) = gauge_field.populate_smooth_links_from_source() {
                return Err(Error::other(e));
            }

            // 3. Extract GM via Gauge Field
            match gauge_field.solve_gm() {
                Ok(gm) => {
                    // Outlier detection requested by user
                    let earth_gm = flt!(3.986004418e14);
                    let diff = (gm - earth_gm).abs();
                    let error = diff / earth_gm;

                    if error > flt!(0.5) {
                        println!(
                            "OUTLIER detected in {}: GM={:e} (Error: {:.1}%)",
                            dataset,
                            gm,
                            error * flt!(100.0)
                        );
                    }

                    Ok(vec![GmDataPoint::new(
                        flt!(0.0),
                        gm,
                        flt!(0.0),
                        flt!(0.0),
                        flt!(0.0),
                    )])
                }
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

/// Calculates and logs the correlation between Satellite Altitude and GM Error
fn analyze_error<T>(data: &[GmDataPoint<T>], label: &str)
where
    T: RealField
        + From<f64>
        + PartialOrd
        + Copy
        + std::fmt::LowerExp
        + std::fmt::Display
        + ToPrimitive,
{
    if data.is_empty() {
        return;
    }

    let earth_gm_ref: T = T::from(EARTH_GM);

    // Calculate Errors
    let errors: Vec<T> = data
        .iter()
        .map(|p| (p.gm - earth_gm_ref) / earth_gm_ref)
        .collect();

    // 1. Mean Error
    let sum_err: T = errors.iter().fold(T::zero(), |acc, &e| acc + e);
    let mean_err = sum_err / T::from(errors.len() as f64);

    // 2. Median Error & GM Statistics
    let mut sorted_indices: Vec<usize> = (0..data.len()).collect();
    sorted_indices.sort_by(|&i, &j| {
        data[i]
            .gm
            .partial_cmp(&data[j].gm)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let len = data.len();
    let two = T::from(2.0);

    let median_gm = if len % 2 == 1 {
        data[sorted_indices[len / 2]].gm
    } else {
        let mid = len / 2;
        (data[sorted_indices[mid - 1]].gm + data[sorted_indices[mid]].gm) / two
    };

    // For median error, we can reuse sorted indices if errors are monotonic with GM
    let median_err = (median_gm - earth_gm_ref) / earth_gm_ref;

    let sum_gm: T = data.iter().fold(T::zero(), |acc, p| acc + p.gm);
    let mean_gm = sum_gm / T::from(len as f64);

    // 3. Pearson Correlation (Altitude vs Signed Error)
    let pearson_signed = calculate_pearson_correlation(data, &errors, |p| p.altitude);

    // 4. Spearman Correlation (Altitude vs Signed Error)
    let spearman_signed = calculate_spearman_correlation(data, &errors, |p| p.altitude);

    // 5. Correlation (Altitude vs |GM - Mean|)
    let abs_dev_mean: Vec<T> = data.iter().map(|p| (p.gm - mean_gm).abs()).collect();
    let pearson_abs_mean = calculate_pearson_correlation(data, &abs_dev_mean, |p| p.altitude);

    // 6. Correlation (Altitude vs |GM - Median|)
    let abs_dev_median: Vec<T> = data.iter().map(|p| (p.gm - median_gm).abs()).collect();
    let pearson_abs_median = calculate_pearson_correlation(data, &abs_dev_median, |p| p.altitude);

    // 7. Radial Velocity Correlations
    let pearson_vr_signed = calculate_pearson_correlation(data, &errors, |p| p.radial_velocity);
    let spearman_vr_signed = calculate_spearman_correlation(data, &errors, |p| p.radial_velocity);
    let pearson_vr_abs = calculate_pearson_correlation(data, &abs_dev_mean, |p| p.radial_velocity);

    // 8. Latitude Correlations
    let pearson_lat_signed = calculate_pearson_correlation(data, &errors, |p| p.latitude);
    let spearman_lat_signed = calculate_spearman_correlation(data, &errors, |p| p.latitude);
    let pearson_lat_abs = calculate_pearson_correlation(data, &abs_dev_mean, |p| p.latitude);

    // New: |Latitude| vs Signed Error (Check for J2/Symmetric bias)
    let pearson_abs_lat_signed = calculate_pearson_correlation(data, &errors, |p| p.latitude.abs());

    // 9. Longitude Correlations
    let pearson_long_signed = calculate_pearson_correlation(data, &errors, |p| p.longitude);
    let spearman_long_signed = calculate_spearman_correlation(data, &errors, |p| p.longitude);
    let pearson_long_abs = calculate_pearson_correlation(data, &abs_dev_mean, |p| p.longitude);

    // New: |Longitude| vs Signed Error (Check for Symmetric longitudinal bias)
    let pearson_abs_long_signed =
        calculate_pearson_correlation(data, &errors, |p| p.longitude.abs());

    println!("\n=== GM ERROR ANALYSIS: {} ===", label);
    println!("  Mean Error:          {:+.4e}", mean_err);
    println!("  Median Error:        {:+.4e}", median_err);
    println!("  Pearson (Alt vs Err):      {:+.4}", pearson_signed);
    println!("  Spearman (Alt vs Err):     {:+.4}", spearman_signed);
    println!("  Pearson (Alt vs |DevMean|):   {:+.4}", pearson_abs_mean);
    println!("  Pearson (Alt vs |DevMedian|): {:+.4}", pearson_abs_median);
    println!("  Pearson (Vr vs Err):       {:+.4}", pearson_vr_signed);
    println!("  Spearman (Vr vs Err):      {:+.4}", spearman_vr_signed);
    println!("  Pearson (Vr vs |DevMean|):    {:+.4}", pearson_vr_abs);
    println!("  Pearson (Lat vs Err):      {:+.4}", pearson_lat_signed);
    println!("  Spearman (Lat vs Err):     {:+.4}", spearman_lat_signed);
    println!("  Pearson (Lat vs |DevMean|):   {:+.4}", pearson_lat_abs);
    println!(
        "  Pearson (|Lat| vs Err):    {:+.4}",
        pearson_abs_lat_signed
    );
    println!("  Pearson (Long vs Err):     {:+.4}", pearson_long_signed);
    println!("  Spearman (Long vs Err):    {:+.4}", spearman_long_signed);
    println!("  Pearson (Long vs |DevMean|):  {:+.4}", pearson_long_abs);
    println!(
        "  Pearson (|Long| vs Err):   {:+.4}",
        pearson_abs_long_signed
    );

    // 10. Latitude Variance Normalization
    analyze_variance_normalization(data, &errors);
}
