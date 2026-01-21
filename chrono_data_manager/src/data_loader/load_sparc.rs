//! SPARC Data Loader
//!
//! Loads galaxy rotation curves from the SPARC (Spitzer Photometry and Accurate
//! Rotation Curves) dataset.
//!
//! Reference: Lelli, McGaugh, Schombert (2016). AJ 152, 157.
//! Data source: https://zenodo.org/records/16284118

use crate::{GalaxyProperties, GalaxyRotationCurve, RotationCurvePoint, SparcGalaxy};
use deep_causality_num::RealField;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Load a single galaxy rotation curve from a .dat file.
///
/// # Arguments
/// * `path` - Path to the rotation curve file (e.g., NGC6503_rotmod.dat)
///
/// # Returns
/// The galaxy rotation curve data
pub fn load_rotation_curve<R, P>(path: P) -> io::Result<GalaxyRotationCurve<R>>
where
    R: RealField + From<f64>,
    P: AsRef<Path>,
{
    let path = path.as_ref();
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Extract galaxy name from filename
    let filename = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("Unknown");
    let name = filename.replace("_rotmod", "");

    let mut distance_mpc_f64 = 0.0;
    let mut points = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        // Parse header for distance
        if trimmed.starts_with("# Distance") {
            // Format: "# Distance = 6.26 Mpc"
            if let Some(dist_str) = trimmed.split('=').nth(1) {
                let dist_str = dist_str.trim().replace("Mpc", "").trim().to_string();
                distance_mpc_f64 = dist_str.parse().unwrap_or(0.0);
            }
            continue;
        }

        // Skip other comment lines
        if trimmed.starts_with('#') || trimmed.is_empty() {
            continue;
        }

        // Parse data line
        // Format: Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 8
            && let (Ok(r), Ok(v), Ok(e), Ok(vg), Ok(vd), Ok(vb), Ok(sd), Ok(sb)) = (
                parts[0].parse::<f64>(),
                parts[1].parse::<f64>(),
                parts[2].parse::<f64>(),
                parts[3].parse::<f64>(),
                parts[4].parse::<f64>(),
                parts[5].parse::<f64>(),
                parts[6].parse::<f64>(),
                parts[7].parse::<f64>(),
            )
        {
            points.push(RotationCurvePoint {
                radius_kpc: R::from(r),
                v_obs: R::from(v),
                v_err: R::from(e),
                v_gas: R::from(vg),
                v_disk: R::from(vd),
                v_bul: R::from(vb),
                sb_disk: R::from(sd),
                sb_bul: R::from(sb),
            });
        }
    }

    Ok(GalaxyRotationCurve {
        name,
        distance_mpc: R::from(distance_mpc_f64),
        points,
    })
}

/// Load all galaxy rotation curves from the Rotmod_LTG directory.
///
/// # Arguments
/// * `base_path` - Path to the SPARC data directory (containing Rotmod_LTG/)
///
/// # Returns
/// Vector of all galaxy rotation curves with properties
pub fn load_all_rotation_curves<R, P>(base_path: P) -> io::Result<Vec<SparcGalaxy<R>>>
where
    R: RealField + From<f64> + Clone,
    P: AsRef<Path>,
{
    let base_path = base_path.as_ref();
    let rotmod_dir = base_path.join("Rotmod_LTG");

    if !rotmod_dir.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Rotmod_LTG directory not found at {:?}", rotmod_dir),
        ));
    }

    // Load galaxy properties from catalog
    let properties_map: HashMap<String, GalaxyProperties<R>> =
        load_galaxy_properties(base_path).unwrap_or_default();

    let mut galaxies = Vec::new();

    for entry in fs::read_dir(&rotmod_dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.extension().and_then(|s| s.to_str()) == Some("dat") {
            match load_rotation_curve(&path) {
                Ok(curve) => {
                    let galaxy_name = curve.name.clone();
                    let props = properties_map.get(&galaxy_name).cloned();
                    galaxies.push(SparcGalaxy {
                        properties: props,
                        rotation_curve: curve,
                    });
                }
                Err(e) => {
                    eprintln!("Warning: Failed to load {:?}: {}", path, e);
                }
            }
        }
    }

    // Sort by galaxy name for consistent ordering
    galaxies.sort_by(|a, b| a.name().cmp(b.name()));

    // Debug: count how many galaxies have properties
    let with_props = galaxies.iter().filter(|g| g.properties.is_some()).count();
    eprintln!(
        "Debug: {} of {} galaxies have properties loaded",
        with_props,
        galaxies.len()
    );

    Ok(galaxies)
}

/// Load galaxy properties from the SPARC catalog file (SPARC_Lelli2016c.mrt).
///
/// # Arguments
/// * `base_path` - Path to the SPARC data directory
///
/// # Returns
/// HashMap of galaxy name -> GalaxyProperties
fn load_galaxy_properties<R, P>(base_path: P) -> io::Result<HashMap<String, GalaxyProperties<R>>>
where
    R: RealField + From<f64>,
    P: AsRef<Path>,
{
    let catalog_path = base_path.as_ref().join("SPARC_Lelli2016c.mrt");

    if !catalog_path.exists() {
        return Ok(HashMap::new());
    }

    let file = File::open(&catalog_path)?;
    let reader = BufReader::new(file);
    let mut properties = HashMap::new();

    for line in reader.lines() {
        let line = line?;

        // Skip header lines (starting with space, letter, or -)
        if line.is_empty()
            || line.starts_with(' ')
                && !line
                    .trim()
                    .chars()
                    .next()
                    .map(|c| c.is_alphanumeric())
                    .unwrap_or(false)
        {
            continue;
        }

        // Skip lines that are too short or header lines
        if line.len() < 97
            || line.contains("Galaxy")
            || line.starts_with('-')
            || line.starts_with('=')
        {
            continue;
        }

        // Safety check for line length
        if line.len() < 99 {
            continue;
        }

        let name = line[0..11].trim().to_string();
        if name.is_empty() {
            continue;
        }

        let hubble_type = line
            .get(11..14)
            .and_then(|s| s.trim().parse::<i32>().ok())
            .unwrap_or(0);
        let distance = line
            .get(14..20)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let distance_err = line
            .get(20..25)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let inclination = line
            .get(26..31)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let inclination_err = line
            .get(31..35)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let luminosity = line
            .get(35..42)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let r_eff = line
            .get(49..54)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let r_disk = line
            .get(62..67)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let m_hi = line
            .get(75..82)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let v_flat = line
            .get(87..92)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let v_flat_err = line
            .get(92..97)
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(0.0);
        let quality = line
            .get(97..100)
            .and_then(|s| s.trim().parse::<i32>().ok())
            .unwrap_or(0);

        properties.insert(
            name.clone(),
            GalaxyProperties {
                name,
                hubble_type,
                distance_mpc: R::from(distance),
                distance_err: R::from(distance_err),
                inclination: R::from(inclination),
                inclination_err: R::from(inclination_err),
                luminosity_3p6: R::from(luminosity),
                r_eff: R::from(r_eff),
                r_disk: R::from(r_disk),
                m_hi: R::from(m_hi),
                v_flat: R::from(v_flat),
                v_flat_err: R::from(v_flat_err),
                quality,
            },
        );
    }

    Ok(properties)
}

/// Load a specific galaxy by name.
///
/// # Arguments
/// * `base_path` - Path to the SPARC data directory
/// * `galaxy_name` - Galaxy name (e.g., "NGC6503")
///
/// # Returns
/// The galaxy rotation curve if found
pub fn load_galaxy_by_name<R, P>(base_path: P, galaxy_name: &str) -> io::Result<SparcGalaxy<R>>
where
    R: RealField + From<f64>,
    P: AsRef<Path>,
{
    let filename = format!("{}_rotmod.dat", galaxy_name);
    let path = base_path.as_ref().join("Rotmod_LTG").join(&filename);

    if !path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Galaxy {} not found at {:?}", galaxy_name, path),
        ));
    }

    let curve = load_rotation_curve(&path)?;
    Ok(SparcGalaxy {
        properties: None,
        rotation_curve: curve,
    })
}

/// Get list of all available galaxy names.
///
/// # Arguments
/// * `base_path` - Path to the SPARC data directory
///
/// # Returns
/// Vector of galaxy names
pub fn list_galaxies<P: AsRef<Path>>(base_path: P) -> io::Result<Vec<String>> {
    let rotmod_dir = base_path.as_ref().join("Rotmod_LTG");

    if !rotmod_dir.exists() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "Rotmod_LTG directory not found",
        ));
    }

    let mut names = Vec::new();

    for entry in fs::read_dir(&rotmod_dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.extension().and_then(|s| s.to_str()) == Some("dat")
            && let Some(name) = path.file_stem().and_then(|s| s.to_str())
        {
            let galaxy_name = name.replace("_rotmod", "");
            names.push(galaxy_name);
        }
    }

    names.sort();
    Ok(names)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::get_sparc_data_path;

    #[test]
    fn test_get_sparc_data_path() {
        let path = get_sparc_data_path();
        assert!(path.to_str().unwrap().contains("sparc"));
    }

    #[test]
    fn test_list_galaxies() {
        let path = get_sparc_data_path();
        let galaxies = list_galaxies(&path).expect("Failed to list galaxies");
        assert_eq!(
            galaxies.len(),
            175,
            "Expected 175 galaxies in SPARC dataset"
        );
        assert!(galaxies.contains(&"NGC6503".to_string()));
    }

    #[test]
    fn test_load_ngc6503() {
        let path = get_sparc_data_path();
        let galaxy: SparcGalaxy<f64> =
            load_galaxy_by_name(&path, "NGC6503").expect("Failed to load NGC6503");

        assert_eq!(galaxy.name(), "NGC6503");
        assert!(
            (galaxy.distance() - 6.26).abs() < 0.01,
            "Distance should be 6.26 Mpc"
        );
        assert!(
            galaxy.curve().len() > 20,
            "NGC6503 should have many data points"
        );

        // Check first data point
        let first = &galaxy.curve().points[0];
        assert!((first.radius_kpc - 0.76).abs() < 0.01);
        assert!((first.v_obs - 77.0).abs() < 0.1);
    }

    #[test]
    fn test_load_all_galaxies() {
        let path = get_sparc_data_path();
        let galaxies: Vec<SparcGalaxy<f64>> =
            load_all_rotation_curves(&path).expect("Failed to load all galaxies");
        assert_eq!(galaxies.len(), 175, "Expected 175 galaxies");

        // Verify each galaxy has data
        for galaxy in &galaxies {
            assert!(
                !galaxy.curve().is_empty(),
                "Galaxy {} has no data",
                galaxy.name()
            );
        }
    }

    #[test]
    fn test_rotation_curve_calculations() {
        let path = get_sparc_data_path();
        let galaxy: SparcGalaxy<f64> =
            load_galaxy_by_name(&path, "NGC6503").expect("Failed to load NGC6503");

        // Test v_flat calculation
        let v_flat = galaxy.curve().v_flat();
        assert!(
            v_flat > 100.0 && v_flat < 120.0,
            "NGC6503 v_flat should be ~115 km/s"
        );

        // Test baryonic velocity
        let v_bar = galaxy.curve().v_baryonic();
        assert!(!v_bar.is_empty());

        // Test velocity anomaly (should be positive for Dark Matter signal)
        let v_anom = galaxy.curve().v_anomaly();
        // At outer radii, observed should exceed baryonic (the DM problem!)
        let outer_anomaly = v_anom.last().unwrap();
        assert!(
            *outer_anomaly > 0.0,
            "Outer velocity anomaly should be positive"
        );
    }
}
