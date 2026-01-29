use crate::GmDataPoint;
use deep_causality_num::{RealField, ToPrimitive};

pub fn analyze_variance_normalization<T>(data: &[GmDataPoint<T>], errors: &[T])
where
    T: RealField
        + From<f64>
        + PartialOrd
        + Copy
        + std::fmt::LowerExp
        + std::fmt::Display
        + ToPrimitive,
{
    // Step 1: Bin data by latitude (18 bins of 10 degrees each)
    let lat_bins = 18;
    // We do calculations in f64 for bin indexing to be safe and simple

    // Collect errors per bin
    let mut bin_errors: Vec<Vec<T>> = vec![Vec::new(); lat_bins];
    for (i, p) in data.iter().enumerate() {
        // Map Lat [-PI/2, PI/2] -> [0, 18)
        let lat_val = p.latitude.to_f64().unwrap();
        let pi_val = std::f64::consts::PI;
        let lat_step_val = pi_val / lat_bins as f64;

        let lat_idx = ((lat_val + pi_val / 2.0) / lat_step_val).floor() as usize;
        let lat_idx = lat_idx.min(lat_bins - 1);
        bin_errors[lat_idx].push(errors[i]);
    }

    // Step 2: Calculate sigma per bin
    let mut bin_sigmas: Vec<T> = Vec::with_capacity(lat_bins);
    for bin in &bin_errors {
        if bin.len() < 2 {
            bin_sigmas.push(T::one()); // Default sigma if insufficient data
        } else {
            let n = T::from(bin.len() as f64);
            let mean = bin.iter().fold(T::zero(), |acc, &x| acc + x) / n;
            let variance = bin
                .iter()
                .fold(T::zero(), |acc, &x| acc + (x - mean) * (x - mean))
                / (n - T::one());
            let sigma = variance.sqrt();
            let epsilon = T::from(1e-20);
            if sigma > epsilon {
                bin_sigmas.push(sigma);
            } else {
                bin_sigmas.push(T::one()); // Prevent division by zero
            }
        }
    }

    // Step 3: Normalize errors
    let mut normalized_errors: Vec<T> = Vec::with_capacity(errors.len());
    for (i, p) in data.iter().enumerate() {
        let lat_val = p.latitude.to_f64().unwrap();
        let pi_val = std::f64::consts::PI;
        let lat_step_val = pi_val / lat_bins as f64;

        let lat_idx = ((lat_val + pi_val / 2.0) / lat_step_val).floor() as usize;
        let lat_idx = lat_idx.min(lat_bins - 1);

        let sigma = bin_sigmas[lat_idx];
        normalized_errors.push(errors[i] / sigma);
    }

    // Step 4: Calculate statistics for normalized errors
    let n = T::from(normalized_errors.len() as f64);
    let mean_norm = normalized_errors.iter().fold(T::zero(), |acc, &x| acc + x) / n;

    let mut sorted = normalized_errors.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median_norm = sorted[sorted.len() / 2];

    // Variance of normalized errors (should be ~1.0 if normalization worked)
    let variance_norm = normalized_errors
        .iter()
        .fold(T::zero(), |acc, &x| acc + (x - mean_norm) * (x - mean_norm))
        / (n - T::one());
    let std_norm = variance_norm.sqrt();

    // Re-check latitude correlation after normalization
    let abs_norm: Vec<T> = normalized_errors.iter().map(|&x| x.abs()).collect();
    let pearson_lat_norm = calculate_pearson_correlation(data, &abs_norm, |p| p.latitude);

    println!("\n--- Latitude Variance Normalization ---");
    println!("  Bin Sigmas (South to North):");
    for (i, sigma) in bin_sigmas.iter().enumerate() {
        let lat_center = -90.0 + (i as f64 + 0.5) * 10.0;
        let count = bin_errors[i].len();
        println!(
            "    Lat {:+6.1}°: σ={:.4e}  (n={})",
            lat_center, sigma, count
        );
    }
    println!("  ---- After Normalization ----");
    println!("  Mean (Norm):         {:+.4}", mean_norm);
    println!("  Median (Norm):       {:+.4}", median_norm);
    println!("  Std Dev (Norm):      {:.4}  (target: 1.0)", std_norm);
    println!(
        "  Pearson(Lat vs |NormErr|): {:+.4}  (target: 0.0)",
        pearson_lat_norm
    );
}

pub fn calculate_pearson_correlation<T, F>(data: &[GmDataPoint<T>], values: &[T], extractor: F) -> T
where
    T: RealField + From<f64> + Copy,
    F: Fn(&GmDataPoint<T>) -> T,
{
    let len = T::from(data.len() as f64);
    let x_mean = data.iter().fold(T::zero(), |acc, p| acc + extractor(p)) / len;
    let y_mean = values.iter().fold(T::zero(), |acc, &v| acc + v) / len;

    let mut numerator = T::zero();
    let mut den_x = T::zero();
    let mut den_y = T::zero();

    for (i, p) in data.iter().enumerate() {
        let dx = extractor(p) - x_mean;
        let dy = values[i] - y_mean;

        numerator += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }

    if den_x > T::zero() && den_y > T::zero() {
        numerator / (den_x.sqrt() * den_y.sqrt())
    } else {
        T::zero()
    }
}

pub fn calculate_spearman_correlation<T, F>(
    data: &[GmDataPoint<T>],
    values: &[T],
    extractor: F,
) -> T
where
    T: RealField + From<f64> + PartialOrd + Copy,
    F: Fn(&GmDataPoint<T>) -> T,
{
    let n = data.len();
    if n < 2 {
        return T::zero();
    }

    let ranks_x = rank_data(&data.iter().map(extractor).collect::<Vec<_>>());
    let ranks_y = rank_data(values);

    let n_flt = T::from(n as f64);
    let mean_rank = (n_flt + T::one()) / T::from(2.0);

    let mut numerator = T::zero();
    let mut den_x = T::zero();
    let mut den_y = T::zero();

    for i in 0..n {
        let dx = ranks_x[i] - mean_rank;
        let dy = ranks_y[i] - mean_rank;

        numerator += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }

    if den_x > T::zero() && den_y > T::zero() {
        numerator / (den_x.sqrt() * den_y.sqrt())
    } else {
        T::zero()
    }
}

pub fn rank_data<T>(data: &[T]) -> Vec<T>
where
    T: RealField + From<f64> + PartialOrd + Copy,
{
    let mut indexed: Vec<(usize, T)> = data.iter().copied().enumerate().collect();
    // Sort by value
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![T::zero(); data.len()];
    let mut i = 0;
    while i < data.len() {
        let mut j = i + 1;
        // Check for ties
        while j < data.len() && (indexed[j].1 - indexed[i].1).abs() < T::from(1e-20) {
            j += 1;
        }

        // Average rank for ties
        let rank_sum = (i..j).map(|k| k as f64 + 1.0).sum::<f64>();
        let avg_rank = T::from(rank_sum / (j - i) as f64);

        for k in i..j {
            ranks[indexed[k].0] = avg_rank;
        }
        i = j;
    }
    ranks
}
