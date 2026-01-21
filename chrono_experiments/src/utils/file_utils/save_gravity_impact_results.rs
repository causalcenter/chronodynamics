/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */
use crate::MonthAnalysis;
use deep_causality_physics::NEWTONIAN_CONSTANT_OF_GRAVITATION as G;
use std::fs::File;
use std::io::{self, Write};

pub fn save_gravity_impact_month_csv(
    analysis: &MonthAnalysis,
    year_out_path: &str,
) -> io::Result<()> {
    let filename = format!(
        "{}/{}_daily_analysis.csv",
        year_out_path,
        analysis.name().to_lowercase()
    );
    let mut file = File::create(&filename)?;

    writeln!(file, "Day,Mean_G,Deviation_Pct,Std_Dev,Min_G,Max_G,Count")?;

    for day in analysis.days() {
        writeln!(
            file,
            "{},{:.10e},{:.6},{:.10e},{:.10e},{:.10e},{}",
            day.day(),
            day.mean(),
            day.deviation_pct(G),
            day.std_dev(),
            day.min(),
            day.max(),
            day.count()
        )?;
    }

    println!("  → Saved: {}", filename);
    Ok(())
}

pub fn save_gravity_impact_summary_csv(
    analyses: &[MonthAnalysis],
    year: &str,
    year_out_path: &str,
) -> io::Result<()> {
    let filename = format!("{}/{}_lunar_tidal_summary.csv", year_out_path, year);
    let mut file = File::create(&filename)?;

    writeln!(
        file,
        "Month,Quarter,Total_Points,Mean_G,Std_Dev,Peak_Deviation_Pct,Peak_Day,Trough_Deviation_Pct,Trough_Day,Amplitude_Pct,Detected_Period_Days,Derived_LunDist_m,Theory_LunDist_m"
    )?;

    for a in analyses {
        writeln!(
            file,
            "{},{},{},{:.10e},{:.10e},{:.6},{},{:.6},{},{:.6},{:.2},{:.4e},{:.4e}",
            a.name(),
            a.quarter(),
            a.total_points(),
            a.overall_mean(),
            a.overall_std_dev(),
            a.peak_deviation(),
            a.peak_day(),
            a.trough_deviation(),
            a.trough_day(),
            a.amplitude(),
            a.detected_period_days(),
            a.derived_lunar_distance(),
            a.theoretical_lunar_distance()
        )?;
    }

    println!("  → Saved: {}", filename);
    Ok(())
}
