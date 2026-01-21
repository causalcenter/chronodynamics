/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */
use super::load_clk::parse_clk;
use super::load_sp3::parse_sp3;
use crate::GnssDataResult;
use deep_causality_num::RealField;

use std::path::Path;

pub fn load_data<R, P>(clk_path: P, sp3_path: P, target_sat: &str) -> GnssDataResult<R>
where
    R: RealField + From<f64>,
    P: AsRef<Path>,
{
    let clocks = parse_clk(clk_path, target_sat)?;
    let orbits = parse_sp3(sp3_path, target_sat)?;
    Ok((clocks, orbits))
}
