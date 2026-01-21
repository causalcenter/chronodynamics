/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

//! Lattice construction utilities for chrono-gauge experiments.
//!
//! This module provides functions to construct ChronoGauge lattices from
//! GNSS satellite data.

use crate::SpaceTimeCoordinate;
use deep_causality_num::{Field, Float, FromPrimitive, RealField, ToPrimitive};
use deep_causality_physics::ChronoGauge;
use deep_causality_topology::{Lattice, LatticeGaugeField};
use std::sync::Arc;

/// Default lattice shape for per-epoch construction.
/// 8^4 = 4096 sites provides good resolution without excessive memory.
pub const DEFAULT_LATTICE_SHAPE: [usize; 4] = [8, 8, 8, 8];

/// Creates a ChronoGauge lattice for a single epoch from satellite data.
///
/// # Arguments
///
/// * `satellites` - Satellites with space-time coordinates for this epoch
///
/// # Returns
///
/// A `ChronoGauge<R>` field initialized with identity configuration.
///
/// # Note
///
/// The lattice is initialized with identity links. Satellite data embedding
/// is handled separately at the observable computation level rather than
/// modifying link variables directly, as the ChronoGaugeOps methods extract
/// physics from the satellite coordinates.
pub fn create_epoch_lattice<R>(_satellites: &[SpaceTimeCoordinate<R>]) -> ChronoGauge<R>
where
    R: RealField
        + Clone
        + From<f64>
        + Into<f64>
        + Field
        + Default
        + std::fmt::Debug
        + Float
        + FromPrimitive
        + ToPrimitive,
{
    let shape = DEFAULT_LATTICE_SHAPE;
    let periodic = [true, true, true, true];
    let lattice = Arc::new(Lattice::new(shape, periodic));

    // Initialize with identity configuration (β = 1.0)
    // The satellite data is used at observable computation time, not
    // embedded into the lattice structure directly. This matches how
    // ChronoGaugeOps methods like mass_density_action() work.
    let gauge: ChronoGauge<R> = LatticeGaugeField::identity(lattice, R::one());

    gauge
}

/// Returns the number of lattice sites for the default shape.
pub fn default_lattice_sites() -> usize {
    DEFAULT_LATTICE_SHAPE.iter().product()
}
