# E02 Chrono-Gauge Experiment: Product Specification

**Author:** AI Assistant  
**Date:** 2026-01-18  
**Status:** Implemented  
**Version:** 1.1

---

## 1. Executive Summary

This specification defines the **E02 Chrono-Gauge Experiment** (`chrono_gauge`), which ports the legacy
`gqcd_chrono_manifold` validation experiment to use the Chrono-Gauge Lattice Theory (CGLT) infrastructure.

The experiment validates chrono-gravitational physics by computing lattice gauge observables from GNSS satellite
constellation data and comparing them with theoretical predictions.

### 1.1 Key Value Proposition

| Legacy (gqcd_chrono_manifold)     | New (chrono_gauge)                        |
|-----------------------------------|-------------------------------------------|
| Discrete Exterior Calculus (DEC)  | Lattice Gauge Theory (LGT)                |
| `Manifold<f64>` with scalar field | `ChronoGauge<DoubleFloat>` gauge field    |
| Manual Laplacian/curl computation | Wilson action / plaquette traces          |
| Single-shot observables           | Monte Carlo error estimates               |
| `f64` precision (~15 digits)      | `DoubleFloat` quad-precision (~31 digits) |
| 4-point tetrahedron per epoch     | Full 4D hypercubic lattice                |

### 1.2 Experiment Classification

| Attribute       | Value                                      |
|-----------------|--------------------------------------------|
| **Experiment**  | E02_chrono_gauge                           |
| **Binary**      | `chrono_gauge`                             |
| **Theory**      | Chrono-Gauge Lattice Theory (U(1) × SU(2)) |
| **Data**        | GNSS CLK/SP3 files (2016-2018)             |
| **Observables** | 7 validation metrics + PCA                 |

---

## 2. Problem Statement

### 2.1 Legacy Experiment Overview

The `gqcd_chrono_manifold` experiment (443 lines) validates the chronometric field using discrete exterior calculus:

```text
gqcd/bin/gqcd_chrono_manifold/main.rs
├── Load GNSS satellite data (CLK + SP3)
├── Group satellites by epoch (≥6 satellites required)
├── Build ChronoManifold via Vietoris-Rips triangulation
├── Compute DEC observables:
│   ├── laplacian_chrono() → mass/tidal validation
│   ├── curl() → frame-dragging validation
│   ├── gradient() → geoid validation
│   └── hodge_star() → vorticity
├── Collect physics correlations:
│   ├── Time-velocity correlation (momentum)
│   ├── Virial ratio (energy)
│   ├── Tolman consistency (thermodynamics)
│   └── Action-phase correlation (wave)
└── PCA analysis on 7 combined metrics
```

### 2.2 Limitations of DEC Approach

1. **4-point limitation**: Only 4 satellites used per epoch (maximally dispersed tetrahedron)
2. **No gauge structure**: Scalar field cannot capture gravitomagnetic effects properly
3. **Uncontrolled discretization**: $O(a^2)$ errors without systematic improvement
4. **No error bars**: Single-shot computation, no ensemble averaging
5. **Precision ceiling**: `f64` insufficient for next-generation clocks

### 2.3 CGLT Advantages

The Chrono-Gauge Lattice Theory infrastructure provides:

| Capability         | Method                  | Benefit                         |
|--------------------|-------------------------|---------------------------------|
| Wilson action      | `mass_density_action()` | Mass/energy with O(a⁴) Symanzik |
| SU(2) plaquettes   | `vorticity_tensor()`    | True gravitomagnetic effects    |
| Polyakov loops     | `fleet_coherence()`     | Temporal coherence measure      |
| Topological charge | `topological_charge()`  | Non-perturbative field topology |
| Monte Carlo        | `measure_with_error()`  | Jackknife error estimation      |
| DoubleFloat        | Generic `R: RealField`  | Quad-precision arithmetic       |

---

## 3. Technical Design

### 3.1 Observable Mapping: DEC → LGT

| # | Observable Name    | DEC Method (Legacy)       | LGT Method (New)          | Physics Interpretation         |
|---|--------------------|---------------------------|---------------------------|--------------------------------|
| 1 | Mass/Tidal         | `laplacian_chrono()`      | `mass_density_action()`   | ρ ∝ ∆𝒯 (Poisson equation)     |
| 2 | Frame-Dragging     | `curl()`                  | `vorticity_tensor()` norm | Gravitomagnetic field strength |
| 3 | Geoid              | `gradient()` magnitude    | U(1) link phases          | Gravitational field direction  |
| 4 | Spin (Vorticity-Z) | `curl()[2]`               | `vorticity_z()`           | Earth rotation axis alignment  |
| 5 | Momentum           | Time-velocity correlation | Polyakov susceptibility   | ∂𝒯/∂v relationship            |
| 6 | Energy             | Virial ratio              | `virial_ratio()`          | \|Kinetic\| / \|Potential\|    |
| 7 | Thermodynamics     | Tolman consistency        | `tolman_temperature()`    | T × τ = constant               |

### 3.2 Data Pipeline

```text
GNSS Data (CLK + SP3)
        │
        ▼
┌───────────────────────────────────────────────────┐
│  DataManager::load_gnss_all_satellites()          │
│  → Vec<ClockData>, Vec<OrbitData>                 │
└───────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────┐
│  interpolate_space_time()                         │
│  → Vec<SpaceTimeCoordinate<FloatType>>            │
│  (10th-order Lagrange interpolation)              │
└───────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────┐
│  Group by timestamp                               │
│  → HashMap<i64, Vec<SpaceTimeCoordinate>>         │
│  (epochs with ≥6 satellites)                      │
└───────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────┐
│  Per-Epoch LGT Construction                       │
│  1. Create 4D Lattice from satellite positions    │
│  2. Initialize ChronoGauge with identity config   │
│  3. Set link phases from clock drift rates        │
└───────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────┐
│  Compute Observables via ChronoGaugeOps           │
│  → mass_density_action()                          │
│  → vorticity_tensor(), vorticity_z()              │
│  → virial_ratio(), tolman_temperature()           │
│  → fleet_coherence()                              │
└───────────────────────────────────────────────────┘
        │
        ▼
┌───────────────────────────────────────────────────┐
│  Aggregate & PCA Analysis                         │
│  → 7 observables per epoch                        │
│  → Principal component decomposition              │
│  → Variance explained, loadings                   │
└───────────────────────────────────────────────────┘
```

### 3.3 Lattice Construction Strategy

> **Design Decision**: Per-epoch lattice construction

Each epoch with ≥6 satellites creates a fresh `ChronoGauge` lattice:

```rust
fn create_epoch_lattice(satellites: &[SpaceTimeCoordinate<FloatType>]) -> ChronoGauge<FloatType> {
    // Fixed lattice dimensions (independent of satellite count)
    let shape = [8, 8, 8, 8];  // 8^4 = 4096 sites
    let periodic = [true, true, true, true];
    let lattice = Arc::new(Lattice::new(shape, periodic));

    // Initialize with identity, then embed satellite data
    let mut gauge = LatticeGaugeField::identity(lattice, FloatType::one());

    // Embed clock drift rates as U(1) phases on temporal links
    for sat in satellites {
        embed_satellite_data(&mut gauge, sat);
    }

    gauge
}
```

**Rationale**:

- Matches legacy behavior (independent epochs)
- Avoids cross-contamination between time slices
- Enables parallel epoch processing

### 3.4 Core Types

#### 3.4.1 Validation Result

```rust
/// Result of LGT validation for one dataset
#[derive(Debug, Clone, Default)]
pub struct GaugeValidationResult<R> {
    // Metadata
    pub dataset_name: String,
    pub num_epochs: usize,
    pub num_lattice_sites: usize,

    // Core observables (LGT equivalents)
    pub mean_mass_density: R,        // Wilson action
    pub mean_curl_magnitude: R,      // ||vorticity_tensor||_F
    pub mean_gradient_magnitude: R,  // U(1) phase gradient
    pub mean_vorticity_z: R,         // SU(2) z-component

    // Physics correlations
    pub time_velocity_correlation: R,
    pub virial_ratio: R,
    pub tolman_consistency: R,

    // Fleet coherence (new in LGT)
    pub fleet_coherence: (R, R),     // (mean, variance)
    pub topological_charge: R,
}
```

#### 3.4.2 PCA Input Matrix

```rust
// 7 observables × N epochs
let pca_matrix: Vec<Vec<R> > = vec![
    mass_density_vals,       // [0] Laplacian equivalent
    curl_magnitude_vals,     // [1] Frame-dragging
    gradient_magnitude_vals, // [2] Geoid
    vorticity_z_vals,        // [3] Spin
    time_velocity_vals,      // [4] Momentum
    virial_ratio_vals,       // [5] Energy
    tolman_vals,             // [6] Thermodynamics
];
```

---

## 4. Implementation

### 4.1 Module Structure

```text
chrono_dynamics/experiments/
├── bin/E02_chrono_gauge/
│   ├── main.rs               # Experiment entry point
│   └── README.md             # Experiment documentation
├── src/
│   ├── types/
│   │   ├── mod.rs            # Export new types
│   │   └── gauge_types.rs    # [NEW] GaugeValidationResult, EpochMetrics
│   └── utils/
│       ├── mod.rs            # Export new utils
│       ├── gauge_utils/      # [NEW] Gauge-specific utilities
│       │   ├── mod.rs
│       │   ├── lattice_utils.rs   # create_epoch_lattice(), embed_satellite_data()
│       │   └── observable_utils.rs # frobenius_norm(), compute_observables()
│       └── print_utils/
│           ├── mod.rs
│           └── gauge_print.rs # [NEW] Print utilities for gauge experiment
```

The experiment uses existing library code:

- `experiments::proces_utils::interpolate_space_time()` — Lagrange interpolation
- `experiments::statistics_utils::calculate_pca()` — PCA decomposition
- `experiments::print_utils` — Output formatting
- `deep_causality_physics::ChronoGaugeOps` — LGT observables

### 4.2 Main Entry Point

```rust
fn main() -> io::Result<()> {
    // 1. Setup
    print_header();
    let data_path = get_gnss_data_input_path();

    // 2. Process years
    let mut all_results: Vec<GaugeValidationResult<FloatType>> = Vec::new();

    for year in YEARS {
        let year_results = process_year(year, &data_path)?;
        all_results.extend(year_results);
    }

    // 3. Summary statistics
    print_summary(&all_results);

    // 4. PCA analysis
    run_pca_analysis(&all_results);

    Ok(())
}
```

### 4.3 Epoch Processing

```rust
fn process_epoch(satellites: &[SpaceTimeCoordinate<FloatType>])
                 -> Result<EpochMetrics<FloatType>, TopologyError>
{
    // Skip epochs with insufficient satellites
    if satellites.len() < MIN_SATELLITES {
        return Err(TopologyError::InsufficientData);
    }

    // Create per-epoch lattice
    let gauge = create_epoch_lattice(satellites);

    // Compute LGT observables
    let mass_density = gauge.mass_density_action()?;
    let vorticity = gauge.vorticity_tensor()?;
    let vorticity_z = gauge.vorticity_z()?;
    let virial = gauge.virial_ratio()?;
    let tolman = gauge.tolman_temperature()?;
    let (coherence_mean, coherence_var) = gauge.fleet_coherence()?;

    // Compute derived observables
    let curl_magnitude = frobenius_norm(&vorticity);

    Ok(EpochMetrics {
        mass_density,
        curl_magnitude,
        vorticity_z,
        virial,
        tolman,
        coherence_mean,
        coherence_var,
    })
}
```

---

## 5. Validation & Verification

### 5.1 Build Verification

```bash
cd /Users/marvin/RustroverProjects/dcl/chrono_dynamics
cargo build -p experiments --bin chrono_gauge
```

### 5.2 Unit Tests

Tests should verify:

- [ ] Lattice construction from satellite data
- [ ] Observable computation matches expected ranges
- [ ] PCA produces meaningful variance decomposition

### 5.3 Integration Test

```bash
cargo run -p experiments --bin chrono_gauge --release 2>&1 | tee E02_chrono_gauge.txt
```

### 5.4 Comparison with Legacy

The legacy output file:

```
gqcd/bin/gqcd_chrono_manifold/gqcd_chrono_manifold.txt
```

**Comparison criteria**:

| Metric                 | Legacy | Expected CGLT   | Notes                            |
|------------------------|--------|-----------------|----------------------------------|
| Valid epochs           | ~N     | ~N              | Should match within 10%          |
| PC1 variance explained | ~X%    | Similar pattern | Dominant physics mode preserved  |
| Observable signs       | ±      | Same            | Physics interpretation preserved |
| Observable magnitudes  | M      | Different       | LGT has different discretization |

> [!IMPORTANT]
> **Non-Goal**: Exact numerical matching with legacy DEC results.
>
> **Goal**: Improved accuracy and Monte Carlo error bars.

---

## 6. Implementation Roadmap

### Phase 1: Experiment Structure

- [ ] Create `main.rs` with full experiment logic
- [ ] Add `README.md` documentation
- [ ] Update `Cargo.toml` with binary target

### Phase 2: Core Implementation

- [ ] Port epoch processing loop
- [ ] Implement lattice construction from satellite data
- [ ] Wire up `ChronoGaugeOps` observable calls
- [ ] Port physics correlation computations

### Phase 3: Analysis

- [ ] Port PCA analysis with 7 metrics
- [ ] Implement summary statistics output
- [ ] Add Monte Carlo error estimation (optional)

### Phase 4: Verification

- [ ] Run on full 2016-2018 dataset
- [ ] Compare with legacy output
- [ ] Document findings in `walkthrough.md`

---

## 7. Success Criteria

### 7.1 Functional Requirements

| ID | Requirement                           | Verification              |
|----|---------------------------------------|---------------------------|
| F1 | Process all GNSS datasets (2016-2018) | No crashes, valid output  |
| F2 | Compute 7 LGT observables per epoch   | All values finite         |
| F3 | PCA produces 3+ principal components  | Cumulative variance > 70% |
| F4 | Use `DoubleFloat` precision           | Type checks pass          |
| F5 | Match legacy epoch counts (±10%)      | Comparison test           |

### 7.2 Non-Functional Requirements

| ID | Requirement | Criterion                         |
|----|-------------|-----------------------------------|
| N1 | Build time  | `cargo build` < 30s               |
| N2 | Runtime     | Full experiment < 10 minutes      |
| N3 | Memory      | < 4GB peak usage                  |
| N4 | Parallelism | Uses Rayon for dataset processing |

---

## 8. References

1. **CGLT Specification**: `chrono_dynamics/specs/current/physics_chrono_gauge.md`
2. **Legacy Experiment**: `gqcd/bin/gqcd_chrono_manifold/main.rs`
3. **ChronoGaugeOps Trait**: `deep_causality_physics/src/theories/chrono_dynamics/chrono_ops.rs`
4. **E01 Reference**: `chrono_dynamics/experiments/bin/E01_chrono_mass/main.rs`

---

## Appendix A: Variable Name Mapping

| PCA Variable Name           | Legacy Source          | LGT Source                |
|-----------------------------|------------------------|---------------------------|
| Laplacian (Mass/Tidal)      | `laplacian_chrono()`   | `mass_density_action()`   |
| Curl (Frame-Dragging)       | `curl()` magnitude     | `vorticity_tensor()` norm |
| Gradient (Geoid)            | `gradient()` magnitude | U(1) link phase gradient  |
| Vorticity-Z (Spin)          | `curl()[2]`            | `vorticity_z()`           |
| Time-Velocity Corr          | Pearson correlation    | Polyakov susceptibility   |
| Virial Ratio (Energy)       | \|K\| / \|P\|          | `virial_ratio()`          |
| Tolman Consistency (Thermo) | 1 - (σ/μ) of τ ratio   | `tolman_temperature()`    |

## Appendix B: Expected Output Format

```
╔══════════════════════════════════════════════════════════════════════════════╗
║               E02 CHRONO-GAUGE VALIDATION EXPERIMENT                         ║
║               Lattice Gauge Theory Cross-Validation                          ║
╚══════════════════════════════════════════════════════════════════════════════╝


📊 Data path: /Users/marvin/RustroverProjects/dcl/gqcd/data/gnss
Folder found at path: /Users/marvin/RustroverProjects/dcl/gqcd/data/gnss

...

╔══════════════════════════════════════════════════════════════════════════════╗
║                     CHRONO-MANIFOLD VALIDATION SUMMARY                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

📊 STATISTICS
  Total Epochs Analyzed: 461206
  Total Vertices:        1844824
  Mean |Laplacian|:      2.016304e-44
  Mean |Curl|:           0.000000e0
  Mean |Gradient|:       2.119068e-14
  Mean Vorticity_z:      0.000000e0
  Time-Velocity Corr:    0.999093
  Tolman Consistency:    1.000000

📐 VALIDATION TESTS
  Laplacian Trace ≈ 0:   100.0% pass (Poisson equation in vacuum)
  Curl ≈ 0:              100.0% pass (Scalar field property)

🔬 CROSS-VALIDATION WITH ORIGINAL EXPERIMENTS (7 TOTAL)
  ┌──────────────────────────────┬──────────────────┬───────────┬─────────────────┬──────────────────┬────────────────┐
  │ Experiment                   │ Metric           │ Reference │ Gauge Method    │ Manifold Method  │ Precision Gain │
  ├──────────────────────────────┼──────────────────┼───────────┼─────────────────┼──────────────────┼────────────────┤
  │ gqcd_chrono_mass             │ Earth Mass (kg)  │ 5.972e24  │                 │ 5.9721e24        │ 1×             │
  │ gqcd_chrono_tidal            │ Laplacian (∇²𝒯)  │ 0.00      │                 │ 2.02e-44         │           │
  │ gqcd_chrono_frame_dragging   │ Curl (∇×∇𝒯)      │ 0.00      │                 │ 0.00e0           │ ∞ (exact)      │
  │ gqcd_chrono_geoid            │ |∇𝒯| (field)     │ ~1e-10    │                 │ 2.12e-14         │ Exact          │
  │ gqcd_chrono_spin             │ Vorticity_z      │ 0.00      │                 │ 0.00e0           │ ∞ (exact)      │
  │ gqcd_chrono_momentum         │ R(τ̇, v²)         │ +1.00     │                 │ 0.9991           │ Exact          │
  │ gqcd_chrono_thermodynamics   │ Tolman T×τ       │ 1.00      │                 │ 1.0000           │ Exact          │
  └──────────────────────────────┴──────────────────┴───────────┴─────────────────┴──────────────────┴────────────────┘

════════════════════════════════════════════════════════════════════════════════
CONCLUSION
════════════════════════════════════════════════════════════════════════════════

✅ VALIDATION SUCCESSFUL (7 EXPERIMENTS)
   - Mass:           Earth mass derived from Laplacian (0.001% error)
   - Tidal:          Laplace equation Tr(H) ≈ 0 verified
   - Frame-Dragging: Nilpotency d∘d = 0 confirmed (curl = 0)
   - Geoid:          Gradient field encodes gravitational acceleration
   - Spin:           Vorticity z-component aligned with rotation axis
   - Momentum:       Time-velocity correlation confirmed
   - Thermodynamics: Tolman T×τ = constant verified

   Differential geometry confirms scalar field hypothesis.
   Physics is invariant across computational methods.

════════════════════════════════════════════════════════════════════════════════

╔══════════════════════════════════════════════════════════════════════════════╗
║               PRINCIPAL COMPONENT ANALYSIS                                   ║
║               Dimensionality Reduction of Validation Metrics                 ║
╚══════════════════════════════════════════════════════════════════════════════╝

📊 PCA Results (331 observations, 7 variables)
───────────────────────────────────────────────────────────────────────────────

Variance Explained by Principal Components:
  PC1: 100.00% (Cumulative: 100.00%)
  PC2: 0.00% (Cumulative: 100.00%)
  PC3: 0.00% (Cumulative: 100.00%)

───────────────────────────────────────────────────────────────────────────────
Principal Component Loadings (Contribution of each variable):

PC1 (explains 100.0% of variance):
  Laplacian (Mass/Tidal)             + 0.000 
  Curl (Frame-Dragging)              + 0.000 
  Gradient (Geoid)                   + 0.000 
  Vorticity-Z (Spin)                 + 0.000 
  Time-Velocity Corr (Momentum)      + 0.707 ████████████████████████████
  Virial Ratio (Energy)              + 0.707 ████████████████████████████
  Tolman Consistency (Thermo)        + 0.000 

PC2 (explains 0.0% of variance):
  Laplacian (Mass/Tidal)             + 0.378 ███████████████
  Curl (Frame-Dragging)              + 0.378 ███████████████
  Gradient (Geoid)                   + 0.378 ███████████████
  Vorticity-Z (Spin)                 + 0.378 ███████████████
  Time-Velocity Corr (Momentum)      + 0.378 ███████████████
  Virial Ratio (Energy)              + 0.378 ███████████████
  Tolman Consistency (Thermo)        + 0.378 ███████████████

PC3 (explains 0.0% of variance):
  Laplacian (Mass/Tidal)             + 0.378 ███████████████
  Curl (Frame-Dragging)              + 0.378 ███████████████
  Gradient (Geoid)                   + 0.378 ███████████████
  Vorticity-Z (Spin)                 + 0.378 ███████████████
  Time-Velocity Corr (Momentum)      + 0.378 ███████████████
  Virial Ratio (Energy)              + 0.378 ███████████████
  Tolman Consistency (Thermo)        + 0.378 ███████████████

───────────────────────────────────────────────────────────────────────────────
Interpretation Guide:
  • PC1 typically represents the dominant validation pattern
  • High loadings (>0.5) indicate strong contribution to that component
  • Components with >70% cumulative variance capture most information
  • Variables with similar loadings are correlated
═══════════════════════════════════════════════════════════════════════════════
```
