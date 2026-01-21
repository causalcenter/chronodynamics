# J2 Oblateness Calculation via Chrono Gauge Field

> [!IMPORTANT]
> This plan rewrites `solve_j2` to use the **Lattice Gauge Field** infrastructure rather than Newtonian scalar
> regression.

## Problem Statement

The current `solve_j2` implementation
in [chrono_hkt.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/deep_causality_physics/src/theories/chrono_dynamics/chrono_hkt.rs#L355-L446)
has a critical flaw:

```rust
pub fn solve_j2<C: SpaceTimeCoord<T>>(_field: &ChronoGauge<T>, data: &[C]) -> T
```

1. The `_field` parameter is **completely unused** (note the underscore prefix)
2. The implementation performs **scalar regression** on raw clock drift rates
3. Uses Newtonian potential formula: `monopole_rate = -GM / (c² r)`
4. Result: **82% error** because J2 is contaminated by systematic noise at the `~10⁻⁷` level

### Why This Approach Fails

| Approach         | What It Measures             | Noise Handling                            |
|------------------|------------------------------|-------------------------------------------|
| Current (Scalar) | Clock drift rates directly   | Noise averages → biased mean              |
| Gauge Field      | Wilson Loop / Plaquette flux | Loop integral cancels gauge-variant noise |

The J2 signal is `~10⁻⁷`, same magnitude as orbital/clock noise. Scalar regression cannot separate them.

---

## Proposed Solution: Gauge-Based J2 Inversion

### Physics Foundation

In Gauge Theory (U(1) × SU(2)), J2 is **not** a scalar bump—it's a **non-zero Field Strength Tensor** component
representing a quadrupole moment in the electric component F₀ᵢ.

**Key Insight**: The Wilson Loop around a plaquette has a non-zero flux that encodes the physical curvature. Systematic
noise integrates to zero around closed loops (gauge invariance), but J2's physical contribution does not.

### Mathematical Framework

1. **Link Variables** $U_{AB}$ between satellites A and B:
   Encode the potential difference (phase) between points:
   $$ U_{AB} = \exp(i \phi_{AB}) \quad \text{where} \quad \phi_{AB} \approx \frac{c^2}{\beta} \int_A^B A_\mu dx^\mu $$
   Physically, this relates to clock drift difference:
   $$ \phi_{AB} \propto (\dot{\tau}_B - \dot{\tau}_A) $$

2. **Wilson Loop** around a satellite triangle A→B→C→A:
   $$ W = U_{AB} \cdot U_{BC} \cdot U_{CA} = \exp(i \Phi_{\text{total}}) $$
   $$ \Phi_{\text{total}} = \phi_{AB} + \phi_{BC} + \phi_{CA} $$

3. **Gauge Invariance**:
   Systematic errors (gauge noise) behave like a derivative $\partial_\mu \Lambda$.
   $$ \oint \nabla \Lambda \cdot dl = 0 $$
   Closed loops naturally filter out this "gradient noise".

4. **J2 Extraction**:
   The surviving flux $\Phi_{\text{total}}$ is proportional to the curvature components $R_{0r0r}$ (electric-like).
   For an oblate spheroid (Earth), this curvature varies with latitude $\phi$:
   $$ F_{0r} \propto J_2 \cdot P_2(\sin \phi) = J_2 \cdot \frac{1}{2}(3\sin^2\phi - 1) $$

---

## Proposed Changes

---

### [MODIFY] [chrono_hkt.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/deep_causality_physics/src/theories/chrono_dynamics/chrono_hkt.rs)

#### Step 1: Add Trait Methods to SpaceTimeCoord (lines 262-274)

```rust
pub trait SpaceTimeCoord<T: RealField> {
    // ... existing methods
    fn satellite_id(&self) -> u32;
    fn latitude(&self) -> T;
    fn longitude(&self) -> T;
    fn timestamp(&self) -> i64;
}
```

#### Step 2: Add SatelliteTriangle Helper Struct

```rust
pub struct SatelliteTriangle<'a, T: RealField, C: SpaceTimeCoord<T>> {
    pub coords: [&'a C; 3],
    _marker: PhantomData<T>,
}

impl<'a, T: RealField + Clone + From<f64>, C: SpaceTimeCoord<T>> SatelliteTriangle<'a, T, C> {
    pub fn average_latitude(&self) -> T { /* ... */ }
    pub fn legendre_p2(&self) -> T {
        // P2(x) = (3x² - 1)/2
        let avg_lat = self.average_latitude();
        let sin_lat = avg_lat.sin();
        let half = T::from(0.5);
        half * (T::from(3.0) * sin_lat * sin_lat - T::one())
    }
}
```

#### Step 3: Rewrite solve_j2 (lines 355-446)

> [!WARNING]
> This is a BREAKING change. The function signature changes from `-> T` to `-> Result<T, TopologyError>`.

```rust
/// Solves for J2 oblateness coefficient using Gauge Field Wilson Loop analysis.
///
/// # Physics Formula
///
/// 1. **Link Phase**: $\phi_{ij} = \beta (\dot{\tau}_j - \dot{\tau}_i)$
/// 2. **Wilson Loop**: $W = \sum_{\Delta} \phi_{ij}$ (Triangle Flux)
/// 3. **Curvature**: $F_{0r} \approx W / \text{Area}$
/// 4. **J2 Regression**: $F_{0r} = k \cdot J_2 \cdot P_2(\sin \theta)$
pub fn solve_j2<C: SpaceTimeCoord<T>>(
    field: &ChronoGauge<T>,
    data: &[C],
) -> Result<T, TopologyError>
where
    T: RealField + Clone + From<f64> + Into<f64>,
{
    // =========================================================
    // 1. DATA GROUPING
    // =========================================================
    // Group satellites by epoch (timestamp) to ensure synchronous triangles.
    // Map: Timestamp -> Vec<&Satellite>
    let mut epoch_groups: std::collections::HashMap<i64, Vec<&C>> = std::collections::HashMap::new();
    for coord in data {
        epoch_groups.entry(coord.timestamp()).or_default().push(coord);
    }

    let gm = T::from(NEWTONIAN_CONSTANT_OF_GRAVITATION) * T::from(EARTH_MASS_KG);
    let c_sq = T::from(SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    let r_earth = T::from(6_378_137.0);

    let mut flux_measurements: Vec<(T, T)> = Vec::new(); // (P2, Flux)

    // =========================================================
    // 2. TRIANGLE FORMATION & WILSON LOOP COMPUTATION
    // =========================================================
    for (_ts, satellites) in epoch_groups.iter() {
        if satellites.len() < 3 { continue; }

        // Iterate all unique triangles (k-combinations)
        for i in 0..satellites.len() {
            for j in (i + 1)..satellites.len() {
                for k in (j + 1)..satellites.len() {
                    let triangle = SatelliteTriangle {
                        coords: [satellites[i], satellites[j], satellites[k]],
                        _marker: PhantomData,
                    };

                    // CALCULATE FLUX (WILSON LOOP)
                    // W = exp(i * sum(phase))
                    // Phase difference driven by potential difference (clock rate)
                    let flux = compute_triangle_wilson_flux(
                        &triangle,
                        field.beta(),
                        c_sq,
                        gm,
                    );

                    let p2 = triangle.legendre_p2();
                    flux_measurements.push((p2, flux));
                }
            }
        }
    }

    // =========================================================
    // 3. J2 EXTRACTION VIA REGRESSION
    // =========================================================
    // Model: Flux = Slope * P2(sin lat)
    // J2 is effectively the slope (normalized by geometry constants)

    // Linear regression on (x=P2, y=Flux)
    // ... (compute slope) ...

    // Normalize slope to J2
    // J2 = slope * (c² r³) / (GM Re²)
    // ...

    Ok(j2_value)
}

/// Computes the Wilson Loop Flux for a single satellite triangle.
fn compute_triangle_wilson_flux<T, C>(
    triangle: &SatelliteTriangle<T, C>,
    beta: &T,
    c_sq: T,
    gm: T,
) -> T
where
    T: RealField + Clone + From<f64>
{
    // 1. Get Clock Drift Rates ($\dot{\tau}$)
    let rate_a = triangle.coords[0].clock_drift_rate();
    let rate_b = triangle.coords[1].clock_drift_rate();
    let rate_c = triangle.coords[2].clock_drift_rate();

    // 2. Remove Monopole (Schwarzschild) Term
    // We only want the Multipole moments (J2, etc)
    // $\delta = \dot{\tau}_{measured} - (-GM/c^2 r)$
    let mono_a = -gm / (c_sq * triangle.coords[0].radius_m());
    let res_a = rate_a - mono_a;
    // ... repeat for b, c

    // 3. Compute Link Phases
    // $\phi_{AB} = \text{Residual}_B - \text{Residual}_A$
    // This represents the gauge potential difference along the path
    let phi_ab = res_b - res_a;
    let phi_bc = res_c - res_b;
    let phi_ca = res_a - res_c;

    // 4. Sum around Loop (Wilson Loop Phase)
    // $\Phi = \oint A \cdot dl = \phi_{AB} + \phi_{BC} + \phi_{CA}$
    let loop_sum = phi_ab + phi_bc + phi_ca;

    // Scale by Beta (Gauge Coupling)
    *beta * loop_sum
}
```

---

### [MODIFY] SpaceTimeCoordinate Implementation

Add new trait methods to `SpaceTimeCoordinate<R>` in `chrono_experiments`:

```rust
impl<R: RealField> SpaceTimeCoord<R> for SpaceTimeCoordinate<R> {
    fn satellite_id(&self) -> u32 { self.sat_id }
    fn latitude(&self) -> R {
        // Geocentric latitude: arcsin(z / r)
        let r = self.radius_m();
        if r.abs() < R::from(1.0) { R::zero() } else { (self.z_m() / r).asin() }
    }
    fn longitude(&self) -> R { self.y_m().atan2(self.x_m()) }
    fn timestamp(&self) -> i64 { self.timestamp as i64 }
}
```

---

### [MODIFY] [E02_chrono_gauge/main.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/chrono_experiments/bin/E02_chrono_gauge/main.rs#L107)

Update call site to handle `Result` and new signature:

```rust
// Previous:
// let derived_j2 = ChronoGaugeWitness::solve_j2(&gauge_field, &all_coords);

// New:
match ChronoGaugeWitness::solve_j2( & gauge_field, & all_coords) {
Ok(derived_j2) => {
// ... print success table ...
}
Err(e) => {
eprintln ! ("J2 Calculation Failed: {}", e);
}
}
```

---

## Verification Plan

### Unit Tests (new file: tests/chrono_hkt_tests.rs)

1. **test_satellite_triangle_average_latitude** - Verify P2 term calculation
2. **test_wilson_flux_closed_loop** - Identical satellites → flux ≈ 0

### Integration Test

```bash
cargo run --release --bin E02_chrono_gauge 2>&1 | grep -A10 "J2 ESTIMATION"
```

**Expected**: Error drops from 82% to < 20%
