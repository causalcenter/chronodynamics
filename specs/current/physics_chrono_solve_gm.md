# Implementation Plan: `solve_gm` via Gauge Field Observables

**Status:** Draft  
**Author:** AI Assistant  
**Date:** 2026-01-23

---

## 1. Objective

Implement `solve_gm(&self) -> Result<R, TopologyError>` that **calculates GM from the gauge field's link variables and lattice observables** - not by wrapping the analytical formula.

The gauge field's `source: Vec<SpaceTimeCoordinate<R>>` provides the clock/orbit data.  
This data is used to **populate the link variables $U_\mu(x)$**, which encode gravitational effects.  
Then GM is **extracted from lattice observables** like Polyakov loops or Wilson action.

---

## 2. Two-Phase Architecture

### Phase 1: Populate Links from Source Data

The source `SpaceTimeCoordinates` contain clock drift rates $\dot{\tau}(r)$.  
These must be encoded into the **temporal link variables** $U_0(x)$:

$$U_0(x) = \exp\left(i \cdot \frac{\dot{\tau}(r(x)) - 1}{\epsilon} \cdot a\right)$$

where:
- $\dot{\tau}(r)$ = clock drift rate at radius $r$ (from source data)
- $\epsilon$ = encoding scale factor (to be determined)
- $a$ = lattice spacing

### Phase 2: Extract GM from Lattice Observables

Once links are populated, use **Polyakov loop** or **Wilson action** to extract GM:

$$|P(r)| = e^{-GM/(rc^2)} \approx 1 - \frac{GM}{rc^2}$$

Therefore:
$$\boxed{GM = c^2 \cdot r \cdot (1 - |P(r)|)}$$

Or from two radii:
$$\boxed{GM = c^2 \cdot \frac{|P(r_b)| - |P(r_a)|}{1/r_a - 1/r_b}}$$

---

## 3. Link Population Strategy

### 3.1 Mapping Source Data to Lattice

Each `SpaceTimeCoordinate` has:
- `r_m`: radius (meters)
- `clock_drift_rate`: $\dot{\tau} = d\tau/dt$

The lattice has D dimensions with shape `[N_t, N_r, N_θ, N_φ]`.

**Mapping radius to lattice index:**
```rust
fn radius_to_lattice_index(&self, r: R) -> usize {
    let r_min = EARTH_RADIUS;          // ~6.37e6 m
    let r_max = GEO_ORBIT_RADIUS;       // ~4.2e7 m
    let shape = self.lattice.shape();
    
    let r_norm = (r - r_min) / (r_max - r_min);
    (r_norm * (shape[1] as f64)).floor() as usize
}
```

### 3.2 Encoding Clock Drift into Links

For each radial shell, compute average clock drift rate from source data, then set:

```rust
fn populate_links_from_source(&mut self) -> Result<(), TopologyError> {
    let data = self.source();
    let shape = self.lattice.shape();
    
    // Bin data by radial index
    let mut radial_bins: Vec<Vec<R>> = vec![Vec::new(); shape[1]];
    for coord in data {
        let r_idx = self.radius_to_lattice_index(coord.r_m);
        radial_bins[r_idx].push(coord.clock_drift_rate);
    }
    
    // For each temporal link at radial index, set phase from avg clock drift
    for cell in self.lattice.cells(1) {
        let pos = cell.position();
        let mu = cell.orientation().trailing_zeros() as usize;
        
        if mu == 0 {  // Temporal link
            let r_idx = pos[1];
            if !radial_bins[r_idx].is_empty() {
                let avg_drift = radial_bins[r_idx].iter().sum::<R>() 
                               / R::from(radial_bins[r_idx].len());
                
                // Encode: U_0 = exp(i * phase)
                // phase ∝ (1 - avg_drift) which is ∝ GM/rc²
                let phase = (R::one() - avg_drift) * ENCODING_SCALE;
                let link = LinkVariable::from_phase(phase)?;
                self.set_link(cell, link);
            }
        }
    }
    
    Ok(())
}
```

---

## 4. GM Extraction via Polyakov Loop

### 4.1 Physics

The Polyakov loop at spatial position $\vec{x}$ is:

$$P(\vec{x}) = \text{Tr}\left[\prod_{t=0}^{N_t-1} U_0(\vec{x}, t)\right]$$

For gravitational encoding:
$$|P(r)| \approx \exp\left(-\frac{GM}{rc^2} \cdot N_t \cdot a_t\right)$$

**Properties:**
- Linear dependence on GM → **high precision** for weak fields
- Sensitive to individual link errors along temporal path

### 4.2 Implementation with Spatial Averaging

Average over angular coordinates for spherically symmetric gravity:

```rust
use rayon::prelude::*;

fn solve_gm_polyakov(&self) -> Result<R, TopologyError> {
    let shape = self.lattice.shape();
    
    // Parallel over radial shells
    let measurements: Vec<(R, R)> = (0..shape[1])
        .into_par_iter()
        .map(|r_idx| {
            let r = self.lattice_index_to_radius(r_idx);
            
            // Inner loop: parallel over angular coordinates
            let angular_sum: R = (0..shape[2])
                .into_par_iter()
                .flat_map(|theta_idx| {
                    (0..shape[3]).into_par_iter().map(move |phi_idx| {
                        let site = [0, r_idx, theta_idx, phi_idx];
                        self.try_polyakov_loop(&site, 0)
                            .map(|p| p.abs())
                            .unwrap_or(R::zero())
                    })
                })
                .sum();
            
            let count = shape[2] * shape[3];
            let avg_p = angular_sum / R::from(count);
            (r, avg_p)
        })
        .collect();
    
    let c_sq = SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    let (slope, _) = linear_regression_inv_r(&measurements)?;
    
    Ok(-slope * c_sq)
}
```

---

## 5. GM Extraction via Wilson Action

### 5.1 Physics

The Wilson action density at radius $r$:

$$s(r) = \beta \left(1 - \frac{1}{N}\text{ReTr}(U_{01}(r))\right)$$

For gravity: $s(r) \propto |\nabla\Phi|^2 \propto (GM/r^2)^2$.

**Properties:**
- Quadratic dependence on GM → **less sensitive** in weak field
- Averages over all plaquettes → **robust** against noise

### 5.2 Implementation

```rust
fn solve_gm_from_action(&self) -> Result<R, TopologyError> {
    let shape = self.lattice.shape();
    
    // Parallel over radial shells
    let gm_estimates: Vec<R> = (0..shape[1])
        .into_par_iter()
        .filter_map(|r_idx| {
            let r = self.lattice_index_to_radius(r_idx);
            
            // Inner loop: parallel over angular coordinates
            let action_sum: R = (0..shape[2])
                .into_par_iter()
                .flat_map(|theta_idx| {
                    (0..shape[3]).into_par_iter().map(move |phi_idx| {
                        let site = [0, r_idx, theta_idx, phi_idx];
                        self.try_plaquette_action(&site, 0, 1)
                            .unwrap_or(R::zero())
                    })
                })
                .sum();
            
            let count = shape[2] * shape[3];
            let avg_action = action_sum / R::from(count);
            
            if avg_action > R::zero() {
                Some(avg_action.sqrt() * r * r)
            } else {
                None
            }
        })
        .collect();
    
    Ok(gm_estimates.par_iter().sum::<R>() / R::from(gm_estimates.len()))
}
```

---

## 6. Hybrid Method (Recommended)

Combine both methods for **precision + robustness**:

```rust
fn solve_gm(&self) -> Result<R, TopologyError> {
    let gm_polyakov = self.solve_gm_polyakov()?;
    let gm_action = self.solve_gm_from_action()?;
    
    // Consistency check
    let tolerance = R::from(0.05);  // 5%
    let relative_diff = ((gm_polyakov - gm_action) / gm_polyakov).abs();
    
    if relative_diff < tolerance {
        // Methods agree: weighted average (favor precision)
        let w_p = R::from(0.7);
        let w_a = R::from(0.3);
        Ok(w_p * gm_polyakov + w_a * gm_action)
    } else {
        // Disagreement: trust Polyakov (more precise for weak field)
        // but could log warning here
        Ok(gm_polyakov)
    }
}
```

---

## 7. Encoding Scale (Resolved)

The clock drift rate $\dot{\tau} \approx 1 - \frac{GM}{rc^2}$ is O($10^{-10}$) for GNSS.

**Encoding:**
```rust
// phase = (1 - clock_drift_rate) * N_t
// This ensures Polyakov loop |P| = exp(-phase * N_t) = exp(-GM/(rc²) * N_t)
let phase = (R::one() - avg_drift) * R::from(shape[0]);  // N_t = temporal extent
```

The factor of $N_t$ ensures the Polyakov loop (product of $N_t$ links) correctly encodes the integrated potential.

---

## 8. Complete Workflow

```rust
// 1. Create field and attach source data
let mut gauge_field = ChronoGauge::identity(lattice, 1.0)
    .with_source(interpolate_space_time(&clocks, &orbits));

// 2. Populate link variables from source (clock drift → link phase)
// This is an EXPLICIT call (not automatic in with_source)
gauge_field.populate_links_from_source()?;

// 3. Extract GM using hybrid method
let gm = gauge_field.solve_gm()?;  // Combines Polyakov + Action
```

---

## 9. Implementation Steps

1. **Add helper methods** to `LatticeGaugeField`:
   - `radius_to_lattice_index(r) -> usize`
   - `lattice_index_to_radius(idx) -> R`

2. **Add link population method**:
   - `populate_links_from_source(&mut self) -> Result<(), TopologyError>`
   - Encoding: `phase = (1 - drift) * N_t`

3. **Implement GM extraction methods**:
   - `solve_gm_polyakov(&self)` - with spatial averaging
   - `solve_gm_from_action(&self)` - with spatial averaging
   - `solve_gm(&self)` - hybrid combining both

4. **Update E01 experiment** to use new workflow

---

## 10. Design Decisions Summary

| Question | Decision | Rationale |
|----------|----------|-----------|
| **Encoding scale** | `(1 - drift) * N_t` | Polyakov loop product naturally accumulates $N_t$ factors |
| **Spatial averaging** | Yes, average over θ,φ | Spherical symmetry → reduces noise |
| **Link initialization** | Explicit `populate_links_from_source()` | Separation of concerns; user can inspect source first |
| **Hybrid method** | 70% Polyakov + 30% Action (if agree) | Precision + robustness |