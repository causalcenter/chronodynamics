# Physics Specification: Chrono-Source Integration

**Status:** Draft
**Context:** Unifying `ChronoGauge` Physics with `LatticeGaugeField` Source Architecture.

## 1. Objective

Refactor `deep_causality_physics` to fully leverage the new `Source` ($S$) capability of `LatticeGaugeField`.
Eliminate the method signature `solve_j2(&self, data: &[C])` entirely in favor of `solve_j2(&self)`, where the data comes directly from `self.source()`.

## 2. Type Architecture

### 2.1 ChronoGauge Alias

The `ChronoGauge` alias will now default to carrying data, or be explicit about it.
To support the "Inverse Problem" (Metrology) directly:

```rust
// Proposed Alias in alias/mod.rs
pub type ChronoGauge<FloatType, FieldSource = Vec<SpaceTimeCoordinate<FloatType>>> =
    LatticeGaugeField<SU2_U1, DIM_4D_SPACE_TIME, Complex<FloatType>, FloatType, FieldSource>;
```

This means:
*   `ChronoGauge<f64>` now implies a data-carrying field (or at least one typed for it).
*   Users can still opt-out with `ChronoGauge<f64, ()>` for pure vacuum if needed, but the default favors the Chrono-Gauge use case (Data -> Geometry).

## 3. Operations Trait (`ChronoGaugeOps`)

### 3.1 Refactoring `solve_j2`

The `ChronoGaugeOps` trait must be updated to remove external data dependencies.

```rust
pub trait ChronoGaugeOps<R: RealField> {
    // ... other methods ...

    // OLD: fn solve_j2<C>(&self, data: &[C]) -> Result<R, TopologyError>;

    // NEW: Derives J2 directly from the internal source S
    fn solve_j2(&self) -> Result<R, TopologyError>;
}
```

### 3.2 The `ChronoObserver` Helper Trait

How do we implement `solve_j2` for generic `S`? We need a trait that extracts observations from `S`.

```rust
/// Helper trait for source data access.
/// This allows ChronoGaugeOps to be implemented generically while delegating to S.
pub trait ChronoObserver<R, C> {
    /// Returns the observation data if this source contains it.
    fn as_observations(&self) -> Option<&[C]>;
}

// Case 1: The Source IS the Data (Vec<C>)
impl<R, C> ChronoObserver<R, C> for Vec<C> {
    fn as_observations(&self) -> Option<&[C]> {
        Some(self)
    }
}

// Case 2: Vacuum Source (())
impl<R, C> ChronoObserver<R, C> for () {
    fn as_observations(&self) -> Option<&[C]> {
        None
    }
}
```

## 4. Implementation Strategy

### 4.1 `chrono_ops_impl.rs`

The blanket implementation unifies everything:

```rust
impl<R, S, C> ChronoGaugeOps<R> for LatticeGaugeField<SU2_U1, 4, Complex<R>, R, S>
where
    R: RealField + ...,
    S: ChronoObserver<R, C>,    // The Source must be an Observer
    C: SpaceTimeCoord<R>,       // The Data must be Coordinates
{
    fn solve_j2(&self) -> Result<R, TopologyError> {
        // 1. Get data from Source
        let data = self.source().as_observations()
            .ok_or_else(|| TopologyError::LatticeGaugeError("Field Source contains no observations for J2 calculation".into()))?;

        // 2. Perform J2 regression on `data`
        // ... implementation ...
    }
}
```

## 5. Migration Steps

1.  **Define `ChronoObserver`**: In `chrono_ops.rs`.
2.  **Update `ChronoGaugeOps`**: Remove `data` argument from `solve_j2`.
3.  **Update `chrono_ops_impl.rs`**: Implement `solve_j2` using `self.source().as_observations()`.
4.  **Refactor `E02 Experiment`**:
    *   Construct field with `with_source(all_coords)`.
    *   Call `field.solve_j2()`.
