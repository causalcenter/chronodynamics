# Observatory Experiment Comparison Report

> **Current**: `chrono_dynamics/experiments/bin/E01_chrono_observatory/E01_observatory.txt` (ChronoGauge-based)
> **Reference**: `gqcd/bin/gqcf_chrono_gravitational_observatory/gqcf_chrono_gravitational_observatory.txt` (Legacy)

---

## Executive Summary

The current ChronoGauge-based implementation exhibits **critical discrepancies** in all three analysis categories across
all years (2016, 2017, 2018). The most severe issues are:

1. **Data counts differ** – Current loads fewer epochs for 2016/2017, but matches closely for 2018
2. **Analysis results show zeroed values** – Standard deviations, correlations, and GR coefficients are all 0.0 in
   Current
3. **Anomaly detection is broken** – PC2 values are incorrectly 0.0%, producing invalid classifications

---

## 1. Data Loading / Processing Comparison

| Year     | Metric                 | Reference     | Current | Match?       |
|----------|------------------------|---------------|---------|--------------|
| **2016** | Total Epochs Collected | 32256         | 20825   | ❌ **-35.4%** |
|          | Valid 5T Measurements  | 29603 / 29350 | 20825   | ❌ **-29.7%** |
|          | J2 Points Analyzed     | 145221        | 104125  | ❌ **-28.3%** |
| **2017** | Total Epochs Collected | 32256         | 30874   | ❌ **-4.3%**  |
|          | Valid 5T Measurements  | 32104 / 32010 | 30874   | ❌ **-3.8%**  |
|          | J2 Points Analyzed     | 160758        | 154370  | ❌ **-4.0%**  |
| **2018** | Total Epochs Collected | 30695         | 30694   | ✅ ≈Match     |
|          | Valid 5T Measurements  | 30695         | 30694   | ✅ ≈Match     |
|          | J2 Points Analyzed     | 153444        | 153470  | ✅ ≈Match     |

### Root Cause Analysis (Source Code Identified)

> [!CAUTION]
> **2016 data has ~35% fewer epochs** – Root causes pinpointed in `data_manager/src/data_loader/load_fleet.rs`

---

#### Issue A: Strict Timestamp Matching Drops Most Clock Data

**File**: [load_fleet.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/data_manager/src/data_loader/load_fleet.rs#L148-L151)  
**Lines**: 148-151

```rust
let orbit = match orbit_map.get(&ts) {
    Some(o) => *o,
    None => continue,  // ← DROPS clock epoch if no exact orbit timestamp match
};
```

**Root Cause**: Clock data is sampled every **30 seconds**, orbit data every **15 minutes** (900s). The code requires an **exact timestamp match**:
- Per 15-minute orbit interval: 30 clock samples
- Only **1 clock sample** matches the orbit timestamp exactly
- **29 out of 30 clock samples are discarded** (~96.7% data loss per interval)

**Evidence**: Comment at lines 135-138 acknowledges this:
```rust
// Note: Orbits are typically every 15 mins, clocks every 30s.
// For exact matching, we need exact timestamps or interpolation.
// ...we just match strictly.
```

---

#### Issue B: Intersection Filtering Removes Partial Epochs

**File**: [load_fleet.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/data_manager/src/data_loader/load_fleet.rs#L86-L91)  
**Lines**: 86-91

```rust
let aligned_epochs: Vec<ObservatoryEpoch> = epochs
    .into_iter()
    .filter(|e| {
        satellite_ids.iter().all(|id| e.satellites.contains_key(*id))
    })
    .collect();
```

**Root Cause**: The code requires **ALL 7 satellites** to be present in every epoch. Any epoch missing even one satellite is discarded entirely.
- If E14 or E18 (eccentric) have gaps, entire epochs are dropped
- 2016 may have more gaps in early constellation data

---

#### Issue C: Why 2016 Loses More Data Than 2017/2018

| Year | Epochs (Reference) | Epochs (Current) | Loss |
|------|-------------------|------------------|------|
| 2016 | 32256 | 20825 | **-35.4%** |
| 2017 | 32256 | 30874 | -4.3% |
| 2018 | 30695 | 30694 | ≈0% |

**Hypothesis**: The Galileo constellation was less complete in 2016:
- Fewer operational satellites → more gaps per epoch
- Intersection filter discards more epochs
- 2017/2018 have better satellite availability → less data loss

---

#### Solution

1. **Implement orbit interpolation** in `merge_satellite_data()` to use nearest or linearly interpolated orbit for clock timestamps
2. **Relax intersection filter** to require minimum N satellites (e.g., 5) instead of ALL
3. **Add logging** to report data loss at each stage

---

## 2. Analysis Results Comparison

### 2.1 Clock Differential Analysis

| Year     | Metric          | Reference    | Current         | Match?       |
|----------|-----------------|--------------|-----------------|--------------|
| **2016** | Mean Clock Diff | +5.8644e6 ns | 5.2883e6 ns     | ❌ **-9.8%**  |
|          | Std Deviation   | 1.4623e6 ns  | **0.0000e0 ns** | ❌ **ZEROED** |
| **2017** | Mean Clock Diff | +4.6115e6 ns | 4.5825e6 ns     | ✅ ~Match     |
|          | Std Deviation   | 8.4918e4 ns  | **0.0000e0 ns** | ❌ **ZEROED** |
| **2018** | Mean Clock Diff | +1.7967e6 ns | 1.7966e6 ns     | ✅ ≈Match     |
|          | Std Deviation   | 6.7280e1 ns  | **0.0000e0 ns** | ❌ **ZEROED** |

> [!IMPORTANT]
> **Standard deviation is always 0.0 in Current** – This indicates the std dev calculation is broken or not being
> computed at all.

### 2.2 Correlation Analysis

| Year     | Metric                | Reference | Current    | Match?       |
|----------|-----------------------|-----------|------------|--------------|
| **2016** | E12-E18 Correlation   | 0.6002    | **0.0000** | ❌ **ZEROED** |
| **2017** | E12-E18 Correlation   | -0.9963   | **0.0000** | ❌ **ZEROED** |
| **2018** | E12-E18 Correlation   | -1.0000   | **0.0000** | ❌ **ZEROED** |
| All      | E14 Temporal Autocorr | 1.0000    | 1.0000     | ✅ Match      |

### 2.3 GR Verification

| Year     | Metric              | Reference        | Current            | Match?       |
|----------|---------------------|------------------|--------------------|--------------|
| **2016** | Observed Shift      | 6.6603e0 ns/s    | **0.0000e0 ns/s**  | ❌ **ZEROED** |
|          | Clock-Altitude Corr | -0.0011          | **0.0000**         | ❌ **ZEROED** |
|          | Measured GR Coeff   | -8.6931e-8 ns/km | **0.0000e0 ns/km** | ❌ **ZEROED** |
|          | Detrended Residual  | 2.5545e-1 ns     | **0.0000e0 ns**    | ❌ **ZEROED** |
| **2017** | Observed Shift      | 4.8022e0 ns/s    | **0.0000e0 ns/s**  | ❌ **ZEROED** |
|          | Clock-Altitude Corr | 0.0048           | **0.0000**         | ❌ **ZEROED** |
| **2018** | Observed Shift      | 1.9511e0 ns/s    | **0.0000e0 ns/s**  | ❌ **ZEROED** |
|          | Clock-Altitude Corr | 0.0145           | **0.0000**         | ❌ **ZEROED** |

> [!CAUTION]
> **All computed metrics are 0.0** – The ChronoGauge implementation is not computing any correlation, GR shift, or
> detrended residual values.

### 2.4 Swarm Covariance Analysis

| Year     | Metric     | Reference  | Current   | Notes              |
|----------|------------|------------|-----------|--------------------|
| **2016** | Quad Power | -3.7526e10 | 2.8278e12 | ❌ Sign + Magnitude |
|          | Mono Power | 6.0490e11  | 1.1903e13 | ❌ ~20x higher      |
| **2017** | Quad Power | 1.1716e8   | 1.6746e10 | ❌ ~143x higher     |
|          | Mono Power | 5.9030e9   | 6.0102e11 | ❌ ~102x higher     |
| **2018** | Quad Power | 7.3476e11  | 1.1938e11 | ✅ Same order       |
|          | Mono Power | -3.1032e9  | 4.3905e12 | ❌ Sign + Magnitude |

### 2.5 J2 Oblateness Inversion

| Year     | Metric      | Reference | Current      | Match?       |
|----------|-------------|-----------|--------------|--------------|
| **2016** | Inverted J2 | 3.3913e-3 | **0.0000e0** | ❌ **ZEROED** |
| **2017** | Inverted J2 | 6.4533e-3 | **0.0000e0** | ❌ **ZEROED** |
| **2018** | Inverted J2 | 5.0235e-2 | **0.0000e0** | ❌ **ZEROED** |

### Root Cause Analysis

> [!CAUTION]
> **Pattern: All derived/computed metrics are zero** – This strongly suggests the ChronoGauge implementation is either:
> 1. Not performing the actual calculations (placeholders?)
> 2. Using incorrect data types that lose precision
> 3. Missing the calculation pipeline entirely

---

## 3. Anomaly Detection Comparison

### 3.1 PCA Analysis

| Year     | PC  | Reference Variance | Current Variance | Reference Dominant | Current Dominant |
|----------|-----|--------------------|------------------|--------------------|------------------|
| **2016** | PC1 | 81.15%             | 80.80%           | E09 (+0.414)       | E14 (+0.556)     |
|          | PC2 | 12.46%             | 18.79%           | E12 (+0.970)       | E14 (+0.754)     |
|          | PC3 | 6.35%              | 0.41%            | E11 (+0.776)       | E12 (+0.669)     |
| **2017** | PC1 | 61.28%             | 97.29%           | E09 (+0.468)       | E12 (+0.940)     |
|          | PC2 | 28.17%             | 2.59%            | E19 (+0.634)       | E14 (+0.758)     |
| **2018** | PC1 | 66.85%             | 97.35%           | E09 (-0.461)       | E11 (+0.935)     |
|          | PC2 | 29.71%             | 2.65%            | E18 (+0.623)       | E12 (+0.874)     |

> [!WARNING]
> **PC variance distributions are vastly different** – Reference shows distributed variance (60-80% PC1), Current shows
> concentrated variance (97%+ PC1 for 2017/2018). This indicates the PCA input data differs significantly.

### 3.2 Detected Anomalies

| Year     | Tier               | Reference | Current   | Match?           |
|----------|--------------------|-----------|-----------|------------------|
| **2016** | Extreme (>99.9%)   | 8 events  | 52 events | ❌ **+550%**      |
|          | Significant (>99%) | 0 events  | 0 events  | ✅                |
|          | Elevated (>90%)    | 5 events  | 9 events  | ❌                |
| **2017** | Extreme            | 1 event   | 0 events  | ❌ **MISSED**     |
|          | Significant        | 2 events  | 0 events  | ❌ **MISSED**     |
|          | Elevated           | 23 events | 0 events  | ❌ **ALL MISSED** |
| **2018** | Extreme            | 10 events | 11 events | ✅ ~Match         |
|          | Significant        | 5 events  | 3 events  | ❌                |
|          | Elevated           | 21 events | 19 events | ✅ ~Match         |

### 3.3 Anomaly Details

**2016 Reference Top Anomalies (correctly identified satellites and PC2 values):**

- 2016-11-30 11:05:00 → E11, PC2=27.2%
- 2016-11-29 11:05:00 → E19, PC2=22.9%
- 2016-09-18 10:40:00 → E11, PC2=21.9%

**2016 Current Anomalies (broken):**

- All show `E??` (satellite not identified)
- All show `PC2=0.0%` (not computed)

> [!CAUTION]
> **Critical Bug: Current implementation shows `E??` and `PC2=0.0%` for all anomalies** – The satellite identification
> and PC2 value extraction are completely broken.

---

## 4. Root Cause Summary (Source Code Analysis)

> [!NOTE]
> Root causes pinpointed from source code in `chrono_dynamics/experiments/src/utils/observatory/` and `main.rs`.

---

### Issue 1: Hardcoded `E??` Satellite Identifier

- **Symptom**: All anomalies display `E??` instead of actual satellite ID
- **File**: [main.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/bin/E01_chrono_observatory/main.rs#L338-L345)
- **Line**: 342

```rust
println!(
    "║   {:<3} | {} | {:<4} | {:>5.1}% | {:<25} ║",
    i + 1,
    dt.format("%Y-%m-%d %H:%M:%S"),
    "E??",  // ← HARDCODED - should use actual satellite ID
    pc2_val * 100.0,
    tier_str
);
```

- **Root Cause**: The satellite identifier is hardcoded as the literal string `"E??"` instead of being looked up from the anomaly data
- **Solution**: Extract the dominant satellite from `AnomalyWindow` or identify from PCA loadings

---

### Issue 2: IsolatedGrResult Fields Hardcoded to Zero

- **Symptom**: `clock_alt_correlation`, `detrended_residual` are always 0.0
- **File**: [analysis.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/observatory/analysis.rs#L189-L196)
- **Lines**: 190, 194

```rust
IsolatedGrResult {
    clock_alt_correlation: 0.0, // TODO: Implement correlation  ← NOT IMPLEMENTED
    measured_coeff: avg_observed_shift / (mean_alt / 1000.0).abs().max(1.0),
    predicted_coeff: 6.6e-3 / 1300.0,
    coeff_agreement_pct: agreement.max(0.0),
    detrended_residual: 0.0,  // ← NOT IMPLEMENTED
},
```

- **Root Cause**: The `clock_alt_correlation` and `detrended_residual` fields have placeholder `0.0` values with a `TODO` comment – the actual calculation was never implemented
- **Solution**: Implement clock-altitude correlation using existing `calculate_correlation` utility; implement detrending

---

### Issue 3: E12-E18 Correlation Returns Zero Due to Data Misalignment

- **Symptom**: E12-E18 correlation is always `0.0000`
- **File**: [analysis.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/observatory/analysis.rs#L280-L316)
- **Lines**: 288-293

```rust
for epoch in epochs {
    if let (Some(s12), Some(s18)) = (epoch.satellites.get("E12"), epoch.satellites.get("E18")) {
        x.push(s12.clock_bias_ns);
        y.push(s18.clock_bias_ns);
    }
}
```

- **Root Cause**: The correlation calculation at lines 296-316 is correct, but the `if let` pattern silently skips epochs where E12 or E18 data is missing. If **all** epochs lack one satellite, `x` and `y` remain empty → returns 0.0
- **Investigation**: Check if `ObservatoryEpoch.satellites` HashMap is correctly populated with E12 and E18 keys
- **Solution**: Add logging to verify satellite data population; ensure data loader includes E12/E18 data

---

### Issue 4: J2 Inversion Returns Zero from ChronoGaugeWitness

- **Symptom**: `Inverted J2 = 0.0000e0` for all years
- **File**: [analysis.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/observatory/analysis.rs#L268-L277)
- **Line**: 269

```rust
let j2_val = ChronoGaugeWitness::solve_j2(field, &data_points);
let j2_f64: f64 = j2_val.into();

J2Result {
    satellites_used: 5,
    points_analyzed: count,
    inverted_j2: j2_f64,  // ← Always 0.0
    error_pct: ((j2_f64 - 1.08e-3).abs() / 1.08e-3) * 100.0,
}
```

- **Root Cause**: `ChronoGaugeWitness::solve_j2()` returns zero. This could be due to:
  1. The gauge field `field` being an identity field (line 58-64 creates it with empty `HashMap`)
  2. The `solve_j2` algorithm not finding a solution
  3. Input `data_points` having incorrect units or values
- **Investigation**: Debug `ChronoGaugeWitness::solve_j2` implementation in `deep_causality_physics`
- **Solution**: Verify gauge field construction; check `solve_j2` algorithm for edge cases

---

### Issue 5: Standard Deviation Calculation Works but Data May Be Constant

- **Symptom**: Std deviation is 0.0 for clock differentials
- **File**: [analysis.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/observatory/analysis.rs#L175-L181)
- **Lines**: 175-181

```rust
let std_diff = if !diffs.is_empty() {
    let mean = diffs.iter().sum::<f64>() / diffs.len() as f64;
    let variance = diffs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / diffs.len() as f64;
    variance.sqrt()
} else {
    0.0
};
```

- **Root Cause**: The std dev calculation logic is correct. However, if all values in `diffs` are identical (or `diffs` is empty), variance = 0 → std\_dev = 0.
- **Hypothesis**: The data loader may be returning constant/identical clock bias values
- **Investigation**: Check `clock_bias_ns` values in loaded `ObservatoryEpoch` data
- **Solution**: Debug data loader to verify clock bias values have variance

---

### Issue 6: PCA Variance Distribution Differs (97% vs 60-80% in PC1)

- **Symptom**: Current shows 97%+ variance in PC1 vs Reference's 60-80%
- **File**: [coherence.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/observatory/coherence.rs#L22-L28)
- **Call**: [statistics_utils/mod.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/src/utils/statistics_utils/mod.rs#L239-L372)

```rust
// coherence.rs:28
let pca = calculate_pca(matrix, 5);
```

- **Root Cause**: The PCA algorithm in `calculate_pca` is correct. The difference is in the **input matrix data**:
  - If satellite data lacks variance (e.g., constant values), PC1 captures nearly all variance
  - If fewer satellites are present, dimensionality is reduced
- **Hypothesis**: Input data to PCA has less variance than legacy implementation
- **Solution**: Debug `fleet_data.matrix` to verify it contains varying clock values across satellites

---

### Issue 7: Anomaly PC2 Values Are Near-Zero

- **Symptom**: PC2 values printed as `0.0%` 
- **File**: [main.rs](file:///Users/marvin/RustroverProjects/dcl/chrono_dynamics/experiments/bin/E01_chrono_observatory/main.rs#L326-L346)

```rust
let pc2_val: f64 = anomaly.pc2_variance.clone().into();
...
println!(
    "║   {:<3} | {} | {:<4} | {:>5.1}% | {:<25} ║",
    i + 1,
    dt.format("%Y-%m-%d %H:%M:%S"),
    "E??",
    pc2_val * 100.0,  // ← Displays as 0.0% if pc2_variance is ~0
    tier_str
);
```

- **Root Cause**: This is a downstream effect of Issue 6. If PC1 captures 97%+ variance, PC2 has only ~2.5% → window-level PCA may compute PC2 even lower, approaching 0.0
- **Solution**: Fix the input data variance issue (Issue 5/6) to restore proper PCA distribution

---

## 5. Recommended Solutions

| Priority | Issue | File | Fix |
|----------|-------|------|-----|
| **P0** | Hardcoded `E??` | `main.rs:342` | Replace `"E??"` with satellite lookup from PCA loadings |
| **P0** | Strict timestamp matching | `load_fleet.rs:148-151` | Implement orbit interpolation (nearest/linear) |
| **P0** | IsolatedGR zeros | `analysis.rs:190-194` | Implement `clock_alt_correlation` and `detrended_residual` |
| **P1** | Intersection filter too strict | `load_fleet.rs:86-91` | Require minimum 5 sats instead of ALL 7 |
| **P1** | E12-E18 correlation 0 | `analysis.rs:288` | Debug data loader → verify E12/E18 in `satellites` HashMap |
| **P1** | J2 returns zero | `analysis.rs:269` | Debug `ChronoGaugeWitness::solve_j2` and gauge field |
| **P1** | Std dev = 0 | Data loader | Verify `clock_bias_ns` values have variance |
| **P1** | PCA 97% in PC1 | Data loader | Debug `fleet_data.matrix` for value variance |

---

## 6. Verification Checklist

After fixes, verify:

- [ ] Standard deviations are non-zero for clock differentials
- [ ] E12-E18 correlations match reference values (±0.01)
- [ ] J2 inversions produce non-zero values
- [ ] GR coefficients and observed shifts are computed
- [ ] Anomalies show correct satellite IDs (not `E??`)
- [ ] Anomalies show non-zero PC2 percentages
- [ ] 2016 epoch count matches reference (~32256)
- [ ] PCA variance distributions are similar to reference
