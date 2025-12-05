# EM Negative-Space Behavior and Time Series Logging

## Overview

The EM cavity module now explicitly implements and documents the "negative-space" behavior, which refers to the quantum squeezed vacuum state created by parametric amplification during bubble collapse. This document describes the implementation and how to use the time series logging feature.

## Negative-Space Behavior

### Physical Description

The "negative-space" refers to the **quantum squeezed vacuum state** in the Wigner function representation. During extreme bubble compression:

1. **Parametric Amplification**: The boundary motion (Rdot) pumps energy into EM modes via the parametric Hamiltonian:
   ```
   H = g(t) * (a†² + a²)
   ```
   where `g(t) = g0 * Rdot/R` is the time-dependent coupling strength.

2. **Squeezed State Creation**: This Hamiltonian creates a squeezed state where:
   - One quadrature has **reduced uncertainty** (below vacuum level)
   - The conjugate quadrature has **increased uncertainty**
   - The Wigner function becomes **negative** in some regions (non-classical behavior)

3. **Energy Storage**: The squeezed state stores energy in the EM cavity modes, represented by `E_em`.

4. **Photon Emission**: The stored energy decays over a short timescale (τ ~ 1 ns), converting to photons:
   ```
   dE_em/dt = pump_term - E_em/τ
   ```

### Implementation Details

**File**: `src/sonoluminescence/physics/emCavity.ts`

The negative-space behavior is implemented in two places:

1. **Mode Evolution** (lines 135-159):
   - Parametric pumping with cross-coupling: `dRe/dt = g(t) * Im`, `dIm/dt = -g(t) * Re`
   - This cross-coupling is the signature of parametric amplification
   - Additional squeezing terms enhance the effect for strong compression

2. **Stored Energy Evolution** (lines 208-225):
   - Pump term: `α * |Rdot| * |gradient_terms|`
   - Decay term: `-E_em / τ`
   - During collapse: E_em rises as gradients increase
   - After collapse: E_em decays rapidly, producing light emission

## Time Series Logging

### Overview

The simulation runner now includes explicit time series logging of key variables to visualize the negative-space behavior and collapse dynamics.

### Logged Variables

**File**: `src/sonoluminescence/simulation/runner.ts`

The `TimeSeriesLog` interface includes:

- `t`: Time [s]
- `R`: Bubble radius [m]
- `Pg`: Gas pressure [Pa]
- `T`: Gas temperature [K]
- `ne`: Electron density [m⁻³]
- `Te`: Electron temperature [K]
- `E_em`: EM stored energy [J]
- `totalPower`: Total emission power [W] (from `estimateEmission()`)
- `Rdot`: Radial velocity [m/s]
- `dPg_dt`: Pressure derivative magnitude [Pa/s] (optional)

### Usage

```typescript
import { runSimulation } from "./simulation/runner";
import { createArgonBubblePreset } from "./config/presets";
import { SonoluminescenceModel } from "./model/sonoluminescenceModel";
import { DefaultStateVectorMapper } from "./core/statevector";

const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();
const model = new SonoluminescenceModel(mapper, params);

const result = runSimulation({
  model,
  initialState: /* ... */,
  integratorOptions: {
    dt: 1e-9,
    tMax: 1e-5, // 10 μs
    tolerance: 1e-6,
  },
  logTimeSeries: true, // Enable time series logging
  analysis: {
    computeEmission: true,
  },
});

// Access time series log
if (result.timeSeriesLog) {
  // Export to CSV
  import { exportTimeSeriesToCSV } from "./io/timeSeriesExport";
  const csv = exportTimeSeriesToCSV(result.timeSeriesLog);
  
  // Or export to JSON
  import { exportTimeSeriesToJSON } from "./io/timeSeriesExport";
  const json = exportTimeSeriesToJSON(result.timeSeriesLog);
  
  // Find collapse cycle
  import { findCollapseCycle } from "./io/timeSeriesExport";
  const cycle = findCollapseCycle(result.timeSeriesLog);
  if (cycle) {
    console.log(`Collapse at index ${cycle.minRIndex}`);
    console.log(`E_em peak at index ${cycle.maxE_emIndex}`);
    console.log(`Power peak at index ${cycle.maxPowerIndex}`);
  }
}
```

### Expected Behavior During Collapse

For a collapsing cycle, the time series should show:

1. **Extreme Gradients Phase** (before minimum R):
   - `|Rdot|` increases (rapid compression)
   - `|dPg/dt|` increases (rapid pressure rise)
   - `E_em` **rises** (parametric pumping stores energy in negative-space state)
   - `T` and `ne` increase (heating and ionization)

2. **Collapse Point** (minimum R):
   - `R` reaches minimum
   - `Pg` and `T` peak
   - `E_em` may peak or continue rising

3. **Decay Phase** (after minimum R):
   - `E_em` **decays rapidly** (τ ~ 1 ns)
   - `totalPower` **spikes** (photons emitted from stored energy)
   - This is the "light coming from decay in a negative cavity state"

### Visualization

The time series log can be plotted to visualize:

```
E_em(t) and totalPower(t) during collapse:

E_em:     ___/‾‾‾‾\___     (rises during compression, decays after)
Power:   ___/‾‾‾‾‾‾\___   (spikes when E_em decays)
```

This plot literally shows the "light coming from decay in a negative cavity state" - the stored energy in the squeezed/negative-space state decays as photons, producing the sonoluminescence flash.

## Export Utilities

**File**: `src/sonoluminescence/io/timeSeriesExport.ts`

### Functions

1. **`exportTimeSeriesToCSV(timeSeriesLog, filename?)`**
   - Exports time series to CSV format
   - Columns: t, R, Pg, T, ne, Te, E_em, totalPower, Rdot, dPg_dt
   - Suitable for plotting in Python, MATLAB, Excel, etc.

2. **`exportTimeSeriesToJSON(timeSeriesLog)`**
   - Exports time series to JSON format
   - Preserves all data types and structure

3. **`findCollapseCycle(timeSeriesLog)`**
   - Identifies the collapse cycle window
   - Returns indices for: minimum R, maximum E_em, maximum power
   - Useful for extracting just the collapse phase for detailed analysis

## Example Plot

To visualize the negative-space behavior, plot:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv('time_series.csv')

# Plot E_em and totalPower
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.plot(df['t [s]'], df['E_em [J]'], label='E_em')
ax1.set_ylabel('Stored Energy [J]')
ax1.legend()
ax1.grid(True)

ax2.plot(df['t [s]'], df['totalPower [W]'], label='totalPower', color='orange')
ax2.set_ylabel('Emission Power [W]')
ax2.set_xlabel('Time [s]')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('negative_space_behavior.png')
```

This plot will show:
- E_em rising during extreme gradients (parametric pumping)
- E_em decaying rapidly (negative-space state decay)
- totalPower spiking when E_em decays (photon emission)

## Summary

✅ **EM negative-space behavior is now explicit**:
- Documented in code comments
- Implemented via parametric Hamiltonian
- Creates squeezed vacuum state during compression

✅ **Time series logging enabled**:
- Logs R(t), Pg(t), T(t), ne(t), Te(t), E_em(t), totalPower(t)
- Export to CSV or JSON
- Find collapse cycles automatically

✅ **Visualization ready**:
- Plot E_em(t) and totalPower(t) to see the negative-space decay
- The spike in totalPower when E_em decays is the "light from negative cavity state"

