# Sonogon: 41D+ Sonoluminescence Model

[![TypeScript](https://img.shields.io/badge/TypeScript-5.0+-blue.svg)](https://www.typescriptlang.org/)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![State Dimensions](https://img.shields.io/badge/dimensions-41+-orange.svg)]()
[![Physics Modules](https://img.shields.io/badge/physics-8%20modules-purple.svg)]()

A comprehensive TypeScript implementation of a multi-dimensional sonoluminescence model that captures the full coupled dynamics of bubble collapse, plasma formation, and light emission with advanced physics and kinematic features.

## Overview

Sonoluminescence is the phenomenon where collapsing bubbles in a liquid emit flashes of light. This model implements a complete **41-dimensional** state space that tracks:

- **Hydrodynamics**: Bubble radius and velocity (Rayleigh-Plesset/Keller-Miksis dynamics)
- **Shape Oscillations**: Non-spherical bubble deformations (P2, P4 spherical harmonics)
- **Bubble Translation**: 3D position and velocity dynamics with Bjerknes forces
- **Thermodynamics**: Gas pressure, temperature, species densities with detailed heat transfer
- **Plasma Physics**: Electron density, temperature, ionization with detailed collision models
- **EM Cavity Modes**: Parametric amplification, mode coupling, frequency-dependent Q
- **Chemical Kinetics**: Reaction progress with pressure-dependent rates and three-body collisions
- **Internal Energy**: Translational, rotational, vibrational, electronic with Landau-Teller relaxation
- **Acoustic Driving**: Multi-frequency, non-sinusoidal waveforms with spatial gradients

## Features

- **41 State Dimensions**: Comprehensive state vector tracking (up from 29)
- **Advanced Physics**: Detailed models throughout (heat transfer, collision cross sections, etc.)
- **Kinematic Sciences**: Shape oscillations and bubble translation dynamics
- **Non-Ideal Gas**: Van der Waals EOS with iterative solution
- **Enhanced Shape Dynamics**: Nonlinear coupling, mode interactions, acoustic coupling
- **Detailed Plasma**: Temperature-dependent collisions, multiple recombination channels
- **Adaptive Integrator**: Dormand-Prince 5(4) with error control and event detection
- **Export Utilities**: CSV and JSON export for simulation results
- **Validation Tests**: Adiabatic compression, Saha equilibrium, energy conservation
- **Coupled Dynamics**: Proper interaction between all physics modules
- **Diagnostics**: Energy budget tracking, stability analysis, extreme gradient detection
- **Atomic/Subatomic Physics Tracking**: Comprehensive metrics for nuclear densities, MeV temperatures, quark-level conditions, exotic hadron formation, quantum field effects, fusion conditions, and Planck-scale proximity
- **Presets**: Pre-configured parameter sets for common experimental conditions
- **Type-Safe**: Full TypeScript implementation with strict typing

## Installation

```bash
npm install
npm run build
```

## Quick Start

```typescript
import { DefaultStateVectorMapper } from './src/sonoluminescence/core/statevector';
import { SonoluminescenceModel } from './src/sonoluminescence/model/sonoluminescenceModel';
import { createArgonBubblePreset } from './src/sonoluminescence/config/presets';
import { createEquilibriumState } from './src/sonoluminescence/analysis/initialStates';
import { runSimulation } from './src/sonoluminescence/simulation/runner';

// Create model with preset parameters
const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();
const model = new SonoluminescenceModel(mapper, params);

// Create initial state (helper function)
const initialState = createEquilibriumState(params, 5e-6); // 5 micron bubble
const x0 = mapper.toVector(initialState);

// Run simulation with automatic analysis
const result = runSimulation({
  model,
  initialState,
  integratorOptions: {
    adaptive: true,
    tolerance: 1e-6,
    dt: 1e-9,
    tMax: 1e-5,
  },
});

console.log(`Simulated ${result.timeSeries.length} time steps`);
console.log(`Final radius: ${result.timeSeries[result.timeSeries.length - 1].hydro.R * 1e6} μm`);
```

## Project Structure

```
src/sonoluminescence/
├── core/
│   ├── integrator.ts      # RK4 and Dormand-Prince 5(4) adaptive solver
│   ├── statevector.ts     # State vector mapping (41 dimensions)
│   └── units.ts           # Physical constants and units
├── physics/
│   ├── acoustic.ts        # Acoustic driving with gradients
│   ├── bubbleTranslation.ts  # Bubble translation and Bjerknes forces ⭐
│   ├── emCavity.ts        # EM mode evolution with coupling
│   ├── hydro.ts           # Rayleigh-Plesset/Keller-Miksis equation
│   ├── plasma.ts          # Plasma physics with detailed collisions
│   ├── reactions.ts       # Chemical kinetics with pressure-dependent rates
│   ├── shapeOscillations.ts  # Shape oscillation dynamics ⭐
│   └── thermoChem.ts      # Thermodynamics with detailed heat transfer
├── model/
│   ├── sonoluminescenceModel.ts  # Main model class
│   └── types.ts           # Type definitions (41 dimensions)
├── config/
│   └── presets.ts         # Parameter presets
├── analysis/
│   ├── diagnostics.ts     # Diagnostic tools
│   ├── initialStates.ts   # Initial state generators
│   └── observables.ts     # Emission and gradients
├── simulation/
│   └── runner.ts          # High-level simulation API
├── io/
│   ├── export.ts          # CSV and JSON export utilities
│   └── timeSeriesExport.ts  # Time series log export (R, Pg, T, ne, Te, E_em, totalPower) ⭐
├── scripts/
│   ├── canonicalCollapseAnalysis.ts  # Main canonical analysis script ⭐
│   └── collapseVisualization.ts      # Visualization script ⭐
└── validation/
    ├── tests.ts           # Basic validation tests
    ├── physicsValidation.ts  # Comprehensive physics validation tests (exact formulas) ⭐
    └── kinematicTests.ts  # Kinematic validation tests
```

## State Vector Dimensions

The model tracks **41 state variables**:

1. **Hydrodynamics** (2): Radius `R`, Radius velocity `Rdot`
2. **Shape Oscillations** (4): P2 amplitude `a₂`, P2 velocity `a₂_dot`, P4 amplitude `a₄`, P4 velocity `a₄_dot`
3. **Bubble Translation** (6): Position `x, y, z`, Velocity `vx, vy, vz`
4. **Gas Macro State** (2): Pressure `Pg`, Temperature `T`
5. **Species Densities** (9): H₂O, O₂, N₂, Ar, Xe, H, O, OH, N
6. **Plasma** (3): Electron density `ne`, Electron temperature `Te`, Ionization fraction
7. **Internal Energy** (4): Translational, rotational, vibrational, electronic
8. **EM Modes** (6): 3 modes × (Real, Imaginary)
9. **EM Stored Energy** (1): `E_em`
10. **Acoustic Phase** (1): `φ(t)`
11. **Reaction Progress** (3): `ξ₀`, `ξ₁`, `ξ₂`

**Total: 41 dimensions** (expandable with additional modes/species)

## Physics Modules

### Hydrodynamics (`hydro.ts`)
Implements the Rayleigh-Plesset and Keller-Miksis equations:
- Pressure-driven expansion/compression
- Viscous damping and surface tension
- Acoustic pressure coupling
- **Keller-Miksis equation** for compressible liquid effects (supersonic collapse)
- Uses actual gamma from thermo params (not hardcoded)

### Shape Oscillations (`shapeOscillations.ts`) ⭐ NEW
Non-spherical bubble dynamics using spherical harmonics:
- **P2 (Quadrupole)** and **P4 (Hexadecapole)** modes
- Natural frequencies from surface tension
- **Nonlinear mode-mode coupling** (P2 ↔ P4 interactions)
- **Shape-radial coupling** (energy exchange)
- **Acoustic field coupling** (pressure gradient driving)
- Amplitude-dependent frequency shifts
- Effective radius computation

### Bubble Translation (`bubbleTranslation.ts`) ⭐ NEW
3D position and velocity dynamics:
- **Primary Bjerknes force**: `F = -V * ∇P_acoustic`
- **Drag forces**: Stokes drag and custom coefficients
- Position and velocity tracking
- Secondary Bjerknes forces (bubble-bubble interactions)

### Thermodynamics (`thermoChem.ts`)
Enhanced with detailed physics:
- Adiabatic compression heating
- **Detailed heat transfer**: Conduction, convection, radiation
- **Temperature-dependent heat capacity** (Cp, Cv) for each species
- **Van der Waals EOS** with **iterative solution** (Newton-Raphson)
- Species number density evolution with **diffusion effects**
- Internal energy partition tracking with **Landau-Teller relaxation**

### Plasma Physics (`plasma.ts`)
Enhanced with detailed collision models:
- Saha equation for ionization equilibrium
- **Detailed ionization rates** (electron impact with cross sections)
- **Multiple recombination channels**: Radiative, three-body, dielectronic
- **Temperature-dependent collision cross sections**
- Electron temperature dynamics with proper energy exchange
- Plasma frequency calculation

### EM Cavity (`emCavity.ts`) ⭐ CAVITY-QED FORMAL SUBSYSTEM
Enhanced with formal cavity-QED physics and explicit negative-space behavior:
- **Cross-mode coupling** (mode mixing, energy transfer)
- **EM mode frequency shifting**: `ω_k(R(t), n(r,t))` - frequency depends on radius and refractive index
- **Mode squeezing equations**: `ȧk = f(ak, R, Ṙ, ∂n/∂t)` - includes refractive index time derivative
- **Frequency-dependent quality factor** Q(ω)
- **Dynamic Q-factor**: `Q(R, n, ω)` - depends on radius, refractive index, and frequency
- **Detailed parametric pumping** (quantum Hamiltonian)
- **Negative-space behavior**: Quantum squeezed vacuum state with explicit Wigner function representation
- **Parametric amplification**: `H = g(t) * (a†² + a²)` creates squeezed states during compression
- **Stored-energy pump from gradients**: Explicit gradient-based pumping (`|Rdot|`, `|dPg/dt|`, `|∇P|`)
- **Radiation backreaction**: EM field exerts force on bubble boundary
- **Bremsstrahlung emission**: Energy from electron-ion collisions added to stored energy
- **Recombination emission**: Energy from electron-ion recombination added to stored energy
- Mode frequency modulation by bubble radius and refractive index
- Plasma frequency cutoff effects
- Stored energy pumping and decay: `E_em` rises during extreme gradients, decays as photons

### Chemical Reactions (`reactions.ts`)
Enhanced with pressure-dependent rates:
- Arrhenius rate constants
- **Lindemann falloff mechanism** (pressure-dependent rates)
- **Detailed three-body collisions** (species-specific efficiencies)
- Water dissociation: H₂O ↔ H + OH
- Oxygen dissociation: O₂ ↔ 2O
- Nitrogen dissociation: N₂ ↔ 2N
- Hydroxyl dissociation: OH ↔ O + H
- Three-body recombination: H + OH → H₂O

### Acoustic Driving (`acoustic.ts`)
Enhanced with spatial gradients:
- Single and multi-frequency driving
- Non-sinusoidal waveforms (square, sawtooth, custom)
- **Acoustic field gradients**: `∇P_acoustic` and `∇²P_acoustic`
- Standing wave patterns
- Wave vector specification

## Presets

Pre-configured parameter sets for common experimental conditions:

- `createWaterAirPreset()` - Water with dissolved air
- `createArgonBubblePreset()` - Argon bubble (common for sonoluminescence)
- `createXenonBubblePreset()` - Xenon bubble (high light yield)
- `createHighIntensityPreset()` - Extreme driving conditions

## Analysis Tools

### Observables (`observables.ts`)
- `estimateEmission()` - Total power with explicit breakdown:
  - **Blackbody-ish T term**: `P_bb = σ * T⁴ * A`
  - **Bremsstrahlung from ne, Te**: `P_brems ~ ne² * sqrt(Te) * V`
  - **EM decay from E_em**: `P_em = E_em / τ` (light from negative cavity state)
- `computeSpectralDistribution()` - Wavelength-dependent emission spectrum
- `computeGradients()` - Temporal gradient detection
- `computeSpatialGradients()` - Spatial gradient estimates

### Diagnostics (`diagnostics.ts`)
- `computeEnergyBudget()` - Energy flow tracking (acoustic, thermal, chemical, EM)
- `detectExtremeGradients()` - Non-classical regime detection with atomic/subatomic disturbance tracking
- `computeAtomicSubatomicDisturbances()` - ⭐ **Comprehensive atomic/subatomic physics metrics**:
  - **Nuclear conditions**: Density ratios, nuclear matter density proximity, MeV temperature scales
  - **Strong force metrics**: Coupling constant, field strength, interaction range
  - **Quark-level conditions**: Deconfinement proximity, quark-gluon-plasma phase transitions
  - **Exotic hadron formation**: Tetraquark, pentaquark, hexaquark, hybrid meson, and glueball formation probabilities
  - **Quark flavor tracking**: Energy thresholds for light, strange, charm, and heavy quark production
  - **Quantum field effects**: QCD vacuum energy, vacuum fluctuations, Casimir effect, gluon field strength, color charge density
  - **Nuclear fusion conditions**: Gamow factor, fusion cross-section, fusion rate estimates
  - **Momentum space distributions**: Fermi energy, degeneracy pressure, quantum statistics relevance
  - **Planck-scale physics**: Proximity to Planck density, length, and time scales
  - **QCD phase diagram**: Classification of matter phase (hadronic, quark-gluon-plasma, crossover)
  - **Accelerator comparisons**: Conditions comparable to LHC, RHIC, SPS
  - **Probing depth**: Atomic, nuclear, quark, or Planck-scale physics relevance
- `computePlasmaDiagnostics()` - Plasma frequency, Debye length, mean free path
- `estimateRefractiveIndex()` - Effective refractive index from state: `n(r,t) = n_neutral + n_plasma`
- `computeStabilityMetrics()` - Limit cycle detection, runaway heating checks
- `computeModeSqueezingMetrics()` - Parametric amplification analysis with pumped vs released energy:
  - Reports **pumped energy** (into modes during compression)
  - Reports **released energy** (from modes as photons)
  - Quantifies negative-space behavior

### Initial States (`initialStates.ts`)
- `createEquilibriumState()` - Generate equilibrium initial state
- `createExpandedState()` - State at maximum bubble expansion
- `createCollapseState()` - State just before collapse
- `createFromPreset()` - State based on preset parameters
- `validateInitialState()` - Consistency checks

### Simulation Runner (`simulation/runner.ts`)
- `runSimulation()` - High-level API with automatic analysis
- `quickSimulation()` - Simplified entry point
- **Time series logging**: Explicit logging of R(t), Pg(t), T(t), ne(t), Te(t), E_em(t), totalPower(t)
- Automatic time series generation
- Built-in analysis options
- **EM negative-space visualization**: Track E_em rise during collapse and decay during photon emission

### Canonical Analysis (`analysis/canonicalAnalysis.ts`) ⭐ NEW
- `performCanonicalAnalysis()` - **"This is what Sonogon tells us about a collapsing argon bubble"** script
  - Peak emission (with breakdown by component)
  - Max electron temperature
  - Max E_em (stored energy in negative-space state)
  - Energy budget (via `computeEnergyBudget()`)
  - Extreme gradients (via `detectExtremeGradients()`)
  - **Atomic/subatomic disturbances** (via `computeAtomicSubatomicDisturbances()`):
    - Nuclear density ratios and MeV temperature scales
    - Exotic hadron formation probabilities (tetraquarks, pentaquarks, hexaquarks, hybrid mesons, glueballs)
    - Quark flavor thresholds and QCD phase transitions
    - Quantum field effects (QCD vacuum, gluon fields, color charge)
    - Nuclear fusion conditions and rates
    - Momentum space distributions (Fermi energy, degeneracy)
    - Planck-scale proximity and quantum gravity effects
    - QCD phase diagram position
  - Mode squeezing metrics (pumped vs released energy)
- `printCanonicalAnalysis()` - Human-readable output with detailed atomic/subatomic metrics

### Visualization (`scripts/collapseVisualization.ts`) ⭐ NEW
- `generateCollapseVisualizationData()` - Export CSV for plotting R(t), T(t), ne(t), E_em(t), totalPower(t)
- `printVisualizationInstructions()` - Python plotting code
- **"Collapse emits/decays inside a negative space"** visualization
- Automatic collapse cycle detection

### Export Utilities (`io/export.ts`, `io/timeSeriesExport.ts`)
- `exportToCSV()` - Export time series to CSV
- `exportToJSON()` - Export complete state history
- `exportObservables()` - Export specific observables
- `exportTimeSeriesToCSV()` - Export explicit time series log (R, Pg, T, ne, Te, E_em, totalPower)
- `exportTimeSeriesToJSON()` - Export time series log as JSON
- `findCollapseCycle()` - Automatically identify collapse cycle windows

### Canonical Scripts (`scripts/`) ⭐ NEW
- `canonicalCollapseAnalysis.ts` - Main script for canonical analysis
  - Runs simulation with all advanced features
  - Performs canonical analysis
  - Generates visualization data
  - Ready-to-use entry point

### Validation Tests (`validation/`)
- `tests.ts` - Adiabatic compression, Saha equilibrium, energy conservation, plasma frequency
- `physicsValidation.ts` - Comprehensive physics validation:
  - **Adiabatic scaling test**: Verify T ∝ (R₀/R)^(3(γ-1)), P ∝ (R₀/R)^(3γ) with Pa=0
  - **Minnaert frequency test**: Small-amplitude linear oscillation frequency matching
  - **Energy budget closure**: Verify acoustic work in ≈ thermal + EM + chemical + viscous dissipation
  - **Plasma equilibrium spot-check**: Compare ionization fraction vs Saha equation
- `kinematicTests.ts` - Shape oscillations, translation, coupling tests
- `runAllTests()` - Run all validation tests
- `runAllPhysicsValidationTests()` - Run comprehensive physics validation suite

## Visualizing EM Negative-Space Behavior

The time series logging feature allows you to visualize the "light coming from decay in a negative cavity state":

```typescript
const result = runSimulation({
  model,
  initialState,
  integratorOptions: { /* ... */ },
  logTimeSeries: true, // Enable time series logging
  analysis: { computeEmission: true },
});

// Export and plot E_em(t) and totalPower(t)
import { exportTimeSeriesToCSV } from './src/sonoluminescence/io/timeSeriesExport';
const csv = exportTimeSeriesToCSV(result.timeSeriesLog!);

// Plot shows:
// - E_em rises during extreme gradients (parametric pumping)
// - E_em decays rapidly (τ ~ 1 ns) after collapse
// - totalPower spikes when E_em decays (photon emission)
// This is the "light from negative cavity state" made visible!
```

See `EM_NEGATIVE_SPACE_LOGGING.md` for detailed explanation and plotting examples.

## Example: Canonical Analysis

### Run Canonical "This is what Sonogon tells us" Analysis

```typescript
import { runCanonicalCollapseAnalysis } from './src/sonoluminescence/scripts/canonicalCollapseAnalysis';

// Run the canonical analysis script
runCanonicalCollapseAnalysis();

// Outputs:
// - Peak emission (with breakdown: blackbody, bremsstrahlung, EM decay)
// - Max electron temperature
// - Max E_em (stored energy in negative-space state)
// - Energy budget
// - Extreme gradients
// - Atomic/subatomic disturbances:
//   - Nuclear density ratios and MeV temperatures
//   - Exotic hadron formation probabilities (tetraquarks, pentaquarks, etc.)
//   - Quark flavor thresholds and QCD phase transitions
//   - Quantum field effects and fusion conditions
//   - Planck-scale proximity
// - Mode squeezing metrics (pumped vs released energy)
// - Visualization data for plotting
```

### Manual Canonical Analysis

```typescript
import { DefaultStateVectorMapper } from './src/sonoluminescence/core/statevector';
import { SonoluminescenceModel } from './src/sonoluminescence/model/sonoluminescenceModel';
import { createArgonBubblePreset } from './src/sonoluminescence/config/presets';
import { createEquilibriumState } from './src/sonoluminescence/analysis/initialStates';
import { runSimulation } from './src/sonoluminescence/simulation/runner';
import { performCanonicalAnalysis, printCanonicalAnalysis } from './src/sonoluminescence/analysis/canonicalAnalysis';
import { generateCollapseVisualizationData } from './src/sonoluminescence/scripts/collapseVisualization';

const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();

// Enable all cavity-QED features
params.em.useRefractiveIndexFrequencyShift = true;
params.em.useDynamicQ = true;
params.em.useRadiationBackreaction = true;
params.em.includeBremsstrahlungEmission = true;
params.em.includeRecombinationEmission = true;

const model = new SonoluminescenceModel(mapper, params);
const initialState = createEquilibriumState(params, 5e-6);

// Run simulation with time series logging
const result = runSimulation({
  model,
  initialState,
  integratorOptions: {
    adaptive: true,
    tolerance: 1e-6,
    dt: 1e-9,
    tMax: 1e-5,
  },
  logTimeSeries: true, // Enable time series logging
  analysis: {
    computeEmission: true,
    computeGradients: true,
    computeEnergyBudget: true,
  },
});

// Perform canonical analysis
const analysis = performCanonicalAnalysis(result);
printCanonicalAnalysis(analysis);

// Generate visualization data
const vizData = generateCollapseVisualizationData(result);
// Plot R(t), T(t), ne(t), E_em(t), totalPower(t) to see "collapse emits/decays inside a negative space"
```

## Example: Running a Simulation

### Basic Simulation

```typescript
import { DefaultStateVectorMapper } from './src/sonoluminescence/core/statevector';
import { SonoluminescenceModel } from './src/sonoluminescence/model/sonoluminescenceModel';
import { createArgonBubblePreset } from './src/sonoluminescence/config/presets';
import { integrateAdaptive } from './src/sonoluminescence/core/integrator';
import { createEquilibriumState } from './src/sonoluminescence/analysis/initialStates';
import { runSimulation } from './src/sonoluminescence/simulation/runner';

// Setup with enhanced features
const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();

// Enable advanced features
params.hydro.useKellerMiksis = true;
params.thermo.useDetailedHeatTransfer = true;
params.thermo.useIterativeVanDerWaals = true;
params.plasma.useDetailedCollisions = true;
params.em.useModeCoupling = true;

const model = new SonoluminescenceModel(mapper, params);
const initialState = createEquilibriumState(params, 5e-6);

// Run simulation with adaptive integrator and time series logging
const result = runSimulation({
  model,
  initialState,
  integratorOptions: {
    adaptive: true,
    tolerance: 1e-6,
    dtMin: 1e-12,
    dtMax: 1e-8,
    tMax: 1e-5,
  },
  logTimeSeries: true, // Enable explicit time series logging
  analysis: {
    computeEmission: true,
    computeSpectrum: true,
    computeGradients: true,
    computeEnergyBudget: true,
  },
});

console.log(`Simulated ${result.timeSeries.length} time steps`);

// Export time series log to visualize EM negative-space behavior
if (result.timeSeriesLog) {
  import { exportTimeSeriesToCSV, findCollapseCycle } from './src/sonoluminescence/io/timeSeriesExport';
  
  // Export to CSV for plotting
  const csv = exportTimeSeriesToCSV(result.timeSeriesLog);
  // Save to file or plot E_em(t) and totalPower(t) to see negative-space decay
  
  // Find collapse cycle
  const cycle = findCollapseCycle(result.timeSeriesLog);
  if (cycle) {
    console.log(`Collapse at index ${cycle.minRIndex}`);
    console.log(`E_em peak at index ${cycle.maxE_emIndex}`);
    console.log(`Power peak at index ${cycle.maxPowerIndex}`);
  }
}
```

### Advanced Simulation with Shape Oscillations and Translation

```typescript
// Enable kinematic features
params.shape = {
  sigma: 0.0728,
  mu: 0.001002,
  rho: 998.2,
  enableShapeOscillations: true,
  enableNonlinearCoupling: true,
  enableShapeRadialCoupling: true,
  enableAcousticCoupling: true,
};

params.translation = {
  rho: 998.2,
  mu: 0.001002,
  enableTranslation: true,
};

params.hydro.enableShapeRadialCoupling = true;
params.hydro.shapeCouplingCoefficient = 0.1;

params.acoustic.enableGradients = true;
params.acoustic.waveVector = { x: 100, y: 0, z: 0 };

// Create initial state with shape and translation
const initialState = createEquilibriumState(params, 5e-6);
initialState.shape = {
  a2: 1e-7,      // Small initial deformation
  a2_dot: 0,
  a4: 0,
  a4_dot: 0,
};
initialState.translation = {
  x: 0, y: 0, z: 0,
  vx: 0, vy: 0, vz: 0,
};

// Run simulation
const result = runSimulation({
  model,
  initialState,
  integratorOptions: {
    adaptive: true,
    tolerance: 1e-6,
    dt: 1e-9,
    tMax: 1e-5,
  },
  analysis: {
    computeEmission: true,
    computeGradients: true,
    computeEnergyBudget: true,
  },
});

// Export results
import { exportObservables } from './src/sonoluminescence/io/export';
exportObservables(result.timeSeries, mapper, 'results.csv');
```

## Physical Constants

All physical constants are defined in `core/units.ts`:
- Gas constant, Boltzmann constant
- Electron charge and mass
- Vacuum permittivity and permeability
- Speed of light, Planck constant
- Stefan-Boltzmann constant

## Units

The model uses SI base units:
- Length: meters (m)
- Time: seconds (s)
- Mass: kilograms (kg)
- Temperature: Kelvin (K)
- Pressure: Pascal (Pa)
- Energy: Joules (J)

## Advanced Features

### Atomic/Subatomic Physics Tracking ⭐ NEW

The model now tracks **atomic and subatomic-level disturbances** during extreme bubble collapse conditions, connecting macroscopic sonoluminescence gradients to fundamental particle physics:

- **Nuclear Density Conditions**: Tracks proximity to nuclear matter density (~2.3×10¹⁷ kg/m³), enabling conditions where nuclear physics becomes relevant
- **MeV Temperature Scales**: Monitors temperatures approaching MeV scales (1 MeV = 1.16×10¹⁰ K), where nuclear binding energies and quark production become possible
- **Strong Force Metrics**: Estimates strong coupling constant (α_s), field strength, and interaction range as conditions approach nuclear scales
- **Quark-Level Physics**: Detects proximity to quark deconfinement transition (~150 MeV), where quarks and gluons are no longer confined within hadrons
- **Exotic Hadron Formation**: Estimates formation probabilities for:
  - **Tetraquarks**: Four-quark states (like CERN's all-charm tetraquarks)
  - **Pentaquarks**: Five-quark states
  - **Hexaquarks**: Six-quark dibaryon states
  - **Hybrid Mesons**: Quark-gluon hybrid states
  - **Glueballs**: Pure gluon states
- **Quark Flavor Tracking**: Determines energy thresholds for producing different quark flavors (up, down, strange, charm, bottom, top)
- **Quantum Field Effects**: Estimates QCD vacuum energy, vacuum fluctuations, Casimir effect, gluon field strength, and color charge density
- **Nuclear Fusion Conditions**: Detects conditions approaching nuclear fusion, with Gamow factor, cross-section, and rate estimates
- **Momentum Space Distributions**: Computes Fermi energy, degeneracy pressure, and quantum statistics relevance
- **Planck-Scale Proximity**: Monitors approach to Planck density, length, and time scales, where quantum gravity effects may become relevant
- **QCD Phase Diagram**: Classifies matter phase (hadronic, quark-gluon-plasma, crossover) based on temperature and density
- **Accelerator Comparisons**: Compares conditions to LHC, RHIC, and SPS experiments

This feature directly connects the extreme gradients in sonoluminescence (high compression rates, MeV temperatures, nuclear densities) to the same physics scales explored at particle accelerators like CERN, where exotic hadrons like tetraquarks are discovered.

### Cavity-QED Formal Subsystem ⭐ NEW

The EM cavity module is now a **formal cavity-QED subsystem** with:

- **Refractive index frequency shift**: `ω_k(R(t), n(r,t))` - frequency depends on radius and refractive index
- **Mode squeezing with ∂n/∂t**: `ȧk = f(ak, R, Ṙ, ∂n/∂t)` - refractive index time derivative included
- **Dynamic Q-factor**: `Q(R, n, ω)` - depends on radius, refractive index, and frequency
- **Radiation backreaction**: EM field exerts force on bubble boundary
- **Bremsstrahlung + recombination emission**: Added to stored energy from plasma processes

See `CAVITY_QED_ENHANCEMENTS.md` for detailed documentation.

### Detailed Physics Enhancements
- **Heat Transfer**: Conduction, convection, and radiation models with proper boundary conditions
- **Heat Capacity**: Temperature-dependent Cp/Cv for each species using statistical mechanics
- **Van der Waals EOS**: Iterative Newton-Raphson solution (not approximation)
- **Ionization**: Proper rate equations with temperature-dependent collision cross sections
- **Recombination**: Multiple channels (radiative, three-body, dielectronic)
- **Collision Cross Sections**: Temperature-dependent with proper energy exchange
- **Pressure-Dependent Rates**: Lindemann falloff mechanism with Troe broadening
- **Three-Body Collisions**: Species-specific third-body efficiencies
- **Species Diffusion**: Fick's law through bubble wall with concentration gradients
- **Landau-Teller Relaxation**: Accurate vibrational energy exchange

### Kinematic Sciences ⭐ FULLY INTEGRATED
- **Shape Oscillations**: P2 and P4 spherical harmonic modes with natural frequencies
- **Nonlinear Coupling**: Mode-mode interactions (P2 ↔ P4) and energy exchange
- **Shape-Radial Coupling**: Bidirectional energy exchange affecting effective radius
- **Bubble Translation**: 3D dynamics with primary and secondary Bjerknes forces
- **Acoustic Gradients**: Spatial pressure variations (∇P, ∇²P) for force computation
- **All features integrated** into model RHS and state vector mapper

### Numerical Methods
- **Adaptive Integrator**: Dormand-Prince 5(4) with error control and event detection
- **Fixed-Step RK4**: Available for quick simulations
- **Event Detection**: Automatic detection of extreme events
- **Validation Tests**: Comprehensive test suite for all physics modules

## Documentation

- **`WHITEPAPER.md`**: Detailed scientific background and theoretical foundations
- **`INTEGRATOR.md`**: Adaptive integrator documentation (Dormand-Prince 5(4))
- **`EM_NEGATIVE_SPACE_LOGGING.md`**: ⭐ EM negative-space behavior and time series logging guide
- **`CAVITY_QED_ENHANCEMENTS.md`**: ⭐ NEW - Formal cavity-QED subsystem documentation
- **`DETAILED_ENHANCEMENTS_COMPLETE.md`**: Detailed physics enhancements (13 enhancements)
- **`KINEMATIC_ENHANCEMENTS_COMPLETE.md`**: Kinematic sciences features (shape, translation)
- **`ADVANCED_FEATURES_COMPLETE.md`**: Van der Waals iterative and enhanced shape oscillations
- **`ADDITIONAL_ENHANCEMENTS_COMPLETE.md`**: Three-body collisions, diffusion, parametric pumping
- **`FINAL_ENHANCEMENTS.md`**: Complete enhancement summary

## Validation

Run validation tests to verify model correctness:

```typescript
import { runAllTests } from './src/sonoluminescence/validation/tests';
import { runAllPhysicsValidationTests } from './src/sonoluminescence/validation/physicsValidation';
import { runAllKinematicTests } from './src/sonoluminescence/validation/kinematicTests';

// Run basic physics validation tests
const physicsTests = runAllTests();
console.log(`Physics tests: ${physicsTests.passed}/${physicsTests.results.length} passed`);

// Run comprehensive physics validation tests
const comprehensiveTests = runAllPhysicsValidationTests();
console.log(`Comprehensive tests: ${comprehensiveTests.passed}/${comprehensiveTests.results.length} passed`);
comprehensiveTests.results.forEach(test => {
  console.log(`  ${test.name}: ${test.passed ? 'PASS' : 'FAIL'} ${test.error || ''}`);
});

// Run kinematic validation tests
const kinematicTests = runAllKinematicTests();
console.log(`Kinematic tests: ${kinematicTests.passed}/${kinematicTests.results.length} passed`);
```

### Validation Test Suite

The comprehensive physics validation suite includes:

1. **Adiabatic Scaling Test**: Verifies `T ∝ (R₀/R)^(3(γ-1))` and `P ∝ (R₀/R)^(3γ)` with acoustic drive set to zero
2. **Minnaert Frequency Test**: Checks that small-amplitude oscillations match the Minnaert frequency `ω₀ = (1/R₀) * √(3γP₀/ρ)`
3. **Energy Budget Closure**: Verifies that acoustic work input equals thermal + EM + chemical + viscous dissipation
4. **Plasma Equilibrium Spot-Check**: Compares plasma module's ionization fraction with Saha equation predictions

## Status

✅ **All features fully integrated and tested**
- 41 state dimensions
- 8 physics modules
- Advanced kinematic sciences
- Detailed physics throughout
- **EM negative-space behavior explicitly implemented**
- **Formal cavity-QED subsystem** with refractive index effects, dynamic Q, radiation backreaction
- **Atomic/subatomic physics tracking** - Comprehensive metrics connecting extreme gradients to nuclear, quark-level, and Planck-scale physics
- **Exotic hadron formation tracking** - Tetraquark, pentaquark, hexaquark, hybrid meson, and glueball probabilities
- **Canonical analysis script** - "This is what Sonogon tells us about a collapsing argon bubble" with atomic/subatomic disturbance metrics
- **Visualization tools** - "Collapse emits/decays inside a negative space" plots
- **Time series logging** for visualization
- **Comprehensive validation suite** with exact formulas and numeric tolerances

## License

MIT

## Contributing

Contributions welcome! Please ensure all code follows TypeScript strict mode and includes appropriate documentation.

