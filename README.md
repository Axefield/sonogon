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
│   └── export.ts          # CSV and JSON export utilities
└── validation/
    └── tests.ts           # Validation tests
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

### EM Cavity (`emCavity.ts`)
Enhanced with advanced features:
- **Cross-mode coupling** (mode mixing, energy transfer)
- **Frequency-dependent quality factor** Q(ω)
- **Detailed parametric pumping** (quantum Hamiltonian)
- Mode frequency modulation by bubble radius
- Plasma frequency cutoff effects
- Stored energy pumping and decay

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
- `estimateEmission()` - Total power (blackbody + bremsstrahlung + EM decay)
- `computeSpectralDistribution()` - Wavelength-dependent emission spectrum
- `computeGradients()` - Temporal gradient detection
- `computeSpatialGradients()` - Spatial gradient estimates

### Diagnostics (`diagnostics.ts`)
- `computeEnergyBudget()` - Energy flow tracking (acoustic, thermal, chemical, EM)
- `detectExtremeGradients()` - Non-classical regime detection
- `computePlasmaDiagnostics()` - Plasma frequency, Debye length, mean free path
- `estimateRefractiveIndex()` - Effective refractive index from state
- `computeStabilityMetrics()` - Limit cycle detection, runaway heating checks
- `computeModeSqueezingMetrics()` - Parametric amplification analysis

### Initial States (`initialStates.ts`)
- `createEquilibriumState()` - Generate equilibrium initial state
- `createExpandedState()` - State at maximum bubble expansion
- `createCollapseState()` - State just before collapse
- `createFromPreset()` - State based on preset parameters
- `validateInitialState()` - Consistency checks

### Simulation Runner (`simulation/runner.ts`)
- `runSimulation()` - High-level API with automatic analysis
- `quickSimulation()` - Simplified entry point
- Automatic time series generation
- Built-in analysis options

### Export Utilities (`io/export.ts`)
- `exportToCSV()` - Export time series to CSV
- `exportToJSON()` - Export complete state history
- `exportObservables()` - Export specific observables

### Validation Tests (`validation/`)
- `tests.ts` - Adiabatic compression, Saha equilibrium, energy conservation
- `kinematicTests.ts` - Shape oscillations, translation, coupling tests
- `runAllTests()` - Run all validation tests

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

// Run simulation with adaptive integrator
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
  analysis: {
    computeEmission: true,
    computeSpectrum: true,
    computeGradients: true,
    computeEnergyBudget: true,
  },
});

console.log(`Simulated ${result.timeSeries.length} time steps`);
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
- **`DETAILED_ENHANCEMENTS_COMPLETE.md`**: Detailed physics enhancements (13 enhancements)
- **`KINEMATIC_ENHANCEMENTS_COMPLETE.md`**: Kinematic sciences features (shape, translation)
- **`ADVANCED_FEATURES_COMPLETE.md`**: Van der Waals iterative and enhanced shape oscillations
- **`ADDITIONAL_ENHANCEMENTS_COMPLETE.md`**: Three-body collisions, diffusion, parametric pumping
- **`FINAL_ENHANCEMENTS.md`**: Complete enhancement summary

## Validation

Run validation tests to verify model correctness:

```typescript
import { runAllTests } from './src/sonoluminescence/validation/tests';
import { runAllKinematicTests } from './src/sonoluminescence/validation/kinematicTests';

// Run physics validation tests
const physicsTests = runAllTests();
console.log(`Physics tests: ${physicsTests.passed}/${physicsTests.results.length} passed`);

// Run kinematic validation tests
const kinematicTests = runAllKinematicTests();
console.log(`Kinematic tests: ${kinematicTests.passed}/${kinematicTests.results.length} passed`);
```

## Status

✅ **All features fully integrated and tested**
- 41 state dimensions
- 8 physics modules
- Advanced kinematic sciences
- Detailed physics throughout
- Comprehensive validation suite

## License

MIT

## Contributing

Contributions welcome! Please ensure all code follows TypeScript strict mode and includes appropriate documentation.

