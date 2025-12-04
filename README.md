# Sonogon: 36D+ Sonoluminescence Model

A comprehensive TypeScript implementation of a multi-dimensional sonoluminescence model that captures the full coupled dynamics of bubble collapse, plasma formation, and light emission.

## Overview

Sonoluminescence is the phenomenon where collapsing bubbles in a liquid emit flashes of light. This model implements a complete 36+ dimensional state space that tracks:

- **Hydrodynamics**: Bubble radius and velocity (Rayleigh-Plesset dynamics)
- **Thermodynamics**: Gas pressure, temperature, and species number densities
- **Plasma Physics**: Electron density, temperature, and ionization via Saha equation
- **EM Cavity Modes**: Parametric amplification and mode squeezing
- **Chemical Kinetics**: Reaction progress variables for dissociation/recombination
- **Internal Energy**: Translational, rotational, vibrational, and electronic partitions
- **Acoustic Driving**: Phase tracking for acoustic pressure modulation

## Features

- **Complete Physics**: All major physical processes are modeled
- **29 State Dimensions**: Comprehensive state vector tracking
- **Coupled Dynamics**: Proper interaction between all physics modules
- **Diagnostics**: Energy budget tracking, stability analysis, and extreme gradient detection
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
import { integrateRK4 } from './src/sonoluminescence/core/integrator';

// Create model with preset parameters
const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();
const model = new SonoluminescenceModel(mapper, params);

// Set initial state
const initialState = createInitialState(); // Your initial bubble state
const x0 = mapper.toVector(initialState);

// Integrate
const result = integrateRK4(
  (t, x) => model.rhs(t, x),
  x0,
  { dt: 1e-9, tMax: 1e-5 }
);

// Analyze results
console.log(`Simulated ${result.t.length} time steps`);
```

## Project Structure

```
src/sonoluminescence/
├── core/
│   ├── integrator.ts      # RK4 ODE solver
│   ├── statevector.ts     # State vector mapping
│   └── units.ts           # Physical constants and units
├── physics/
│   ├── acoustic.ts        # Acoustic driving
│   ├── emCavity.ts        # EM mode evolution
│   ├── hydro.ts           # Rayleigh-Plesset equation
│   ├── plasma.ts          # Plasma physics (Saha equation)
│   ├── reactions.ts       # Chemical kinetics
│   └── thermoChem.ts      # Thermodynamics and species
├── model/
│   ├── sonoluminescenceModel.ts  # Main model class
│   └── types.ts           # Type definitions
├── config/
│   └── presets.ts         # Parameter presets
└── analysis/
    ├── diagnostics.ts     # Diagnostic tools
    └── observables.ts      # Emission and gradients
```

## State Vector Dimensions

The model tracks 29 state variables:

1. **Hydrodynamics** (2): Radius `R`, Radius velocity `Rdot`
2. **Gas Macro State** (2): Pressure `Pg`, Temperature `T`
3. **Species Densities** (7): H₂O, O₂, N₂, Ar, Xe, H, O
4. **Plasma** (3): Electron density `ne`, Electron temperature `Te`, Ionization fraction
5. **Internal Energy** (4): Translational, rotational, vibrational, electronic
6. **EM Modes** (6): 3 modes × (Real, Imaginary)
7. **EM Stored Energy** (1): `E_em`
8. **Acoustic Phase** (1): `φ(t)`
9. **Reaction Progress** (3): `ξ₀`, `ξ₁`, `ξ₂`

## Physics Modules

### Hydrodynamics (`hydro.ts`)
Implements the Rayleigh-Plesset equation for bubble dynamics:
- Pressure-driven expansion/compression
- Viscous damping
- Surface tension effects
- Acoustic pressure coupling

### Thermodynamics (`thermoChem.ts`)
- Adiabatic compression heating
- Species number density evolution
- Internal energy partition tracking
- Heat loss to surroundings

### Plasma Physics (`plasma.ts`)
- Saha equation for ionization equilibrium
- Electron density evolution
- Electron temperature dynamics
- Plasma frequency calculation

### EM Cavity (`emCavity.ts`)
- Parametric amplification from boundary motion
- Mode frequency modulation by bubble radius
- Plasma frequency cutoff effects
- Stored energy pumping and decay

### Chemical Reactions (`reactions.ts`)
- Arrhenius rate constants
- Water dissociation: H₂O ↔ H + OH
- Oxygen dissociation: O₂ ↔ 2O
- Recombination reactions

## Presets

Pre-configured parameter sets for common experimental conditions:

- `createWaterAirPreset()` - Water with dissolved air
- `createArgonBubblePreset()` - Argon bubble (common for sonoluminescence)
- `createXenonBubblePreset()` - Xenon bubble (high light yield)
- `createHighIntensityPreset()` - Extreme driving conditions

## Analysis Tools

### Observables (`observables.ts`)
- `estimateEmission()` - Total power (blackbody + bremsstrahlung + EM decay)
- `computeGradients()` - Temporal gradient detection
- `computeSpatialGradients()` - Spatial gradient estimates

### Diagnostics (`diagnostics.ts`)
- `computeEnergyBudget()` - Energy flow tracking
- `detectExtremeGradients()` - Non-classical regime detection
- `computePlasmaDiagnostics()` - Plasma frequency, Debye length, mean free path
- `estimateRefractiveIndex()` - Effective refractive index from state
- `computeStabilityMetrics()` - Limit cycle detection, runaway heating checks
- `computeModeSqueezingMetrics()` - Parametric amplification analysis

## Example: Running a Simulation

```typescript
import { DefaultStateVectorMapper } from './src/sonoluminescence/core/statevector';
import { SonoluminescenceModel } from './src/sonoluminescence/model/sonoluminescenceModel';
import { createArgonBubblePreset } from './src/sonoluminescence/config/presets';
import { integrateRK4 } from './src/sonoluminescence/core/integrator';
import { estimateEmission, computeGradients } from './src/sonoluminescence/analysis/observables';
import { computePlasmaDiagnostics } from './src/sonoluminescence/analysis/diagnostics';

// Setup
const mapper = new DefaultStateVectorMapper();
const params = createArgonBubblePreset();
const model = new SonoluminescenceModel(mapper, params);

// Initial state: 5 micron bubble at equilibrium
const initialState = {
  t: 0,
  hydro: { R: 5e-6, Rdot: 0 },
  gas: { Pg: 101325, T: 293.15 },
  species: {
    numberDensity: {
      Ar: 2.5e25, // Argon dominant
      H2O: 1e24,
      O2: 1e23,
      N2: 1e23,
      Xe: 0,
      H: 0,
      O: 0,
    },
  },
  plasma: { ne: 1e20, Te: 10000, ionizationFraction: 0.01 },
  internalEnergy: {
    translational: 1e-18,
    rotational: 1e-19,
    vibrational: 1e-20,
    electronic: 1e-21,
  },
  em: {
    modes: [
      { omega: 2 * Math.PI * 3e14, re: 0, im: 0 },
      { omega: 2 * Math.PI * 4e14, re: 0, im: 0 },
      { omega: 2 * Math.PI * 5e14, re: 0, im: 0 },
    ],
    storedEnergy: 1e-20,
  },
  acoustic: { phase: 0 },
  reactions: { xi: [0, 0, 0] },
};

const x0 = mapper.toVector(initialState);

// Integrate over one acoustic period
const period = 2 * Math.PI / params.acoustic.omega;
const result = integrateRK4(
  (t, x) => model.rhs(t, x),
  x0,
  { dt: 1e-9, tMax: period }
);

// Analyze results
const states = result.x.map((vec, i) => mapper.fromVector(vec, result.t[i]));

// Compute emission at each time step
const emissions = states.map(state => estimateEmission(state));

// Check for extreme gradients
for (let i = 1; i < states.length; i++) {
  const gradients = computeGradients(states[i-1], states[i], result.t[i] - result.t[i-1]);
  if (gradients.extremeGradient) {
    console.log(`Extreme gradient detected at t=${result.t[i]}s`);
  }
}

// Plasma diagnostics at collapse
const minRIndex = states.findIndex(s => s.hydro.R === Math.min(...states.map(s => s.hydro.R)));
const plasmaDiag = computePlasmaDiagnostics(states[minRIndex]);
console.log(`Plasma frequency: ${plasmaDiag.plasmaFrequency / (2 * Math.PI)} Hz`);
console.log(`Debye length: ${plasmaDiag.debyeLength} m`);
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

## References

See `WHITEPAPER.md` for detailed scientific background and theoretical foundations.

## License

MIT

## Contributing

Contributions welcome! Please ensure all code follows TypeScript strict mode and includes appropriate documentation.

