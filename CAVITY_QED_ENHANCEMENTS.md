# Cavity-QED Formal Subsystem Enhancements

## Overview

The EM cavity module has been expanded into a formal cavity-QED style subsystem with comprehensive physics including refractive index effects, dynamic Q-factors, radiation backreaction, and emission additions.

## New Features

### 1. EM Mode Frequency Shifting: ω_k(R(t), n(r,t))

**Implementation**: `computeModeFrequency()` in `emCavity.ts`

The mode frequency now explicitly depends on both bubble radius and refractive index:

```
ω_k(t) = ω_k0 * (R0 / R(t)) * (1 / n_effective)
```

Where `n_effective = n_neutral + n_plasma` combines:
- **Neutral gas refractive index**: `n_neutral ≈ 1 + α * ρ` (density-dependent)
- **Plasma refractive index**: `n²_plasma = 1 - (ωp² / ω²)` (frequency-dependent)

**Usage**:
```typescript
params.em.useRefractiveIndexFrequencyShift = true;
```

### 2. Mode Squeezing Equations: ȧk = f(ak, R, Ṙ, ∂n/∂t)

**Implementation**: Enhanced parametric pumping in `computeEmDerivatives()`

The mode evolution now includes refractive index time derivative:

```
ȧk = -i*ωk*ak + g(t)*ak + (∂ωk/∂n)*(∂n/∂t)*ak + coupling - damping*ak
```

Where:
- `g(t) = g0 * Rdot/R` is the parametric coupling
- `∂ωk/∂n = -ωk/n²` is the frequency shift from refractive index
- `∂n/∂t` is computed from state changes

**Usage**:
```typescript
params.em.useRefractiveIndexTimeDerivative = true;
// Note: Requires statePrev and dt to be passed to computeEmDerivatives()
```

### 3. Stored-Energy Pump from Gradients

**Implementation**: Enhanced pump term calculation

The stored energy pump explicitly uses gradients:

```
dE_em/dt = α * (|Rdot| + |dPg/dt| + |∇P_acoustic|) - E_em/τ
```

Where:
- `|Rdot|`: Radial velocity gradient
- `|dPg/dt|`: Pressure time derivative
- `|∇P_acoustic|`: Acoustic pressure gradient (if available)

This is already implemented and documented in the code.

### 4. Dynamic Q-Factor & Cavity Losses

**Implementation**: `useDynamicQ` parameter

The quality factor Q now depends on:
- **Frequency**: `Q(ω) = Q0 * (ω0/ω)^α`
- **Radius**: `Q(R) ~ R` (larger bubble = higher Q)
- **Refractive index**: `Q(n) ~ n` (higher n = better confinement)

```
Q(R, n, ω) = Q0 * (ω0/ω)^α * (R/R_ref) * (n/n_ref)
```

**Usage**:
```typescript
params.em.useDynamicQ = true;
params.em.Q0 = 1000; // Base Q
```

### 5. Radiation Backreaction Term

**Implementation**: `useRadiationBackreaction` parameter

The EM field exerts a force back on the bubble boundary:

```
F_radiation = -E_em / R * backreaction_coeff
```

This represents radiation pressure from the stored EM energy.

**Usage**:
```typescript
params.em.useRadiationBackreaction = true;
params.em.radiationBackreactionCoeff = 1e-15; // [N·s/J]
```

**Note**: The backreaction force is returned in `EmDerivativesWithState.radiationBackreaction` and can be applied to bubble translation dynamics.

### 6. Bremsstrahlung + Recombination Emission Additions

**Implementation**: `includeBremsstrahlungEmission` and `includeRecombinationEmission` parameters

These add energy to the stored EM energy from plasma processes:

**Bremsstrahlung**:
```
P_brems = ne² * sqrt(Te) * V * prefactor
dE_em/dt += P_brems * dt
```

**Recombination**:
```
P_recomb = ne * n_ions * α_recomb * E_ionization * V
dE_em/dt += P_recomb * dt
```

**Usage**:
```typescript
params.em.includeBremsstrahlungEmission = true;
params.em.includeRecombinationEmission = true;
```

## Complete Stored Energy Equation

The full stored energy evolution is now:

```
dE_em/dt = pump_term - decay_term + bremsstrahlung + recombination
```

Where:
- **pump_term**: `α * |Rdot| * |gradient_terms|` (parametric amplification)
- **decay_term**: `-E_em / τ` (photon emission)
- **bremsstrahlung**: Energy from electron-ion collisions
- **recombination**: Energy from electron-ion recombination

## API Changes

### Enhanced Function Signature

```typescript
export function computeEmDerivatives(
  state: BubbleFullState,
  params: EmCavityParams & { thermoGamma?: number },
  statePrev?: BubbleFullState, // Optional: for computing ∂n/∂t
  dt?: number // Optional: for computing time derivatives
): EmDerivatives | EmDerivativesWithState
```

### New Return Type

```typescript
export interface EmDerivativesWithState {
  dModes: { dRe: number; dIm: number }[];
  dStoredEnergyDt: number;
  radiationBackreaction?: { fx: number; fy: number; fz: number };
  dynamicQ?: number[]; // Q-factor for each mode
}
```

## Usage Example

```typescript
const params = createArgonBubblePreset();

// Enable all cavity-QED features
params.em.useRefractiveIndexFrequencyShift = true;
params.em.useRefractiveIndexTimeDerivative = true;
params.em.useDynamicQ = true;
params.em.useRadiationBackreaction = true;
params.em.includeBremsstrahlungEmission = true;
params.em.includeRecombinationEmission = true;

// For full features, need state history
let statePrev: BubbleFullState | undefined;
for (let i = 1; i < states.length; i++) {
  const dt = times[i] - times[i - 1];
  statePrev = states[i - 1];
  
  const emDeriv = computeEmDerivatives(
    states[i],
    params.em,
    statePrev, // For ∂n/∂t
    dt // For time derivatives
  );
  
  // Access radiation backreaction if enabled
  if ('radiationBackreaction' in emDeriv && emDeriv.radiationBackreaction) {
    // Apply force to bubble translation
    const F = emDeriv.radiationBackreaction;
    // ... apply to bubble dynamics
  }
  
  // Access dynamic Q if enabled
  if ('dynamicQ' in emDeriv && emDeriv.dynamicQ) {
    const Q_values = emDeriv.dynamicQ;
    // ... use for analysis
  }
}
```

## Physical Significance

### Refractive Index Frequency Shift

The refractive index changes during collapse:
- **Compression**: Density increases → n increases → ω decreases
- **Ionization**: Plasma forms → n_plasma changes → ω shifts
- **Temperature**: Affects polarizability → n changes

This creates a **chirped frequency** during collapse, important for parametric amplification.

### Mode Squeezing with ∂n/∂t

The refractive index time derivative contributes to mode evolution:
- **Fast compression**: Large ∂n/∂t → strong frequency modulation
- **Parametric resonance**: When ∂n/∂t matches mode frequency, strong amplification
- **Negative-space state**: Enhanced squeezing from combined Rdot and ∂n/∂t effects

### Dynamic Q-Factor

The quality factor changes during collapse:
- **Small R**: Lower Q (more losses)
- **High n**: Higher Q (better confinement)
- **High ω**: Lower Q (more radiation losses)

This affects the decay time and stored energy.

### Radiation Backreaction

The EM field exerts pressure on the bubble:
- **High E_em**: Strong backreaction force
- **Small R**: Stronger force (F ~ 1/R)
- **Feedback loop**: Backreaction affects R → affects E_em → affects backreaction

### Bremsstrahlung + Recombination

Plasma processes add energy to EM modes:
- **Bremsstrahlung**: Electron-ion collisions → photons → stored energy
- **Recombination**: Electron-ion recombination → photons → stored energy
- **Energy source**: Additional energy beyond parametric pumping

## Summary

✅ **EM mode frequency shifting**: ω_k(R(t), n(r,t)) implemented
✅ **Mode squeezing equations**: ȧk = f(ak, R, Ṙ, ∂n/∂t) implemented
✅ **Stored-energy pump from gradients**: Explicit gradient-based pumping
✅ **Dynamic Q-factor & cavity losses**: Q(R, n, ω) implemented
✅ **Radiation backreaction term**: F = -E_em/R implemented
✅ **Bremsstrahlung + recombination emission**: Added to stored energy

The EM model is now a **formal cavity-QED subsystem** with comprehensive physics!

