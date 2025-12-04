# Logical Expansions and Future Enhancements

This document outlines planned logical expansions and enhancements to the sonoluminescence model.

## Current TODOs and Placeholders

### 1. Spectral Emission Calculation ⚠️
**Status**: Placeholder only  
**Location**: `analysis/observables.ts`  
**Current**: `spectrum?: number[]` is undefined

**Implementation Plan**:
- Compute wavelength-dependent emission spectrum
- Blackbody spectrum: Planck's law `B(λ, T) = (2πhc²/λ⁵) / (exp(hc/(λkT)) - 1)`
- Bremsstrahlung spectrum: Frequency-dependent power density
- EM mode contributions: Discrete spectral lines from cavity modes
- Combine into total spectrum with proper normalization

**Use Case**: Compare with experimental spectra, identify emission mechanisms

### 2. Chemical Energy Calculation ⚠️
**Status**: TODO - Currently returns 0  
**Location**: `analysis/diagnostics.ts:91`

**Implementation Plan**:
- Track reaction enthalpies for each reaction
- Compute energy change: `ΔE_chem = Σ(ΔH_i · dξ_i/dt)`
- Include dissociation energies (H₂O → H + OH: ~500 kJ/mol)
- Include recombination energy release
- Properly account for temperature-dependent enthalpies

**Use Case**: Complete energy budget, understand chemical contribution to heating

### 3. EM Mode Frequency Computation ⚠️
**Status**: Placeholder comment  
**Location**: `core/statevector.ts:170`

**Implementation Plan**:
- Compute mode frequencies from R(t) in `fromVector()`
- Use: `ω_k(t) = ω_k0 · (R₀/R(t)) · √(1 - (ω_p/ω_k0)²)`
- Update mode frequencies dynamically based on current state
- Consider plasma frequency cutoff effects

**Use Case**: Accurate EM mode tracking, proper cavity dynamics

## High-Priority Expansions

### 4. Initial State Helpers
**Priority**: High  
**Benefit**: Easier to use, less error-prone

**Features**:
- `createEquilibriumState()`: Bubble at acoustic equilibrium
- `createCollapseState()`: State just before collapse
- `createExpandedState()`: State at maximum radius
- `createFromPreset()`: Generate state matching preset parameters
- Validation of physical consistency

**Implementation**:
```typescript
// analysis/initialStates.ts
export function createEquilibriumState(
  params: SonoluminescenceParams,
  R_equilibrium: number
): BubbleFullState {
  // Compute equilibrium state from parameters
  // - Hydro: R = R_equilibrium, Rdot = 0
  // - Gas: Pg from acoustic balance
  // - Species: From preset or defaults
  // - Plasma: Low ionization at equilibrium
  // - EM: Ground state modes
  // - Acoustic: Phase = 0
}
```

### 5. Simulation Runner API
**Priority**: High  
**Benefit**: High-level interface, easier workflow

**Features**:
- `runSimulation()`: Complete simulation with analysis
- Automatic result processing
- Event detection configuration
- Progress callbacks
- Result aggregation

**Implementation**:
```typescript
// simulation/runner.ts
export interface SimulationConfig {
  model: SonoluminescenceModel;
  initialState: BubbleFullState;
  integratorOptions: IntegratorOptions;
  events?: EventFunction[];
  analysis?: {
    computeEmission?: boolean;
    computeGradients?: boolean;
    computeDiagnostics?: boolean;
  };
}

export interface SimulationResult {
  timeSeries: IntegrationResult;
  analysis?: {
    emissions?: EmissionSnapshot[];
    gradients?: GradientMetrics[];
    diagnostics?: any;
  };
  events?: Array<{ time: number; state: BubbleFullState; eventIndex: number }>;
}
```

### 6. Export Utilities
**Priority**: Medium  
**Benefit**: Data analysis, visualization, sharing

**Features**:
- CSV export: Time series data
- JSON export: Complete state history
- HDF5 export: Large datasets (future)
- Selective export: Choose which variables
- Compression: For large datasets

**Implementation**:
```typescript
// io/export.ts
export function exportToCSV(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  filename: string,
  variables?: DimensionId[]
): void;

export function exportToJSON(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  filename: string
): void;
```

## Medium-Priority Expansions

### 7. Expanded Reaction Set
**Priority**: Medium  
**Benefit**: More realistic chemistry

**Additional Reactions**:
- OH dissociation: OH ↔ O + H
- Nitrogen reactions: N₂ ↔ 2N
- Radical recombination: H + OH → H₂O
- Three-body recombination: H + H + M → H₂ + M
- More complex reaction networks

**Implementation**: Extend `reactions.ts` with more reaction progress variables and rate constants

### 8. Validation and Benchmarking
**Priority**: Medium  
**Benefit**: Confidence in model accuracy

**Tests**:
- Adiabatic compression: Compare with analytical solution
- Saha equilibrium: Validate ionization at known T, n
- Rayleigh-Plesset: Compare with known solutions
- Energy conservation: Check over full cycle
- Experimental comparison: Match published data

**Implementation**:
```typescript
// validation/tests.ts
export function testAdiabaticCompression(): boolean;
export function testSahaEquilibrium(): boolean;
export function testEnergyConservation(): boolean;
export function benchmarkAgainstExperiment(): ComparisonResult;
```

### 9. More EM Modes
**Priority**: Low  
**Benefit**: Better spectral resolution

**Expansion**:
- Increase from 3 to 10+ modes
- Support configurable number of modes
- Higher frequency modes for UV emission

**Implementation**: Make mode count configurable in `EmCavityParams`

## Advanced Expansions (Future)

### 10. Spatial Resolution
**Priority**: Low (Major undertaking)  
**Benefit**: Radial gradients, non-uniform conditions

**Features**:
- 1D radial model: R(r, t) with radial dependence
- 2D axisymmetric: Full spatial model
- Non-uniform temperature/density
- Radial pressure gradients

**Challenges**: Requires PDE solver, much more complex

### 11. Quantum Corrections
**Priority**: Low  
**Benefit**: More accurate at extreme conditions

**Features**:
- Quantum pressure corrections
- Quantum tunneling effects
- Zero-point energy
- Quantum statistics for electrons

### 12. Radiation Transport
**Priority**: Low  
**Benefit**: Accurate spectral emission

**Features**:
- Radiative transfer equation
- Absorption/emission coefficients
- Scattering effects
- Detailed spectral lines

### 13. Multi-Bubble Interactions
**Priority**: Low  
**Benefit**: Realistic experimental conditions

**Features**:
- Multiple bubbles
- Bubble-bubble interactions
- Acoustic field coupling
- Collective effects

## Implementation Priority

### Phase 1 (Immediate)
1. ✅ Spectral emission calculation
2. ✅ Chemical energy calculation
3. ✅ EM mode frequency computation
4. ✅ Initial state helpers

### Phase 2 (Short-term)
5. ✅ Simulation runner API
6. ✅ Export utilities
7. ✅ Validation tests

### Phase 3 (Medium-term)
8. Expanded reaction set
9. More EM modes (configurable)

### Phase 4 (Long-term)
10. Spatial resolution
11. Quantum corrections
12. Radiation transport
13. Multi-bubble interactions

## Notes

- All expansions maintain backward compatibility
- New features are optional (don't break existing code)
- Performance is considered (e.g., spectral calculation only if requested)
- Type safety is maintained throughout

