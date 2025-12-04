// analysis/diagnostics.ts
// Diagnostic tools for monitoring bubble dynamics and detecting extreme conditions

import { BubbleFullState } from "../model/types";
import {
  computeGradients,
  computeSpatialGradients,
  estimateEmission,
  GradientMetrics,
} from "./observables";
import { Constants, Calculations } from "../core/units";

/**
 * Reaction enthalpies [J/mol] at standard conditions
 * Positive = endothermic (requires energy)
 * Negative = exothermic (releases energy)
 */
const REACTION_ENTHALPIES = {
  // Reaction 0: H2O ↔ H + OH
  H2O_dissociation: 500e3 * Constants.R_gas, // ~500 kJ/mol
  
  // Reaction 1: O2 ↔ 2O
  O2_dissociation: 498e3 * Constants.R_gas, // ~498 kJ/mol
  
  // Reaction 2: H + O ↔ OH (recombination, exothermic)
  H_O_recombination: -428e3 * Constants.R_gas, // ~-428 kJ/mol (releases energy)
} as const;

/**
 * Compute chemical energy change from reaction progress
 * 
 * Chemical energy change = Σ(ΔH_i · dξ_i/dt · dt)
 * where ΔH_i is the reaction enthalpy and dξ_i/dt is the reaction rate
 */
function computeChemicalEnergy(
  state: BubbleFullState,
  statePrev: BubbleFullState,
  dt: number
): number {
  const { reactions } = state;
  const { reactions: reactionsPrev } = statePrev;
  const { T } = state.gas;
  
  // Compute reaction progress changes
  const dxi0 = reactions.xi[0] - reactionsPrev.xi[0];
  const dxi1 = reactions.xi.length > 1 ? reactions.xi[1] - reactionsPrev.xi[1] : 0;
  const dxi2 = reactions.xi.length > 2 ? reactions.xi[2] - reactionsPrev.xi[2] : 0;
  
  // Convert reaction progress to number of moles (simplified: use volume)
  const { R } = state.hydro;
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  const volumePrev = (4.0 / 3.0) * Math.PI * statePrev.hydro.R * statePrev.hydro.R * statePrev.hydro.R;
  const avgVolume = (volume + volumePrev) / 2.0;
  
  // Approximate number of moles (using ideal gas law)
  const n_moles = (state.gas.Pg * avgVolume) / (Constants.R_gas * T);
  
  // Energy change from each reaction
  // Reaction 0: H2O ↔ H + OH (endothermic)
  const E0 = REACTION_ENTHALPIES.H2O_dissociation * dxi0 * n_moles;
  
  // Reaction 1: O2 ↔ 2O (endothermic)
  const E1 = REACTION_ENTHALPIES.O2_dissociation * dxi1 * n_moles;
  
  // Reaction 2: H + O ↔ OH (exothermic, releases energy)
  const E2 = REACTION_ENTHALPIES.H_O_recombination * dxi2 * n_moles;
  
  // Total chemical energy change (positive = energy absorbed, negative = energy released)
  return E0 + E1 + E2;
}

export interface EnergyBudget {
  acousticWorkIn: number; // Acoustic work input [J]
  compressionWork: number; // Compression work [J]
  heatLoss: number; // Heat loss to surroundings [J]
  lightEmission: number; // Light emission [J]
  viscousDissipation: number; // Viscous dissipation [J]
  chemicalEnergy: number; // Chemical reaction energy [J]
  emStoredEnergy: number; // EM stored energy [J]
  totalEnergy: number; // Total energy (conservation check)
}

export interface StabilityMetrics {
  isStable: boolean; // System appears stable
  limitCycleDetected: boolean; // Periodic orbit detected
  runawayHeating: boolean; // Temperature increasing without bound
  energyConserved: boolean; // Energy conservation check
  periodEstimate?: number; // Estimated oscillation period [s]
}

export interface RefractiveIndexEstimate {
  n_effective: number; // Effective refractive index
  n_plasma: number; // Plasma contribution
  n_neutral: number; // Neutral gas contribution
  n_gradient: number; // Radial gradient magnitude
}

export interface PlasmaDiagnostics {
  plasmaFrequency: number; // ωp [rad/s]
  debyeLength: number; // λD [m]
  opticalFrequency: number; // Typical optical frequency [rad/s]
  plasmaCutoff: boolean; // True if ωp > ω_optical (mode cutoff)
  meanFreePath: number; // Electron mean free path [m]
  collisional: boolean; // True if mean free path << bubble size
}

/**
 * Compute energy budget for the bubble system
 * 
 * Tracks energy flows to detect conservation and identify
 * dominant energy pathways.
 */
export function computeEnergyBudget(
  state: BubbleFullState,
  statePrev: BubbleFullState,
  dt: number,
  acousticPower: number
): EnergyBudget {
  const { R, Rdot } = state.hydro;
  const { Pg, T } = state.gas;
  const { ne, Te } = state.plasma;
  const { em } = state;

  // Acoustic work input (simplified: acoustic power * time)
  const acousticWorkIn = acousticPower * dt;

  // Compression work: P * dV
  // dV = 4*π*R²*dR = 4*π*R²*Rdot*dt
  const dV = 4.0 * Math.PI * R * R * Rdot * dt;
  const compressionWork = Pg * dV;

  // Heat loss (simplified: proportional to surface area and temperature difference)
  const T_ambient = 293.15; // 20°C
  const heatLossCoeff = 100; // W/(m²·K) - approximate
  const surfaceArea = 4.0 * Math.PI * R * R;
  const heatLoss = heatLossCoeff * (T - T_ambient) * surfaceArea * dt;

  // Light emission
  const emission = estimateEmission(state);
  const lightEmission = emission.totalPower * dt;

  // Viscous dissipation (simplified: ~mu * (Rdot/R)² * volume)
  const mu = 1e-3; // Water viscosity [Pa·s]
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  const strainRate = Math.abs(Rdot / Math.max(R, 1e-10));
  const viscousDissipation = mu * strainRate * strainRate * volume * dt;

  // Chemical energy: compute from reaction progress and enthalpies
  const chemicalEnergy = computeChemicalEnergy(state, statePrev, dt);

  // EM stored energy
  const emStoredEnergy = em.storedEnergy;

  // Total energy (should be approximately conserved over cycles)
  const totalEnergy =
    acousticWorkIn -
    compressionWork -
    heatLoss -
    lightEmission -
    viscousDissipation -
    chemicalEnergy -
    emStoredEnergy;

  return {
    acousticWorkIn,
    compressionWork,
    heatLoss,
    lightEmission,
    viscousDissipation,
    chemicalEnergy,
    emStoredEnergy,
    totalEnergy,
  };
}

/**
 * Detect extreme gradient events
 * 
 * Identifies when the system enters non-classical regimes
 * with extreme spatial or temporal gradients.
 */
export function detectExtremeGradients(
  statePrev: BubbleFullState,
  stateNext: BubbleFullState,
  dt: number
): {
  detected: boolean;
  metrics: GradientMetrics;
  spatialGradients: ReturnType<typeof computeSpatialGradients>;
} {
  const metrics = computeGradients(statePrev, stateNext, dt);
  const spatialGradients = computeSpatialGradients(stateNext);

  // Additional spatial gradient thresholds
  const extremePressureGradient = spatialGradients.pressureGradient > 1e18; // Pa/m
  const extremeDensityGradient = spatialGradients.densityGradient > 1e40; // 1/m⁴

  const detected = metrics.extremeGradient || extremePressureGradient || extremeDensityGradient;

  return {
    detected,
    metrics,
    spatialGradients,
  };
}

/**
 * Estimate refractive index from bubble state
 * 
 * Combines contributions from:
 * - Neutral gas (density-dependent)
 * - Plasma (frequency-dependent, can be negative)
 * - Ionization effects
 */
export function estimateRefractiveIndex(
  state: BubbleFullState,
  opticalFrequency: number
): RefractiveIndexEstimate {
  const { R } = state.hydro;
  const { T } = state.gas;
  const { ne } = state.plasma;
  const { numberDensity } = state.species;

  // Neutral gas refractive index (simplified: n ≈ 1 + α * ρ)
  // α is polarizability, ρ is density
  const totalNeutralDensity = Object.values(numberDensity).reduce(
    (sum, n) => sum + n,
    0
  );
  const alpha = 1e-29; // Approximate polarizability [m³]
  const n_neutral = 1.0 + alpha * totalNeutralDensity;

  // Plasma refractive index: n² = 1 - (ωp² / ω²)
  // For ωp > ω, n becomes imaginary (cutoff)
  const omega_p = Calculations.plasmaFrequency(ne);
  const n_plasma_squared = 1.0 - (omega_p * omega_p) / (opticalFrequency * opticalFrequency);
  const n_plasma = n_plasma_squared > 0 ? Math.sqrt(n_plasma_squared) : 0;

  // Effective refractive index (simplified combination)
  const n_effective = n_neutral + (n_plasma - 1.0) * 0.5; // Weighted average

  // Radial gradient (simplified: ~n / R)
  const R_safe = Math.max(R, 1e-10);
  const n_gradient = Math.abs(n_effective) / R_safe;

  return {
    n_effective,
    n_plasma,
    n_neutral,
    n_gradient,
  };
}

/**
 * Compute plasma diagnostics
 * 
 * Analyzes plasma properties and their relationship to
 * optical frequencies and bubble size.
 */
export function computePlasmaDiagnostics(
  state: BubbleFullState,
  opticalFrequency: number = 2 * Math.PI * 3e14 // ~500 nm
): PlasmaDiagnostics {
  const { R } = state.hydro;
  const { T } = state.gas;
  const { ne, Te } = state.plasma;

  const plasmaFrequency = Calculations.plasmaFrequency(ne);
  const debyeLength = Calculations.debyeLength(ne, T);

  // Electron mean free path (simplified)
  // λ_mfp ≈ v_th / ν_coll, where v_th = sqrt(k_B*Te/m_e)
  const v_th = Math.sqrt((Constants.k_B * Te) / Constants.m_e);
  const nu_coll = 1e13; // Approximate collision frequency [1/s]
  const meanFreePath = v_th / nu_coll;

  const plasmaCutoff = plasmaFrequency >= opticalFrequency;
  const collisional = meanFreePath < R * 0.1; // Mean free path << bubble size

  return {
    plasmaFrequency,
    debyeLength,
    opticalFrequency,
    plasmaCutoff,
    meanFreePath,
    collisional,
  };
}

/**
 * Compute stability metrics
 * 
 * Detects limit cycles, runaway heating, and energy conservation
 * to assess system stability.
 */
export function computeStabilityMetrics(
  states: BubbleFullState[],
  times: number[]
): StabilityMetrics {
  if (states.length < 10) {
    // Need enough data points
    return {
      isStable: false,
      limitCycleDetected: false,
      runawayHeating: false,
      energyConserved: false,
    };
  }

  // Check for runaway heating
  const temperatures = states.map((s) => s.gas.T);
  const maxTemp = Math.max(...temperatures);
  const minTemp = Math.min(...temperatures);
  const tempTrend = temperatures[temperatures.length - 1] - temperatures[0];
  const runawayHeating = tempTrend > 1000 && maxTemp > 1e5; // Rapid temperature increase

  // Check for limit cycle (periodic behavior)
  // Simplified: look for repeating patterns in radius
  const radii = states.map((s) => s.hydro.R);
  const maxR = Math.max(...radii);
  const minR = Math.min(...radii);
  const amplitude = maxR - minR;

  // Count zero crossings of (R - R_mean) to estimate period
  const R_mean = radii.reduce((sum, r) => sum + r, 0) / radii.length;
  let zeroCrossings = 0;
  for (let i = 1; i < radii.length; i++) {
    if (
      (radii[i - 1] - R_mean) * (radii[i] - R_mean) < 0
    ) {
      zeroCrossings++;
    }
  }

  const periodEstimate =
    zeroCrossings > 0
      ? (times[times.length - 1] - times[0]) / (zeroCrossings / 2)
      : undefined;

  const limitCycleDetected =
    zeroCrossings >= 4 && amplitude > 0.1 * R_mean; // At least 2 full cycles

  // Energy conservation check (simplified)
  // Energy should be approximately constant over a cycle
  const energies = states.map((s) => {
    const volume = (4.0 / 3.0) * Math.PI * s.hydro.R * s.hydro.R * s.hydro.R;
    return s.gas.Pg * volume + s.em.storedEnergy; // Simplified energy
  });
  const energyVariance =
    energies.reduce((sum, e) => {
      const mean = energies.reduce((s, v) => s + v, 0) / energies.length;
      return sum + (e - mean) * (e - mean);
    }, 0) / energies.length;
  const energyMean = energies.reduce((sum, e) => sum + e, 0) / energies.length;
  const energyConserved = energyVariance / (energyMean * energyMean) < 0.1; // <10% variance

  const isStable =
    !runawayHeating && limitCycleDetected && energyConserved;

  return {
    isStable,
    limitCycleDetected,
    runawayHeating,
    energyConserved,
    periodEstimate,
  };
}

/**
 * Compute mode squeezing metrics
 * 
 * Analyzes EM mode amplitudes to detect parametric amplification
 * and mode squeezing effects.
 */
export function computeModeSqueezingMetrics(
  state: BubbleFullState
): {
  modeAmplitudes: number[];
  totalModeEnergy: number;
  squeezingDetected: boolean;
  pumpEfficiency: number;
} {
  const { modes, storedEnergy } = state.em;

  const modeAmplitudes = modes.map(
    (m) => Math.sqrt(m.re * m.re + m.im * m.im)
  );
  const totalModeEnergy = modes.reduce(
    (sum, m) => sum + (m.re * m.re + m.im * m.im),
    0
  );

  // Squeezing detected if mode amplitudes are large relative to stored energy
  // (indicating parametric amplification)
  const squeezingDetected =
    totalModeEnergy > storedEnergy * 0.1 && storedEnergy > 1e-20;

  // Pump efficiency: ratio of stored energy to input (simplified)
  const pumpEfficiency = storedEnergy > 0 ? totalModeEnergy / storedEnergy : 0;

  return {
    modeAmplitudes,
    totalModeEnergy,
    squeezingDetected,
    pumpEfficiency,
  };
}

