// physics/thermoChem.ts
import { BubbleFullState, SpeciesId } from "../model/types";
import { Constants } from "../core/units";

export interface ThermoParams {
  gamma: number;  // effective polytropic exponent (Cp/Cv)
  R0: number;     // reference radius [m]
  Pg0: number;    // reference gas pressure [Pa]
  T0: number;     // reference temperature [K]
  heatLossCoeff: number; // heat loss coefficient [W/(m²·K)] (legacy, use detailedHeatTransfer)
  
  // Non-ideal gas effects
  useNonIdealGas?: boolean; // Use van der Waals EOS
  vanDerWaals_a?: number;   // van der Waals parameter a [Pa·m⁶/mol²]
  vanDerWaals_b?: number;    // van der Waals parameter b [m³/mol]
  
  // Temperature-dependent properties
  useTemperatureDependentGamma?: boolean; // Gamma varies with temperature
  gammaFunction?: (T: number) => number;   // Custom gamma(T) function
  
  // Energy relaxation times [s]
  tau_trans?: number;  // Translational relaxation (default: 1e-12)
  tau_rot?: number;    // Rotational relaxation (default: 1e-11)
  tau_vib?: number;     // Vibrational relaxation (default: 1e-9)
  tau_elec?: number;   // Electronic relaxation (default: 1e-6)
  
  // Enhanced energy exchange
  useDetailedEnergyExchange?: boolean; // Use temperature-dependent relaxation
  useLandauTellerRelaxation?: boolean; // Use Landau-Teller model for vibrational relaxation
  
  // Detailed heat transfer
  useDetailedHeatTransfer?: boolean; // Use detailed heat transfer model
  thermalConductivity?: number; // Gas thermal conductivity [W/(m·K)]
  liquidThermalConductivity?: number; // Liquid thermal conductivity [W/(m·K)]
  convectiveCoeff?: number; // Convective heat transfer coefficient [W/(m²·K)]
  includeRadiation?: boolean; // Include radiative heat transfer
  
  // Detailed heat capacity
  useDetailedHeatCapacity?: boolean; // Use temperature-dependent heat capacity
  
  // Species diffusion
  useSpeciesDiffusion?: boolean; // Include diffusion through bubble wall
  diffusionCoefficients?: Record<string, number>; // Diffusion coefficients [m²/s] for each species
  liquidConcentrations?: Record<string, number>; // Liquid-side concentrations [mol/m³]
  
  // Van der Waals iterative solution
  useIterativeVanDerWaals?: boolean; // Use iterative solution instead of approximation
  vdwIterationTolerance?: number; // Convergence tolerance (default: 1e-6)
  vdwMaxIterations?: number; // Maximum iterations (default: 100)
}

export interface ThermoDerivatives {
  dPgdt: number;
  dTdt: number;
}

export interface SpeciesDerivatives {
  dNdt: Record<SpeciesId, number>;
}

export interface InternalEnergyDerivatives {
  dE_trans_dt: number;
  dE_rot_dt: number;
  dE_vib_dt: number;
  dE_elec_dt: number;
}

/**
 * Van der Waals equation of state
 * (P + a*n²/V²) * (V - n*b) = n*R_gas*T
 * 
 * Solving for P: P = n*R_gas*T/(V - n*b) - a*n²/V²
 */
function vanDerWaalsPressure(
  n: number,  // number of moles
  V: number,  // volume
  T: number,  // temperature
  a: number,   // van der Waals parameter a
  b: number    // van der Waals parameter b
): number {
  const V_corrected = V - n * b;
  if (V_corrected <= 0) {
    // Volume too small, use ideal gas as fallback
    return (n * Constants.R_gas * T) / V;
  }
  const ideal_term = (n * Constants.R_gas * T) / V_corrected;
  const correction_term = (a * n * n) / (V * V);
  return ideal_term - correction_term;
}

/**
 * Iterative van der Waals solution
 * 
 * For given P, T, n, solve for V using Newton-Raphson iteration
 * 
 * Van der Waals equation: (P + a*n²/V²) * (V - n*b) = n*R_gas*T
 * 
 * Rearranged: f(V) = P*V³ - (P*n*b + n*R_gas*T)*V² + a*n²*V - a*n³*b = 0
 * 
 * Or: f(V) = V³ - (n*b + n*R_gas*T/P)*V² + (a*n²/P)*V - (a*n³*b/P) = 0
 * 
 * Newton-Raphson: V_new = V_old - f(V_old) / f'(V_old)
 */
function solveVanDerWaalsVolume(
  P: number,  // Pressure [Pa]
  T: number,  // Temperature [K]
  n: number,  // Number of moles
  a: number,   // van der Waals parameter a
  b: number,   // van der Waals parameter b
  tolerance: number = 1e-6,
  maxIterations: number = 100
): number {
  // Initial guess: ideal gas volume
  let V = (n * Constants.R_gas * T) / P;
  const V_min = n * b * 1.001; // Minimum volume (just above n*b)
  
  for (let iter = 0; iter < maxIterations; iter++) {
    const V_corrected = V - n * b;
    if (V_corrected <= 0) {
      V = V_min;
      continue;
    }
    
    // f(V) = (P + a*n²/V²) * (V - n*b) - n*R_gas*T
    const P_corr = P + (a * n * n) / (V * V);
    const f = P_corr * V_corrected - n * Constants.R_gas * T;
    
    // f'(V) = d/dV[f(V)]
    // f'(V) = P + a*n²/V² - 2*a*n²*(V - n*b)/V³ - (P + a*n²/V²)
    // Simplified: f'(V) = -2*a*n²*(V - n*b)/V³ + a*n²/V²
    const df_dV = -2.0 * a * n * n * V_corrected / (V * V * V) + a * n * n / (V * V) + P;
    
    if (Math.abs(df_dV) < 1e-10) {
      // Derivative too small, use fallback
      break;
    }
    
    const V_new = V - f / df_dV;
    
    // Check convergence
    if (Math.abs(V_new - V) < tolerance * V) {
      return Math.max(V_new, V_min);
    }
    
    V = Math.max(V_new, V_min);
  }
  
  return V;
}

/**
 * Compute van der Waals pressure derivative using iterative solution
 * 
 * For van der Waals: P = n*R_gas*T/(V - n*b) - a*n²/V²
 * 
 * dP/dt = (∂P/∂V) * (dV/dt) + (∂P/∂T) * (dT/dt)
 * 
 * Where:
 * ∂P/∂V = -n*R_gas*T/(V - n*b)² + 2*a*n²/V³
 * ∂P/∂T = n*R_gas/(V - n*b)
 */
function computeVanDerWaalsDerivative(
  state: BubbleFullState,
  params: ThermoParams,
  volumeChangeRate: number,
  dTdt: number
): number {
  const { R } = state.hydro;
  const { Pg, T } = state.gas;
  const { vanDerWaals_a = 0.136, vanDerWaals_b = 3.2e-5 } = params;
  
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  const n_total = Object.values(state.species.numberDensity).reduce(
    (sum, n) => sum + n,
    0
  ) * volume;
  const n_moles = n_total / Constants.N_A;
  
  const V_corrected = volume - n_moles * vanDerWaals_b;
  if (V_corrected <= 0) {
    // Fallback to ideal gas
    return -params.gamma * Pg * volumeChangeRate;
  }
  
  // Partial derivatives
  const dP_dV = -(n_moles * Constants.R_gas * T) / (V_corrected * V_corrected) 
                + (2.0 * vanDerWaals_a * n_moles * n_moles) / (volume * volume * volume);
  const dP_dT = (n_moles * Constants.R_gas) / V_corrected;
  
  // dV/dt = volume * volumeChangeRate
  const dV_dt = volume * volumeChangeRate;
  
  // Total derivative
  const dPgdt = dP_dV * dV_dt + dP_dT * dTdt;
  
  return dPgdt;
}

/**
 * Compute heat capacity at constant pressure Cp for a species
 * 
 * Uses statistical mechanics:
 * - Translational: (3/2) * R
 * - Rotational: R per rotational mode (diatomic: 2 modes, linear: 2, nonlinear: 3)
 * - Vibrational: R * sum over modes of (θ_v/T)² * exp(θ_v/T) / (exp(θ_v/T) - 1)²
 * 
 * Simplified: Cp = Cp_trans + Cp_rot + Cp_vib
 */
function computeHeatCapacityCp(
  species: SpeciesId,
  T: number
): number {
  // Base translational contribution (always 3/2 * R)
  const Cp_trans = (3.0 / 2.0) * Constants.R_gas;
  
  // Rotational contribution
  let Cp_rot = 0;
  // Monatomic: no rotation
  // Diatomic/linear: 2 rotational modes
  // Nonlinear: 3 rotational modes
  if (species === 'Ar' || species === 'Xe' || species === 'H' || species === 'O' || species === 'N') {
    // Monatomic: no rotation
    Cp_rot = 0;
  } else if (species === 'H2O') {
    // Nonlinear triatomic: 3 rotational modes
    Cp_rot = 3.0 * Constants.R_gas;
  } else {
    // Diatomic (O2, N2, OH): 2 rotational modes
    Cp_rot = 2.0 * Constants.R_gas;
  }
  
  // Vibrational contribution
  let Cp_vib = 0;
  if (species === 'H2O') {
    // H2O has 3 vibrational modes
    const theta_v1 = 2290; // Symmetric stretch [K]
    const theta_v2 = 5250; // Asymmetric stretch [K]
    const theta_v3 = 5400; // Bending [K]
    const modes = [theta_v1, theta_v2, theta_v3];
    for (const theta_v of modes) {
      const x = theta_v / T;
      if (x > 0.01) {
        const exp_x = Math.exp(x);
        Cp_vib += Constants.R_gas * x * x * exp_x / ((exp_x - 1) * (exp_x - 1));
      }
    }
  } else if (species === 'O2' || species === 'N2') {
    // Diatomic: 1 vibrational mode
    const theta_v = species === 'O2' ? 2270 : 3390; // [K]
    const x = theta_v / T;
    if (x > 0.01) {
      const exp_x = Math.exp(x);
      Cp_vib = Constants.R_gas * x * x * exp_x / ((exp_x - 1) * (exp_x - 1));
    }
  } else if (species === 'OH') {
    // Diatomic: 1 vibrational mode
    const theta_v = 5360; // [K]
    const x = theta_v / T;
    if (x > 0.01) {
      const exp_x = Math.exp(x);
      Cp_vib = Constants.R_gas * x * x * exp_x / ((exp_x - 1) * (exp_x - 1));
    }
  }
  // Monatomic: no vibration
  
  return Cp_trans + Cp_rot + Cp_vib; // J/(mol·K)
}

/**
 * Compute heat capacity at constant volume Cv
 * 
 * Cv = Cp - R (for ideal gas)
 */
function computeHeatCapacityCv(
  species: SpeciesId,
  T: number
): number {
  const Cp = computeHeatCapacityCp(species, T);
  return Cp - Constants.R_gas; // J/(mol·K)
}

/**
 * Compute mixture heat capacity (weighted average)
 */
function computeMixtureHeatCapacity(
  state: BubbleFullState,
  T: number,
  useDetailed: boolean
): { Cp: number; Cv: number } {
  if (!useDetailed) {
    // Simplified: average diatomic
    const Cp = (5.0 / 2.0) * Constants.R_gas;
    const Cv = Cp - Constants.R_gas;
    return { Cp, Cv };
  }
  
  // Detailed: weighted average over species
  const volume = (4.0 / 3.0) * Math.PI * Math.pow(state.hydro.R, 3);
  const species = state.species.numberDensity;
  
  let totalCp = 0;
  let totalCv = 0;
  let totalMoles = 0;
  
  for (const [speciesId, n] of Object.entries(species) as [SpeciesId, number][]) {
    const n_moles = (n * volume) / Constants.N_A;
    if (n_moles > 0) {
      const Cp_species = computeHeatCapacityCp(speciesId, T);
      const Cv_species = computeHeatCapacityCv(speciesId, T);
      totalCp += n_moles * Cp_species;
      totalCv += n_moles * Cv_species;
      totalMoles += n_moles;
    }
  }
  
  if (totalMoles > 0) {
    return {
      Cp: totalCp / totalMoles, // J/(mol·K)
      Cv: totalCv / totalMoles,
    };
  } else {
    // Fallback
    const Cp = (5.0 / 2.0) * Constants.R_gas;
    const Cv = Cp - Constants.R_gas;
    return { Cp, Cv };
  }
}

/**
 * Compute effective gamma (Cp/Cv) as function of temperature
 * 
 * For diatomic gases: gamma decreases with temperature as vibrational modes activate
 * gamma(T) ≈ 1 + 2/(f(T)) where f(T) is effective degrees of freedom
 */
function computeTemperatureDependentGamma(T: number, gamma0: number = 1.4): number {
  // Simplified model: gamma decreases as vibrational modes activate
  // At low T: gamma ≈ 1.4 (diatomic, rotation only)
  // At high T: gamma → 1.33 (vibrational modes active)
  const T_vib = 3000; // Characteristic vibrational temperature [K]
  const gamma_high = 1.33;
  const transition = Math.exp(-T_vib / T);
  return gamma0 * (1 - transition) + gamma_high * transition;
}

/**
 * Compute thermodynamic derivatives (pressure and temperature)
 * 
 * For adiabatic compression/expansion:
 * - dPg/dt from volume change: dPg/dt = -gamma * Pg * (3*Rdot/R)
 * - dT/dt from adiabatic relation: dT/dt = (gamma-1) * T * (3*Rdot/R) + heat_loss
 * 
 * Enhanced with:
 * - Non-ideal gas effects (van der Waals)
 * - Temperature-dependent gamma
 */
export function computeThermoDerivatives(
  state: BubbleFullState,
  params: ThermoParams
): ThermoDerivatives {
  const { R, Rdot } = state.hydro;
  const { Pg, T } = state.gas;
  const {
    gamma: gamma0,
    heatLossCoeff,
    useNonIdealGas = false,
    vanDerWaals_a = 0.136, // Argon: 0.136 Pa·m⁶/mol²
    vanDerWaals_b = 3.2e-5, // Argon: 3.2×10⁻⁵ m³/mol
    useTemperatureDependentGamma = false,
    gammaFunction,
  } = params;

  // Avoid division by zero
  const R_safe = Math.max(R, 1e-10);
  const volume = (4.0 / 3.0) * Math.PI * R_safe * R_safe * R_safe;

  // Compute effective gamma
  let gamma = gamma0;
  if (useTemperatureDependentGamma) {
    if (gammaFunction) {
      gamma = gammaFunction(T);
    } else {
      gamma = computeTemperatureDependentGamma(T, gamma0);
    }
  }

  // Volume change rate: dV/dt = 4*π*R²*Rdot
  // Relative volume change: (1/V) * dV/dt = 3*Rdot/R
  const volumeChangeRate = 3.0 * Rdot / R_safe;

  // Pressure change from adiabatic compression/expansion
  if (useNonIdealGas && vanDerWaals_a > 0 && vanDerWaals_b > 0) {
    // Temperature change (needed for iterative solution)
    const adiabaticHeating = (gamma - 1.0) * T * volumeChangeRate;
    const heatLossRate = computeHeatLossRate(state, params, R_safe);
    const dTdt = adiabaticHeating + heatLossRate;
    
    let dPgdt: number;
    
    if (params.useIterativeVanDerWaals) {
      // Iterative van der Waals solution
      // Use proper partial derivatives
      dPgdt = computeVanDerWaalsDerivative(
        state,
        params,
        volumeChangeRate,
        dTdt
      );
    } else {
      // Approximate solution (original)
      const dPgdt_ideal = -gamma * Pg * volumeChangeRate;
      
      // Correction term from van der Waals (simplified)
      // dP/dt correction ≈ -2*a*n²/V³ * dV/dt
      const n_total = Object.values(state.species.numberDensity).reduce(
        (sum, n) => sum + n,
        0
      ) * volume;
      const n_moles = n_total / Constants.N_A;
      const correction = -(2.0 * vanDerWaals_a * n_moles * n_moles) / (volume * volume * volume) * (4.0 * Math.PI * R_safe * R_safe * Rdot);
      
      dPgdt = dPgdt_ideal + correction;
    }
    
    return { dPgdt, dTdt };
  } else {
    // Ideal gas (original implementation)
    const dPgdt = -gamma * Pg * volumeChangeRate;
    
    const adiabaticHeating = (gamma - 1.0) * T * volumeChangeRate;
    const heatLossRate = computeHeatLossRate(state, params, R_safe);
    const dTdt = adiabaticHeating + heatLossRate;
    
    return { dPgdt, dTdt };
  }
}

/**
 * Compute detailed heat loss rate
 * 
 * Includes:
 * - Conduction through gas and liquid
 * - Convection at bubble wall
 * - Radiation (if enabled)
 */
function computeHeatLossRate(
  state: BubbleFullState,
  params: ThermoParams,
  R: number
): number {
  const { T } = state.gas;
  const T_ambient = params.T0;
  const surfaceArea = 4.0 * Math.PI * R * R;
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  
  if (params.useDetailedHeatTransfer) {
    // Detailed heat transfer model
    
    // 1. Conduction through gas
    // Heat flux: q = -k * dT/dr
    // Simplified: dT/dr ≈ (T - T_wall) / R
    // where T_wall is wall temperature (intermediate between T and T_ambient)
    const k_gas = params.thermalConductivity || 0.02; // W/(m·K) for air at room temp
    const T_wall = (T + T_ambient) / 2; // Simplified: linear temperature profile
    const conductionGas = -k_gas * (T - T_wall) / R * surfaceArea;
    
    // 2. Conduction through liquid boundary layer
    // Boundary layer thickness: δ ~ sqrt(alpha * t) where alpha is thermal diffusivity
    // Simplified: use characteristic boundary layer thickness
    const k_liquid = params.liquidThermalConductivity || 0.6; // W/(m·K) for water
    const delta_boundary = 1e-6; // ~1 micron boundary layer (simplified)
    const conductionLiquid = -k_liquid * (T_wall - T_ambient) / delta_boundary * surfaceArea;
    
    // 3. Convection
    const h_conv = params.convectiveCoeff || 1000; // W/(m²·K) - high for small bubbles
    const convection = -h_conv * (T_wall - T_ambient) * surfaceArea;
    
    // 4. Radiation (if enabled)
    let radiation = 0;
    if (params.includeRadiation) {
      // Stefan-Boltzmann: q = ε * σ * (T⁴ - T_ambient⁴)
      const epsilon = 0.9; // Emissivity (simplified)
      const sigma_SB = Constants.sigma_SB;
      radiation = -epsilon * sigma_SB * (T * T * T * T - T_ambient * T_ambient * T_ambient * T_ambient) * surfaceArea;
    }
    
    // Total heat loss rate
    // Convert to dT/dt: dT/dt = Q_loss / (m * Cp)
    // where m is mass and Cp is heat capacity
    const n_total = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0) * volume;
    const n_moles = n_total / Constants.N_A;
    
    // Heat capacity (detailed or simplified)
    const { Cp } = computeMixtureHeatCapacity(state, T, params.useDetailedHeatCapacity || false);
    const heatCapacity = n_moles * Cp; // J/K
    
    const totalHeatLoss = conductionGas + conductionLiquid + convection + radiation;
    const heatLossRate = totalHeatLoss / Math.max(heatCapacity, Constants.k_B * 1e20);
    
    return heatLossRate;
  } else {
    // Simplified heat loss (original)
    const heatLossRate = -params.heatLossCoeff * (T - T_ambient) * surfaceArea / (Constants.k_B * 1e20);
    return heatLossRate;
  }
}

/**
 * Compute species number density derivatives
 * 
 * Species evolve due to:
 * 1. Compression/expansion (conservation of mass)
 * 2. Chemical reactions (dissociation, recombination)
 * 3. Diffusion through bubble wall (if enabled)
 */
export function computeSpeciesDerivatives(
  state: BubbleFullState,
  params?: ThermoParams
): SpeciesDerivatives {
  const { R, Rdot } = state.hydro;
  const R_safe = Math.max(R, 1e-10);
  const volume = (4.0 / 3.0) * Math.PI * R_safe * R_safe * R_safe;
  const surfaceArea = 4.0 * Math.PI * R_safe * R_safe;
  
  // Relative volume change rate
  const volumeChangeRate = 3.0 * Rdot / R_safe;

  // For each species, compression/expansion changes number density
  // dN/dt = -N * (1/V) * dV/dt = -N * volumeChangeRate
  // (plus reaction terms which will be added in reactions module)
  
  const dNdt: Record<SpeciesId, number> = {
    H2O: -state.species.numberDensity.H2O * volumeChangeRate,
    O2: -state.species.numberDensity.O2 * volumeChangeRate,
    N2: -state.species.numberDensity.N2 * volumeChangeRate,
    Ar: -state.species.numberDensity.Ar * volumeChangeRate,
    Xe: -state.species.numberDensity.Xe * volumeChangeRate,
    H: -state.species.numberDensity.H * volumeChangeRate,
    O: -state.species.numberDensity.O * volumeChangeRate,
    OH: -(state.species.numberDensity.OH || 0) * volumeChangeRate,
    N: -(state.species.numberDensity.N || 0) * volumeChangeRate,
  };

  // Add diffusion terms if enabled
  if (params?.useSpeciesDiffusion && params.diffusionCoefficients) {
    const diffusionCoeffs = params.diffusionCoefficients;
    const liquidConcs = params.liquidConcentrations || {};
    
    for (const species of Object.keys(dNdt) as SpeciesId[]) {
      const D = diffusionCoeffs[species];
      if (D && D > 0) {
        // Fick's law: flux = -D * dC/dr
        // Simplified: dC/dr ≈ (C_liquid - C_gas) / R
        // where C_gas is gas-phase concentration
        const n_gas = state.species.numberDensity[species] || 0;
        const C_gas = n_gas; // Number density [1/m³]
        const C_liquid = (liquidConcs[species] || 0) * Constants.N_A; // Convert to number density
        
        // Diffusion flux through surface
        // J = -D * (C_liquid - C_gas) / R * surfaceArea
        // dN/dt from diffusion = J / volume
        const concentrationGradient = (C_liquid - C_gas) / R_safe;
        const flux = -D * concentrationGradient * surfaceArea;
        const diffusionRate = flux / Math.max(volume, 1e-30); // Rate of change of number density
        
        dNdt[species] += diffusionRate;
      }
    }
  }

  // Note: Chemical reaction terms will be added in the reactions module
  // and combined here or in the model integration step

  return { dNdt };
}

/**
 * Compute temperature-dependent relaxation times
 * 
 * Relaxation times typically decrease with temperature (faster collisions)
 * tau(T) = tau0 * (T0/T)^alpha
 */
function computeRelaxationTime(
  tau0: number,
  T: number,
  T_ref: number = 1000,
  alpha: number = 0.5
): number {
  if (T <= 0) return tau0;
  return tau0 * Math.pow(T_ref / T, alpha);
}

/**
 * Landau-Teller relaxation time for vibrational energy
 * 
 * The Landau-Teller model gives:
 * tau_vib(T) = tau0 * exp((T_vib/T)^(1/3))
 * 
 * where T_vib is the characteristic vibrational temperature
 * This is more accurate than simple power-law for vibrational relaxation
 */
function computeLandauTellerRelaxationTime(
  tau0: number,
  T: number,
  T_vib: number = 3000
): number {
  if (T <= 0) return tau0;
  const exponent = Math.pow(T_vib / T, 1.0 / 3.0);
  return tau0 * Math.exp(exponent);
}

/**
 * Compute internal energy partition derivatives
 * 
 * Tracks energy in different modes:
 * - Translational: (3/2) * n * k_B * T (for monatomic)
 * - Rotational: n * k_B * T (for diatomic, per rotational mode)
 * - Vibrational: quantum harmonic oscillator energy
 * - Electronic: excited state populations
 * 
 * Energy exchange between modes via relaxation processes
 * Enhanced with temperature-dependent relaxation times
 */
export function computeInternalEnergyDerivatives(
  state: BubbleFullState,
  params: ThermoParams
): InternalEnergyDerivatives {
  const { R, Rdot } = state.hydro;
  const { T } = state.gas;
  const R_safe = Math.max(R, 1e-10);
  
  // Total number density (sum of all species)
  const totalN = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0);
  const volume = (4.0 / 3.0) * Math.PI * R_safe * R_safe * R_safe;
  const n_total = totalN * volume;
  
  // Volume change rate
  const volumeChangeRate = 3.0 * Rdot / R_safe;

  // Relaxation times (configurable or use defaults)
  const tau_trans_base = params.tau_trans || 1e-12;
  const tau_rot_base = params.tau_rot || 1e-11;
  const tau_vib_base = params.tau_vib || 1e-9;
  const tau_elec_base = params.tau_elec || 1e-6;

  // Use temperature-dependent relaxation if enabled
  let tau_trans = tau_trans_base;
  let tau_rot = tau_rot_base;
  let tau_vib = tau_vib_base;
  let tau_elec = tau_elec_base;

  if (params.useDetailedEnergyExchange) {
    // Temperature-dependent relaxation: faster at higher T
    tau_trans = computeRelaxationTime(tau_trans_base, T, 1000, 0.3);
    tau_rot = computeRelaxationTime(tau_rot_base, T, 1000, 0.5);
    
    // Vibrational relaxation: use Landau-Teller if enabled
    if (params.useLandauTellerRelaxation) {
      const T_vib = 3000; // Characteristic vibrational temperature [K]
      tau_vib = computeLandauTellerRelaxationTime(tau_vib_base, T, T_vib);
    } else {
      tau_vib = computeRelaxationTime(tau_vib_base, T, 2000, 0.7); // Power-law fallback
    }
    
    tau_elec = computeRelaxationTime(tau_elec_base, T, 5000, 1.0); // Electronic very slow
  }

  // Equilibrium energy per mode
  // Translational: (3/2) * k_B * T per particle
  const E_trans_eq = (3.0 / 2.0) * Constants.k_B * T * n_total;
  
  // Rotational: k_B * T per rotational mode (assuming diatomic, 2 modes)
  const E_rot_eq = 2.0 * Constants.k_B * T * n_total;
  
  // Vibrational: quantum harmonic oscillator
  // At high T: E_vib ≈ k_B * T per mode
  // At low T: E_vib ≈ 0 (ground state)
  const T_vib = 3000; // Characteristic vibrational temperature [K]
  const vib_factor = T > T_vib ? 1.0 : Math.exp(-T_vib / T);
  const E_vib_eq = Constants.k_B * T * n_total * vib_factor;
  
  // Electronic: excited states (Boltzmann distribution)
  // E_elec ≈ k_B * T * exp(-E_excited / (k_B * T))
  const E_excited = 10 * Constants.e; // ~10 eV excited state
  const elec_factor = Math.exp(-E_excited / (Constants.k_B * T));
  const E_elec_eq = Constants.k_B * T * n_total * elec_factor * 0.1;

  // Current energy values
  const E_trans = state.internalEnergy.translational;
  const E_rot = state.internalEnergy.rotational;
  const E_vib = state.internalEnergy.vibrational;
  const E_elec = state.internalEnergy.electronic;

  // Energy change from compression/expansion (adiabatic work)
  // Plus relaxation toward equilibrium
  const dE_trans_dt = 
    -E_trans * volumeChangeRate + // compression work
    (E_trans_eq - E_trans) / tau_trans; // relaxation

  const dE_rot_dt = 
    -E_rot * volumeChangeRate +
    (E_rot_eq - E_rot) / tau_rot;

  const dE_vib_dt = 
    -E_vib * volumeChangeRate +
    (E_vib_eq - E_vib) / tau_vib;

  const dE_elec_dt = 
    -E_elec * volumeChangeRate +
    (E_elec_eq - E_elec) / tau_elec;

  return {
    dE_trans_dt,
    dE_rot_dt,
    dE_vib_dt,
    dE_elec_dt,
  };
}
