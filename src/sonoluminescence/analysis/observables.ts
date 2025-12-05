// analysis/observables.ts
import { BubbleFullState } from "../model/types";
import { Constants, Calculations } from "../core/units";

export interface EmissionSnapshot {
  totalPower: number; // Total emitted power [W]
  blackbodyPower: number; // Blackbody component [W]
  bremsstrahlungPower: number; // Bremsstrahlung component [W]
  emDecayPower: number; // EM stored energy decay [W]
  spectrum?: SpectralDistribution; // Wavelength-dependent spectrum
}

export interface SpectralDistribution {
  wavelengths: number[]; // Wavelengths [m]
  powerDensity: number[]; // Power per unit wavelength [W/m]
  blackbody: number[]; // Blackbody contribution [W/m]
  bremsstrahlung: number[]; // Bremsstrahlung contribution [W/m]
  emModes: number[]; // EM mode contributions [W/m]
}

export interface GradientMetrics {
  dRdt: number; // |∂R/∂t| [m/s]
  dRdotdt: number; // |∂Rdot/∂t| [m/s²]
  dPgdt: number; // |∂Pg/∂t| [Pa/s]
  dTdt: number; // |∂T/∂t| [K/s]
  dnedt: number; // |∂ne/∂t| [1/(m³·s)]
  maxGradient: number; // Maximum gradient magnitude
  extremeGradient: boolean; // True if gradients exceed thresholds
}

/**
 * Compute Planck's law for blackbody spectral radiance
 * B(λ, T) = (2πhc²/λ⁵) / (exp(hc/(λkT)) - 1) [W/(m³·sr)]
 * 
 * Integrated over solid angle and converted to power per unit wavelength
 */
function planckSpectrum(wavelength: number, T: number): number {
  if (wavelength <= 0 || T <= 0) return 0;
  
  const hc_over_lambda_kT = (Constants.h * Constants.c) / (wavelength * Constants.k_B * T);
  if (hc_over_lambda_kT > 100) return 0; // Avoid overflow
  
  const exp_term = Math.exp(hc_over_lambda_kT);
  if (exp_term === Infinity) return 0;
  
  const numerator = 2.0 * Math.PI * Constants.h * Constants.c * Constants.c;
  const denominator = wavelength * wavelength * wavelength * wavelength * wavelength * (exp_term - 1.0);
  
  return numerator / denominator; // [W/(m³·sr)]
}

/**
 * Compute bremsstrahlung spectral power density
 * Simplified: P(λ) ~ ne² * sqrt(Te) * exp(-hc/(λ*k_B*Te)) / λ²
 */
function bremsstrahlungSpectrum(
  wavelength: number,
  ne: number,
  Te: number
): number {
  if (wavelength <= 0 || ne <= 0 || Te <= 0) return 0;
  
  const hc_over_lambda_kT = (Constants.h * Constants.c) / (wavelength * Constants.k_B * Te);
  if (hc_over_lambda_kT > 100) return 0;
  
  const exp_term = Math.exp(-hc_over_lambda_kT);
  const prefactor = 1e-50; // Approximate constant [W·m⁵/(K^0.5)]
  
  return prefactor * ne * ne * Math.sqrt(Te) * exp_term / (wavelength * wavelength);
}

/**
 * Compute EM mode spectral contributions
 * Each mode contributes a Lorentzian line shape centered at its frequency
 */
function emModeSpectrum(
  wavelength: number,
  modes: Array<{ omega: number; re: number; im: number }>,
  decayTime: number
): number {
  let total = 0;
  const c = Constants.c;
  
  for (const mode of modes) {
    if (mode.omega <= 0) continue;
    
    const modeWavelength = (2.0 * Math.PI * c) / mode.omega;
    const modePower = (mode.re * mode.re + mode.im * mode.im) / decayTime;
    
    // Lorentzian line shape with width ~ 1/decayTime
    const deltaLambda = (c / (mode.omega * mode.omega)) * (1.0 / decayTime);
    const detuning = (wavelength - modeWavelength) / deltaLambda;
    const lorentzian = (1.0 / (Math.PI * deltaLambda)) / (1.0 + detuning * detuning);
    
    total += modePower * lorentzian;
  }
  
  return total;
}

/**
 * Compute spectral emission distribution
 * 
 * Returns power per unit wavelength across specified wavelength range
 */
export function computeSpectralDistribution(
  state: BubbleFullState,
  wavelengthMin: number = 100e-9, // 100 nm (UV)
  wavelengthMax: number = 1000e-9, // 1000 nm (near-IR)
  nBins: number = 100
): SpectralDistribution {
  const { gas, em, plasma } = state;
  const { R } = state.hydro;
  
  const wavelengths: number[] = [];
  const blackbody: number[] = [];
  const bremsstrahlung: number[] = [];
  const emModes: number[] = [];
  const powerDensity: number[] = [];
  
  const dLambda = (wavelengthMax - wavelengthMin) / nBins;
  const surfaceArea = 4.0 * Math.PI * R * R;
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  const decayTime = 1e-9; // 1 ns
  
  for (let i = 0; i < nBins; i++) {
    const lambda = wavelengthMin + (i + 0.5) * dLambda;
    wavelengths.push(lambda);
    
    // Blackbody: integrate over solid angle (4π) and convert to power per wavelength
    const B_lambda = planckSpectrum(lambda, gas.T);
    const blackbody_contrib = B_lambda * surfaceArea * 4.0 * Math.PI; // [W/m]
    blackbody.push(blackbody_contrib);
    
    // Bremsstrahlung: volume emission
    const brems_contrib = bremsstrahlungSpectrum(lambda, plasma.ne, plasma.Te) * volume; // [W/m]
    bremsstrahlung.push(brems_contrib);
    
    // EM modes: discrete contributions
    const em_contrib = emModeSpectrum(lambda, em.modes, decayTime); // [W/m]
    emModes.push(em_contrib);
    
    // Total
    powerDensity.push(blackbody_contrib + brems_contrib + em_contrib);
  }
  
  return {
    wavelengths,
    powerDensity,
    blackbody,
    bremsstrahlung,
    emModes,
  };
}

/**
 * Estimate total emission power from the bubble
 * 
 * EXPLICITLY COMBINES THREE COMPONENTS:
 * 
 * 1. BLACKBODY-ISH T TERM:
 *    P_bb = σ * T⁴ * A
 *    Where σ is Stefan-Boltzmann constant, T is gas temperature, A is bubble surface area.
 *    This represents thermal radiation from the hot gas.
 * 
 * 2. BREMSSTRAHLUNG FROM ne, Te:
 *    P_brems ~ ne² * sqrt(Te) * V
 *    Where ne is electron density, Te is electron temperature, V is bubble volume.
 *    This represents free-free radiation from electron-ion collisions in the plasma.
 * 
 * 3. EM DECAY FROM E_em:
 *    P_em = E_em / τ
 *    Where E_em is stored EM energy in the negative-space squeezed state, τ is decay time (~1 ns).
 *    This represents photons emitted from the decay of the squeezed vacuum state created
 *    by parametric amplification during collapse. This is the "light from negative cavity state."
 * 
 * Optionally computes spectral distribution if computeSpectrum is true
 */
export function estimateEmission(
  state: BubbleFullState,
  computeSpectrum: boolean = false
): EmissionSnapshot {
  const { gas, em, plasma } = state;
  const { R } = state.hydro;

  // Bubble surface area
  const surfaceArea = 4.0 * Math.PI * R * R;

  // 1. BLACKBODY-ISH T TERM: P = σ * T⁴ * A
  // Thermal radiation from hot gas (blackbody approximation)
  const blackbodyPowerDensity = Calculations.blackbodyPowerDensity(gas.T);
  const blackbodyPower = blackbodyPowerDensity * surfaceArea;

  // 2. BREMSSTRAHLUNG FROM ne, Te: P ~ ne² * sqrt(Te) * V
  // Free-free radiation from electron-ion collisions in plasma
  const volume = (4.0 / 3.0) * Math.PI * R * R * R;
  const bremsstrahlungPowerDensity = Calculations.bremsstrahlungPowerDensity(
    plasma.ne,
    plasma.Te
  );
  const bremsstrahlungPower = bremsstrahlungPowerDensity * volume;

  // 3. EM DECAY FROM E_em: P = E_em / τ
  // Photons from decay of negative-space squeezed state
  // The stored energy E_em was pumped into the cavity during extreme compression
  // via parametric amplification, creating a squeezed vacuum state. This energy
  // decays over timescale τ (~1 ns) as photons are emitted.
  const decayTime = 1e-9; // 1 ns (typical decay time for EM cavity)
  const emDecayPower = em.storedEnergy / decayTime;

  // TOTAL POWER: Sum of all three components
  const totalPower = blackbodyPower + bremsstrahlungPower + emDecayPower;

  // Compute spectrum if requested
  const spectrum = computeSpectrum
    ? computeSpectralDistribution(state, 100e-9, 1000e-9, 100)
    : undefined;

  return {
    totalPower,
    blackbodyPower,
    bremsstrahlungPower,
    emDecayPower,
    spectrum,
  };
}

/**
 * Compute gradient metrics between two states
 * 
 * Calculates time derivatives and detects extreme gradients
 * that indicate non-classical behavior.
 */
export function computeGradients(
  statePrev: BubbleFullState,
  stateNext: BubbleFullState,
  dt: number
): GradientMetrics {
  if (dt <= 0) {
    throw new Error("Time step dt must be positive");
  }

  // Compute gradients as finite differences
  const dRdt = Math.abs((stateNext.hydro.R - statePrev.hydro.R) / dt);
  const dRdotdt = Math.abs(
    (stateNext.hydro.Rdot - statePrev.hydro.Rdot) / dt
  );
  const dPgdt = Math.abs((stateNext.gas.Pg - statePrev.gas.Pg) / dt);
  const dTdt = Math.abs((stateNext.gas.T - statePrev.gas.T) / dt);
  const dnedt = Math.abs(
    (stateNext.plasma.ne - statePrev.plasma.ne) / dt
  );

  // Find maximum gradient
  const maxGradient = Math.max(dRdt, dRdotdt, dPgdt, dTdt, dnedt);

  // Extreme gradient thresholds (heuristic values)
  // These indicate when the system enters non-classical regime
  const threshold_dRdt = 1000; // m/s (supersonic)
  const threshold_dRdotdt = 1e12; // m/s² (extreme acceleration)
  const threshold_dPgdt = 1e15; // Pa/s (extreme pressure change)
  const threshold_dTdt = 1e12; // K/s (extreme heating)
  const threshold_dnedt = 1e35; // 1/(m³·s) (rapid ionization)

  const extremeGradient =
    dRdt > threshold_dRdt ||
    dRdotdt > threshold_dRdotdt ||
    dPgdt > threshold_dPgdt ||
    dTdt > threshold_dTdt ||
    dnedt > threshold_dnedt;

  return {
    dRdt,
    dRdotdt,
    dPgdt,
    dTdt,
    dnedt,
    maxGradient,
    extremeGradient,
  };
}

/**
 * Compute spatial gradient estimates (simplified)
 * 
 * Estimates spatial gradients like |∇P|, |∇n| using
 * bubble geometry and state variables.
 */
export function computeSpatialGradients(
  state: BubbleFullState
): {
  pressureGradient: number; // |∇P| [Pa/m]
  densityGradient: number; // |∇n| [1/m⁴]
  temperatureGradient: number; // |∇T| [K/m]
} {
  const { R } = state.hydro;
  const { Pg, T } = state.gas;
  const { ne } = state.plasma;

  // Simplified: assume gradients scale as ~variable / R
  // (characteristic length scale is bubble radius)
  const R_safe = Math.max(R, 1e-10);

  const pressureGradient = Math.abs(Pg) / R_safe;
  const densityGradient = Math.abs(ne) / R_safe;
  const temperatureGradient = Math.abs(T) / R_safe;

  return {
    pressureGradient,
    densityGradient,
    temperatureGradient,
  };
}
