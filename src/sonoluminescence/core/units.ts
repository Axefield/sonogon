// core/units.ts
// Physical constants and unit system for sonoluminescence model

/**
 * Fundamental physical constants (SI units)
 */
export const Constants = {
  // Gas constant
  R_gas: 8.314462618, // J/(mol·K)
  
  // Boltzmann constant
  k_B: 1.380649e-23, // J/K
  
  // Avogadro's number
  N_A: 6.02214076e23, // 1/mol
  
  // Electron charge
  e: 1.602176634e-19, // C
  
  // Electron mass
  m_e: 9.1093837015e-31, // kg
  
  // Proton mass
  m_p: 1.67262192369e-27, // kg
  
  // Vacuum permittivity
  epsilon_0: 8.8541878128e-12, // F/m
  
  // Vacuum permeability
  mu_0: 1.25663706212e-6, // H/m
  
  // Speed of light
  c: 2.99792458e8, // m/s
  
  // Planck constant
  h: 6.62607015e-34, // J·s
  
  // Reduced Planck constant
  hbar: 1.054571817e-34, // J·s
  
  // Stefan-Boltzmann constant
  sigma_SB: 5.670374419e-8, // W/(m²·K⁴)
} as const;

/**
 * Common physical properties (SI units)
 */
export const Properties = {
  // Water properties at 20°C
  water: {
    density: 998.2, // kg/m³
    viscosity: 1.002e-3, // Pa·s
    surfaceTension: 0.0728, // N/m
    vaporPressure: 2339, // Pa (at 20°C)
    speedOfSound: 1482, // m/s
  },
  
  // Air properties at 20°C, 1 atm
  air: {
    density: 1.204, // kg/m³
    viscosity: 1.825e-5, // Pa·s
  },
  
  // Argon properties
  argon: {
    molecularWeight: 39.948e-3, // kg/mol
    gamma: 1.67, // adiabatic index
  },
  
  // Xenon properties
  xenon: {
    molecularWeight: 131.293e-3, // kg/mol
    gamma: 1.67, // adiabatic index
  },
} as const;

/**
 * Reference values for normalization
 */
export const ReferenceValues = {
  // Length
  R0_ref: 1e-6, // 1 micron (typical bubble size)
  
  // Time
  t0_ref: 1e-6, // 1 microsecond (acoustic period scale)
  
  // Pressure
  P0_ref: 101325, // 1 atm (standard atmospheric pressure)
  
  // Temperature
  T0_ref: 293.15, // 20°C (room temperature)
  
  // Density
  rho0_ref: 998.2, // water density at 20°C
  
  // Energy
  E0_ref: 1e-18, // 1 attojoule (typical bubble energy scale)
} as const;

/**
 * Unit conversion utilities
 */
export const Units = {
  /**
   * Convert pressure from Pa to normalized units
   */
  pressureToNormalized(P: number): number {
    return P / ReferenceValues.P0_ref;
  },
  
  /**
   * Convert pressure from normalized units to Pa
   */
  pressureFromNormalized(P_norm: number): number {
    return P_norm * ReferenceValues.P0_ref;
  },
  
  /**
   * Convert temperature from K to normalized units
   */
  temperatureToNormalized(T: number): number {
    return T / ReferenceValues.T0_ref;
  },
  
  /**
   * Convert temperature from normalized units to K
   */
  temperatureFromNormalized(T_norm: number): number {
    return T_norm * ReferenceValues.T0_ref;
  },
  
  /**
   * Convert length from m to normalized units
   */
  lengthToNormalized(L: number): number {
    return L / ReferenceValues.R0_ref;
  },
  
  /**
   * Convert length from normalized units to m
   */
  lengthFromNormalized(L_norm: number): number {
    return L_norm * ReferenceValues.R0_ref;
  },
  
  /**
   * Convert time from s to normalized units
   */
  timeToNormalized(t: number): number {
    return t / ReferenceValues.t0_ref;
  },
  
  /**
   * Convert time from normalized units to s
   */
  timeFromNormalized(t_norm: number): number {
    return t_norm * ReferenceValues.t0_ref;
  },
  
  /**
   * Convert number density from 1/m³ to normalized units
   */
  numberDensityToNormalized(n: number): number {
    const n0_ref = ReferenceValues.P0_ref / (Constants.k_B * ReferenceValues.T0_ref);
    return n / n0_ref;
  },
  
  /**
   * Convert number density from normalized units to 1/m³
   */
  numberDensityFromNormalized(n_norm: number): number {
    const n0_ref = ReferenceValues.P0_ref / (Constants.k_B * ReferenceValues.T0_ref);
    return n_norm * n0_ref;
  },
};

/**
 * Helper functions for common calculations
 */
export const Calculations = {
  /**
   * Ideal gas law: P = n * k_B * T
   */
  idealGasPressure(n: number, T: number): number {
    return n * Constants.k_B * T;
  },
  
  /**
   * Number density from pressure and temperature: n = P / (k_B * T)
   */
  numberDensityFromPressure(P: number, T: number): number {
    return P / (Constants.k_B * T);
  },
  
  /**
   * Plasma frequency: ωp = sqrt(ne * e² / (ε0 * me))
   */
  plasmaFrequency(ne: number): number {
    return Math.sqrt((ne * Constants.e * Constants.e) / (Constants.epsilon_0 * Constants.m_e));
  },
  
  /**
   * Debye length: λD = sqrt(ε0 * k_B * T / (ne * e²))
   */
  debyeLength(ne: number, T: number): number {
    return Math.sqrt((Constants.epsilon_0 * Constants.k_B * T) / (ne * Constants.e * Constants.e));
  },
  
  /**
   * Blackbody power density: P = σ * T⁴
   */
  blackbodyPowerDensity(T: number): number {
    return Constants.sigma_SB * T * T * T * T;
  },
  
  /**
   * Bremsstrahlung power density (simplified): P ~ ne² * sqrt(Te)
   */
  bremsstrahlungPowerDensity(ne: number, Te: number): number {
    // Simplified bremsstrahlung formula
    // More accurate: P ~ ne² * sqrt(Te) * Z² * g_ff
    const prefactor = 1e-40; // W·m³/(K^0.5) - approximate
    return prefactor * ne * ne * Math.sqrt(Te);
  },
};

