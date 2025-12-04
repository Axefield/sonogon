"use strict";
// core/units.ts
// Physical constants and unit system for sonoluminescence model
Object.defineProperty(exports, "__esModule", { value: true });
exports.Calculations = exports.Units = exports.ReferenceValues = exports.Properties = exports.Constants = void 0;
/**
 * Fundamental physical constants (SI units)
 */
exports.Constants = {
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
};
/**
 * Common physical properties (SI units)
 */
exports.Properties = {
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
};
/**
 * Reference values for normalization
 */
exports.ReferenceValues = {
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
};
/**
 * Unit conversion utilities
 */
exports.Units = {
    /**
     * Convert pressure from Pa to normalized units
     */
    pressureToNormalized(P) {
        return P / exports.ReferenceValues.P0_ref;
    },
    /**
     * Convert pressure from normalized units to Pa
     */
    pressureFromNormalized(P_norm) {
        return P_norm * exports.ReferenceValues.P0_ref;
    },
    /**
     * Convert temperature from K to normalized units
     */
    temperatureToNormalized(T) {
        return T / exports.ReferenceValues.T0_ref;
    },
    /**
     * Convert temperature from normalized units to K
     */
    temperatureFromNormalized(T_norm) {
        return T_norm * exports.ReferenceValues.T0_ref;
    },
    /**
     * Convert length from m to normalized units
     */
    lengthToNormalized(L) {
        return L / exports.ReferenceValues.R0_ref;
    },
    /**
     * Convert length from normalized units to m
     */
    lengthFromNormalized(L_norm) {
        return L_norm * exports.ReferenceValues.R0_ref;
    },
    /**
     * Convert time from s to normalized units
     */
    timeToNormalized(t) {
        return t / exports.ReferenceValues.t0_ref;
    },
    /**
     * Convert time from normalized units to s
     */
    timeFromNormalized(t_norm) {
        return t_norm * exports.ReferenceValues.t0_ref;
    },
    /**
     * Convert number density from 1/m³ to normalized units
     */
    numberDensityToNormalized(n) {
        const n0_ref = exports.ReferenceValues.P0_ref / (exports.Constants.k_B * exports.ReferenceValues.T0_ref);
        return n / n0_ref;
    },
    /**
     * Convert number density from normalized units to 1/m³
     */
    numberDensityFromNormalized(n_norm) {
        const n0_ref = exports.ReferenceValues.P0_ref / (exports.Constants.k_B * exports.ReferenceValues.T0_ref);
        return n_norm * n0_ref;
    },
};
/**
 * Helper functions for common calculations
 */
exports.Calculations = {
    /**
     * Ideal gas law: P = n * k_B * T
     */
    idealGasPressure(n, T) {
        return n * exports.Constants.k_B * T;
    },
    /**
     * Number density from pressure and temperature: n = P / (k_B * T)
     */
    numberDensityFromPressure(P, T) {
        return P / (exports.Constants.k_B * T);
    },
    /**
     * Plasma frequency: ωp = sqrt(ne * e² / (ε0 * me))
     */
    plasmaFrequency(ne) {
        return Math.sqrt((ne * exports.Constants.e * exports.Constants.e) / (exports.Constants.epsilon_0 * exports.Constants.m_e));
    },
    /**
     * Debye length: λD = sqrt(ε0 * k_B * T / (ne * e²))
     */
    debyeLength(ne, T) {
        return Math.sqrt((exports.Constants.epsilon_0 * exports.Constants.k_B * T) / (ne * exports.Constants.e * exports.Constants.e));
    },
    /**
     * Blackbody power density: P = σ * T⁴
     */
    blackbodyPowerDensity(T) {
        return exports.Constants.sigma_SB * T * T * T * T;
    },
    /**
     * Bremsstrahlung power density (simplified): P ~ ne² * sqrt(Te)
     */
    bremsstrahlungPowerDensity(ne, Te) {
        // Simplified bremsstrahlung formula
        // More accurate: P ~ ne² * sqrt(Te) * Z² * g_ff
        const prefactor = 1e-40; // W·m³/(K^0.5) - approximate
        return prefactor * ne * ne * Math.sqrt(Te);
    },
};
//# sourceMappingURL=units.js.map