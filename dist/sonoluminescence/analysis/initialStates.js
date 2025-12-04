"use strict";
// analysis/initialStates.ts
// Helper functions to generate initial states for common scenarios
Object.defineProperty(exports, "__esModule", { value: true });
exports.createEquilibriumState = createEquilibriumState;
exports.createExpandedState = createExpandedState;
exports.createCollapseState = createCollapseState;
exports.createFromPreset = createFromPreset;
exports.validateInitialState = validateInitialState;
const units_1 = require("../core/units");
/**
 * Create an equilibrium initial state
 *
 * Bubble at acoustic equilibrium: R = R_eq, Rdot = 0
 * Gas pressure balances acoustic pressure
 */
function createEquilibriumState(params, R_equilibrium) {
    const { acoustic, thermo, hydro } = params;
    // Equilibrium: Pg = P0 + Pacoustic (at phase = 0, Pacoustic = 0)
    // For simplicity, use P0 as equilibrium pressure
    const Pg_eq = hydro.P0;
    // Temperature at equilibrium (ambient)
    const T_eq = thermo.T0;
    // Compute number density from ideal gas law
    const n_total = units_1.Calculations.numberDensityFromPressure(Pg_eq, T_eq);
    // Distribute species (simplified: assume mostly Ar for argon preset)
    const speciesDensities = {
        H2O: n_total * 0.01, // 1% water vapor
        O2: n_total * 0.02, // 2% oxygen
        N2: n_total * 0.02, // 2% nitrogen
        Ar: n_total * 0.90, // 90% argon (dominant)
        Xe: n_total * 0.0, // No xenon
        H: n_total * 0.0, // No atomic hydrogen
        O: n_total * 0.0, // No atomic oxygen
        OH: n_total * 0.0, // No hydroxyl
        N: n_total * 0.0, // No atomic nitrogen
    };
    return {
        t: 0,
        hydro: {
            R: R_equilibrium,
            Rdot: 0, // At equilibrium, no motion
        },
        gas: {
            Pg: Pg_eq,
            T: T_eq,
        },
        species: {
            numberDensity: speciesDensities,
        },
        plasma: {
            ne: 1e18, // Low electron density at equilibrium
            Te: T_eq, // Electron temperature equals gas temperature
            ionizationFraction: 1e-6, // Very low ionization
        },
        internalEnergy: {
            // Equilibrium energy partitions
            translational: (3.0 / 2.0) * units_1.Constants.k_B * T_eq * n_total * (4.0 / 3.0) * Math.PI * R_equilibrium * R_equilibrium * R_equilibrium,
            rotational: 2.0 * units_1.Constants.k_B * T_eq * n_total * (4.0 / 3.0) * Math.PI * R_equilibrium * R_equilibrium * R_equilibrium,
            vibrational: units_1.Constants.k_B * T_eq * n_total * (4.0 / 3.0) * Math.PI * R_equilibrium * R_equilibrium * R_equilibrium,
            electronic: 0.01 * units_1.Constants.k_B * T_eq * n_total * (4.0 / 3.0) * Math.PI * R_equilibrium * R_equilibrium * R_equilibrium,
        },
        em: {
            modes: [
                {
                    omega: params.em.modeFrequencies0[0] || 2 * Math.PI * 3e14,
                    re: 0,
                    im: 0,
                },
                {
                    omega: params.em.modeFrequencies0[1] || 2 * Math.PI * 4e14,
                    re: 0,
                    im: 0,
                },
                {
                    omega: params.em.modeFrequencies0[2] || 2 * Math.PI * 5e14,
                    re: 0,
                    im: 0,
                },
            ],
            storedEnergy: 1e-25, // Very small initial stored energy
        },
        acoustic: {
            phase: 0, // Start at phase 0
        },
        reactions: {
            xi: [0, 0, 0], // No reaction progress initially
        },
    };
}
/**
 * Create initial state at maximum expansion
 *
 * Bubble at maximum radius (typically 10-50x equilibrium)
 * Rdot = 0 (turning point)
 */
function createExpandedState(params, R_max) {
    const equilibrium = createEquilibriumState(params, params.thermo.R0);
    // At maximum expansion:
    // - Radius is maximum
    // - Velocity is zero (turning point)
    // - Pressure is lower (expanded)
    // - Temperature is lower (adiabatic expansion)
    const R_eq = params.thermo.R0;
    const gamma = params.thermo.gamma;
    // Adiabatic expansion: P * V^gamma = constant
    const volume_ratio = Math.pow(R_max / R_eq, 3);
    const Pg_max = equilibrium.gas.Pg / Math.pow(volume_ratio, gamma);
    // Adiabatic: T * V^(gamma-1) = constant
    const T_max = equilibrium.gas.T / Math.pow(volume_ratio, gamma - 1);
    // Number densities scale with volume
    const n_scale = 1.0 / volume_ratio;
    const speciesDensities = {
        H2O: equilibrium.species.numberDensity.H2O * n_scale,
        O2: equilibrium.species.numberDensity.O2 * n_scale,
        N2: equilibrium.species.numberDensity.N2 * n_scale,
        Ar: equilibrium.species.numberDensity.Ar * n_scale,
        Xe: equilibrium.species.numberDensity.Xe * n_scale,
        H: equilibrium.species.numberDensity.H * n_scale,
        O: equilibrium.species.numberDensity.O * n_scale,
        OH: (equilibrium.species.numberDensity.OH || 0) * n_scale,
        N: (equilibrium.species.numberDensity.N || 0) * n_scale,
    };
    return {
        ...equilibrium,
        hydro: {
            R: R_max,
            Rdot: 0, // Turning point
        },
        gas: {
            Pg: Pg_max,
            T: T_max,
        },
        species: {
            numberDensity: speciesDensities,
        },
        plasma: {
            ne: equilibrium.plasma.ne * n_scale,
            Te: T_max,
            ionizationFraction: equilibrium.plasma.ionizationFraction,
        },
    };
}
/**
 * Create initial state just before collapse
 *
 * Bubble starting to collapse: R = R_min_start, Rdot < 0
 */
function createCollapseState(params, R_min_start, Rdot_collapse = -100 // m/s, negative for collapse
) {
    const equilibrium = createEquilibriumState(params, params.thermo.R0);
    // Before collapse, bubble is compressed
    const R_eq = params.thermo.R0;
    const gamma = params.thermo.gamma;
    // Adiabatic compression
    const volume_ratio = Math.pow(R_min_start / R_eq, 3);
    const Pg_compressed = equilibrium.gas.Pg * Math.pow(volume_ratio, gamma);
    const T_compressed = equilibrium.gas.T * Math.pow(volume_ratio, gamma - 1);
    // Number densities increase with compression
    const n_scale = 1.0 / volume_ratio;
    const speciesDensities = {
        H2O: equilibrium.species.numberDensity.H2O * n_scale,
        O2: equilibrium.species.numberDensity.O2 * n_scale,
        N2: equilibrium.species.numberDensity.N2 * n_scale,
        Ar: equilibrium.species.numberDensity.Ar * n_scale,
        Xe: equilibrium.species.numberDensity.Xe * n_scale,
        H: equilibrium.species.numberDensity.H * n_scale,
        O: equilibrium.species.numberDensity.O * n_scale,
        OH: (equilibrium.species.numberDensity.OH || 0) * n_scale,
        N: (equilibrium.species.numberDensity.N || 0) * n_scale,
    };
    // Higher temperature may start ionization
    const n_total = Object.values(speciesDensities).reduce((sum, n) => sum + n, 0);
    const ionizationPotential = 15.76 * units_1.Constants.e; // Argon
    const ne_estimate = units_1.Calculations.numberDensityFromPressure(Pg_compressed, T_compressed) * 0.01; // Rough estimate
    return {
        ...equilibrium,
        hydro: {
            R: R_min_start,
            Rdot: Rdot_collapse,
        },
        gas: {
            Pg: Pg_compressed,
            T: T_compressed,
        },
        species: {
            numberDensity: speciesDensities,
        },
        plasma: {
            ne: ne_estimate,
            Te: T_compressed,
            ionizationFraction: ne_estimate / n_total,
        },
    };
}
/**
 * Create initial state from preset with custom radius
 */
function createFromPreset(params, R_initial, Rdot_initial = 0) {
    const state = createEquilibriumState(params, R_initial);
    return {
        ...state,
        hydro: {
            R: R_initial,
            Rdot: Rdot_initial,
        },
    };
}
/**
 * Validate initial state for physical consistency
 */
function validateInitialState(state) {
    const errors = [];
    // Check radius
    if (state.hydro.R <= 0) {
        errors.push("Bubble radius must be positive");
    }
    if (state.hydro.R < 1e-9) {
        errors.push("Bubble radius too small (< 1 nm)");
    }
    if (state.hydro.R > 1e-3) {
        errors.push("Bubble radius too large (> 1 mm)");
    }
    // Check pressure
    if (state.gas.Pg <= 0) {
        errors.push("Gas pressure must be positive");
    }
    // Check temperature
    if (state.gas.T <= 0) {
        errors.push("Temperature must be positive");
    }
    if (state.gas.T > 1e6) {
        errors.push("Temperature unreasonably high (> 1 MK)");
    }
    // Check species densities
    const totalSpecies = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0);
    if (totalSpecies <= 0) {
        errors.push("Total species number density must be positive");
    }
    // Check plasma
    if (state.plasma.ne < 0) {
        errors.push("Electron density cannot be negative");
    }
    if (state.plasma.Te <= 0) {
        errors.push("Electron temperature must be positive");
    }
    if (state.plasma.ionizationFraction < 0 || state.plasma.ionizationFraction > 1) {
        errors.push("Ionization fraction must be between 0 and 1");
    }
    // Check EM modes
    if (state.em.modes.length === 0) {
        errors.push("At least one EM mode must be specified");
    }
    return {
        valid: errors.length === 0,
        errors,
    };
}
//# sourceMappingURL=initialStates.js.map