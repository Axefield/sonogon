// config/presets.ts
// Pre-configured parameter sets for common experimental conditions

import { SonoluminescenceParams } from "../model/sonoluminescenceModel";
import { HydroParams } from "../physics/hydro";
import { ThermoParams } from "../physics/thermoChem";
import { PlasmaParams } from "../physics/plasma";
import { EmCavityParams } from "../physics/emCavity";
import { AcousticParams } from "../physics/acoustic";
import { ReactionParams } from "../physics/reactions";
import { Properties, Constants } from "../core/units";

/**
 * Default water properties at 20°C
 */
const water20C: HydroParams = {
  rho: Properties.water.density,
  mu: Properties.water.viscosity,
  sigma: Properties.water.surfaceTension,
  Pv: Properties.water.vaporPressure,
  P0: 101325, // Will be overridden in presets
  c: Properties.water.speedOfSound, // Speed of sound for Keller-Miksis
};

/**
 * Typical acoustic driving parameters
 */
const typicalAcoustic: AcousticParams = {
  Pa: 1.3e5, // 1.3 atm acoustic amplitude
  omega: 2 * Math.PI * 20e3, // 20 kHz driving frequency
};

/**
 * Default EM cavity parameters
 */
const defaultEmCavity: EmCavityParams = {
  modeFrequencies0: [
    2 * Math.PI * 3e14, // ~500 nm optical mode
    2 * Math.PI * 4e14, // ~400 nm optical mode
    2 * Math.PI * 5e14, // ~300 nm optical mode
  ],
  pumpCoefficient: 1e-15, // J/(m·s) - parametric pumping strength
  decayTime: 1e-9, // 1 ns decay time
  couplingStrength: 1e6, // coupling between modes and boundary
  refractiveIndexBase: 1.33, // water-like refractive index
};

/**
 * Default reaction parameters (Arrhenius form)
 */
const defaultReactions: ReactionParams = {
  // Reaction 0: H2O ↔ H + OH
  reaction0_A: 1e13, // pre-exponential [1/s]
  reaction0_Ea: 5e5 * Constants.R_gas, // activation energy [J/mol] (~500 kJ/mol)

  // Reaction 1: O2 ↔ 2O
  reaction1_A: 1e14,
  reaction1_Ea: 5e5 * Constants.R_gas, // ~500 kJ/mol

  // Reaction 2: H + O ↔ OH
  reaction2_A: 1e11,
  reaction2_Ea: 1e5 * Constants.R_gas, // ~100 kJ/mol

  // Reaction 3: OH ↔ O + H
  reaction3_A: 1e12,
  reaction3_Ea: 4.5e5 * Constants.R_gas, // ~450 kJ/mol

  // Reaction 4: N2 ↔ 2N
  reaction4_A: 1e15,
  reaction4_Ea: 9.5e5 * Constants.R_gas, // ~950 kJ/mol (very high)

  // Reaction 5: H + OH → H2O (three-body)
  reaction5_A: 1e10,
  reaction5_Ea: 0.5e5 * Constants.R_gas, // ~50 kJ/mol
};

/**
 * Default plasma parameters
 */
const defaultPlasma: PlasmaParams = {
  ionizationPotential_Ar: 15.76 * Constants.e, // 15.76 eV
  ionizationPotential_Xe: 12.13 * Constants.e, // 12.13 eV
  ionizationPotential_H: 13.6 * Constants.e, // 13.6 eV
  ionizationPotential_O: 13.62 * Constants.e, // 13.62 eV
  recombinationCoeff: 1e-13, // m³/s
  electronCollisionFreq: 1e13, // 1/s
};

/**
 * Preset: Water with dissolved air (typical experimental setup)
 */
export function createWaterAirPreset(): SonoluminescenceParams {
  return {
    hydro: {
      ...water20C,
      P0: 101325, // 1 atm ambient pressure
    },
    thermo: {
      gamma: 1.4, // air-like (diatomic)
      R0: 1e-6, // 1 micron reference radius
      Pg0: 101325, // 1 atm reference pressure
      T0: 293.15, // 20°C reference temperature
      heatLossCoeff: 100, // W/(m²·K) - moderate heat loss
    },
    plasma: defaultPlasma,
    em: defaultEmCavity,
    acoustic: typicalAcoustic,
    reactions: defaultReactions,
  };
}

/**
 * Preset: Argon bubble in water (common for sonoluminescence)
 */
export function createArgonBubblePreset(): SonoluminescenceParams {
  return {
    hydro: {
      ...water20C,
      P0: 101325, // 1 atm
    },
    thermo: {
      gamma: Properties.argon.gamma, // 1.67 for monatomic
      R0: 1e-6, // 1 micron
      Pg0: 101325,
      T0: 293.15,
      heatLossCoeff: 50, // Lower heat loss (argon is inert)
    },
    plasma: {
      ...defaultPlasma,
      // Argon is the dominant species, so use its ionization potential primarily
    },
    em: defaultEmCavity,
    acoustic: typicalAcoustic,
    reactions: defaultReactions,
  };
}

/**
 * Preset: Xenon bubble (high light yield)
 */
export function createXenonBubblePreset(): SonoluminescenceParams {
  return {
    hydro: {
      ...water20C,
      P0: 101325,
    },
    thermo: {
      gamma: Properties.xenon.gamma, // 1.67 for monatomic
      R0: 1e-6,
      Pg0: 101325,
      T0: 293.15,
      heatLossCoeff: 50,
    },
    plasma: {
      ...defaultPlasma,
      // Xenon has lower ionization potential, easier to ionize
    },
    em: {
      ...defaultEmCavity,
      pumpCoefficient: 2e-15, // Higher pumping for xenon
    },
    acoustic: typicalAcoustic,
    reactions: defaultReactions,
  };
}

/**
 * Preset: High-intensity driving (extreme conditions)
 */
export function createHighIntensityPreset(): SonoluminescenceParams {
  return {
    hydro: {
      ...water20C,
      P0: 101325,
    },
    thermo: {
      gamma: 1.4,
      R0: 1e-6,
      Pg0: 101325,
      T0: 293.15,
      heatLossCoeff: 200, // Higher heat loss at high intensity
    },
    plasma: defaultPlasma,
    em: {
      ...defaultEmCavity,
      pumpCoefficient: 5e-15, // Stronger pumping
      decayTime: 5e-10, // Faster decay
    },
    acoustic: {
      Pa: 2.0e5, // Higher acoustic amplitude (2 atm)
      omega: 2 * Math.PI * 25e3, // 25 kHz
    },
    reactions: defaultReactions,
  };
}

/**
 * Preset: Pistol Shrimp (Alpheidae) 
 * 
 * The pistol shrimp creates cavitation bubbles through rapid claw closure,
 * ejecting a high-speed water jet (30-60 m/s) that forms a cavitation bubble.
 * The bubble collapses with extreme violence, generating:
 * - Temperatures up to 8,000°F (4,427°C) - comparable to sun's surface
 * - Sound intensity up to 218 dB (louder than gunshot/jet engine)
 * - Sonoluminescence (light flash)
 * - Shockwaves that stun or kill prey
 * 
 * Typical conditions:
 * - Initial bubble radius: 1-5 mm
 * - Jet velocity: 30-60 m/s
 * - Jet duration: ~0.1-1 ms
 * - Collapse velocity: can exceed sound speed (supersonic)
 */
export function createPistolShrimpPreset(): SonoluminescenceParams {
  return {
    hydro: {
      ...water20C,
      P0: 101325, // 1 atm (underwater)
      useKellerMiksis: true, // Use Keller-Miksis for supersonic collapse
    },
    thermo: {
      gamma: 1.4, // Air-like (bubble contains air/water vapor)
      R0: 2e-3, // 2 mm initial bubble (typical for pistol shrimp)
      Pg0: 101325, // 1 atm
      T0: 293.15, // 20°C (ambient water temperature)
      heatLossCoeff: 1000, // Very high heat loss (water is good conductor)
    },
    plasma: {
      ...defaultPlasma,
      // Pistol shrimp conditions can reach extreme temperatures
      // Ionization will occur during collapse
    },
    em: {
      ...defaultEmCavity,
      pumpCoefficient: 1e-14, // Strong parametric pumping from rapid collapse
      decayTime: 1e-9, // 1 ns decay (fast emission)
    },
    acoustic: {
      // Jet-driven cavitation (pistol shrimp mechanism)
      jetDriven: true,
      jetVelocity: 45, // 45 m/s (typical pistol shrimp jet velocity)
      jetDuration: 5e-4, // 0.5 ms jet duration
      jetNozzleArea: 1e-6, // ~1 mm² effective nozzle area
      jetDirection: { x: 1, y: 0, z: 0 }, // Jet direction (normalized)
      // Also include some acoustic driving (optional, can be zero)
      Pa: 0, // No additional acoustic driving (pure jet-driven)
      omega: 0,
    },
    reactions: defaultReactions,
  };
}

/**
 * Validate parameter set for physical reasonableness
 */
export function validateParams(params: SonoluminescenceParams): {
  valid: boolean;
  errors: string[];
} {
  const errors: string[] = [];

  // Hydro validation
  if (params.hydro.rho <= 0) {
    errors.push("Fluid density must be positive");
  }
  if (params.hydro.mu < 0) {
    errors.push("Viscosity cannot be negative");
  }
  if (params.hydro.sigma < 0) {
    errors.push("Surface tension cannot be negative");
  }

  // Thermo validation
  if (params.thermo.gamma < 1.0 || params.thermo.gamma > 2.0) {
    errors.push("Adiabatic index gamma should be between 1 and 2");
  }
  if (params.thermo.T0 <= 0) {
    errors.push("Reference temperature must be positive");
  }

  // Acoustic validation
  if (params.acoustic.Pa !== undefined && params.acoustic.Pa < 0) {
    errors.push("Acoustic amplitude cannot be negative");
  }
  if (params.acoustic.omega !== undefined && params.acoustic.omega <= 0) {
    errors.push("Acoustic frequency must be positive");
  }
  if (!params.acoustic.Pa && !params.acoustic.frequencies) {
    errors.push("Either Pa/omega or frequencies must be specified");
  }

  // EM validation
  if (params.em.decayTime <= 0) {
    errors.push("EM decay time must be positive");
  }
  if (params.em.modeFrequencies0.length === 0) {
    errors.push("At least one EM mode frequency must be specified");
  }

  // Plasma validation
  if (params.plasma.recombinationCoeff < 0) {
    errors.push("Recombination coefficient cannot be negative");
  }

  return {
    valid: errors.length === 0,
    errors,
  };
}

/**
 * Get default preset (water with air)
 */
export function getDefaultPreset(): SonoluminescenceParams {
  return createWaterAirPreset();
}

