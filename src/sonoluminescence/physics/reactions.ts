// physics/reactions.ts
import { BubbleFullState } from "../model/types";
import { Constants } from "../core/units";

export interface ReactionParams {
  // Reaction rate constants (Arrhenius form: k = A * exp(-Ea / (R*T)))
  // Reaction 0: H2O ↔ H + OH (water dissociation)
  reaction0_A: number;  // pre-exponential factor [1/s]
  reaction0_Ea: number; // activation energy [J/mol]
  
  // Reaction 1: O2 ↔ 2O (oxygen dissociation)
  reaction1_A: number;
  reaction1_Ea: number;
  
  // Reaction 2: H + O ↔ OH (recombination)
  reaction2_A: number;
  reaction2_Ea: number;
  
  // Reaction 3: OH ↔ O + H (hydroxyl dissociation)
  reaction3_A?: number;
  reaction3_Ea?: number;
  
  // Reaction 4: N2 ↔ 2N (nitrogen dissociation)
  reaction4_A?: number;
  reaction4_Ea?: number;
  
  // Reaction 5: H + OH → H2O (three-body recombination)
  reaction5_A?: number;
  reaction5_Ea?: number;
  
  // Pressure-dependent rates (Lindemann mechanism)
  usePressureDependentRates?: boolean; // Enable Lindemann falloff
  reaction0_k0?: number; // Low-pressure limit rate constant [m³/(mol·s)]
  reaction0_kInf?: number; // High-pressure limit rate constant [1/s]
  reaction0_Fc?: number; // Broadening factor (Troe parameter)
  
  // Detailed three-body reactions
  useDetailedThreeBody?: boolean; // Use species-specific third-body efficiencies
  thirdBodyEfficiencies?: Record<string, number>; // Efficiency factors for each species
}

export interface ReactionDerivatives {
  dXidt: number[]; // reaction progress derivatives
}

/**
 * Arrhenius rate constant: k = A * exp(-Ea / (R*T))
 */
function arrheniusRate(A: number, Ea: number, T: number): number {
  if (T <= 0) return 0;
  return A * Math.exp(-Ea / (Constants.R_gas * T));
}

/**
 * Lindemann falloff rate constant
 * 
 * For unimolecular/bimolecular reactions at high pressure:
 * k(P) = k_inf * (Pr / (1 + Pr)) * F
 * 
 * where Pr = k0 * [M] / k_inf is the reduced pressure
 * and F is the broadening factor (Troe parameter)
 * 
 * At low pressure: k → k0 * [M]
 * At high pressure: k → k_inf
 */
function lindemannRate(
  k0: number,      // Low-pressure limit [m³/(mol·s)]
  kInf: number,     // High-pressure limit [1/s]
  M: number,        // Third-body concentration [mol/m³]
  Fc: number = 0.6  // Broadening factor (Troe parameter)
): number {
  if (M <= 0) return kInf; // No third body, use high-pressure limit
  
  const Pr = (k0 * M) / kInf; // Reduced pressure
  const F = Math.pow(10, Math.log10(Pr) / (1.0 + Math.pow(Math.log10(Pr), 2.0))); // Broadening
  const F_troe = Fc * F; // Troe broadening
  
  return kInf * (Pr / (1.0 + Pr)) * F_troe;
}

/**
 * Compute reaction progress derivatives
 * 
 * Reaction progress variables xi represent the extent of each reaction.
 * For reaction: aA + bB ↔ cC + dD
 * 
 * Forward rate: r_f = k_f * [A]^a * [B]^b
 * Reverse rate: r_r = k_r * [C]^c * [D]^d
 * 
 * dxi/dt = r_f - r_r
 * 
 * The reactions affect species number densities:
 * - Reaction 0: H2O ↔ H + OH
 * - Reaction 1: O2 ↔ 2O
 * - Reaction 2: H + O ↔ OH (recombination)
 */
export function computeReactionDerivatives(
  state: BubbleFullState,
  params: ReactionParams
): ReactionDerivatives {
  const { T } = state.gas;
  const { numberDensity } = state.species;
  const { xi } = state.reactions;

  // Ensure we have at least 3 reaction progress variables
  const nReactions = Math.max(xi.length, 3);
  const dXidt = new Array(nReactions).fill(0);

  // Reaction 0: H2O ↔ H + OH
  // Forward: H2O → H + OH
  // Reverse: H + OH → H2O
  if (nReactions > 0) {
    let k_f0 = arrheniusRate(params.reaction0_A, params.reaction0_Ea, T);
    
    // Apply pressure-dependent rate if enabled
    if (params.usePressureDependentRates && params.reaction0_k0 && params.reaction0_kInf) {
      const n_total = Object.values(numberDensity).reduce((sum, n) => sum + n, 0);
      const volume = (4.0 / 3.0) * Math.PI * Math.pow(state.hydro.R, 3);
      const M = (n_total * volume) / Constants.N_A; // Third-body concentration [mol/m³]
      k_f0 = lindemannRate(
        params.reaction0_k0,
        params.reaction0_kInf || k_f0,
        M,
        params.reaction0_Fc || 0.6
      );
    }
    
    const k_r0 = arrheniusRate(params.reaction0_A * 0.1, params.reaction0_Ea * 0.9, T); // Reverse typically slower
    
    const r_f0 = k_f0 * numberDensity.H2O;
    const r_r0 = k_r0 * numberDensity.H * numberDensity.O; // Simplified: treating O as OH
    
    dXidt[0] = r_f0 - r_r0;
  }

  // Reaction 1: O2 ↔ 2O
  // Forward: O2 → 2O
  // Reverse: 2O → O2
  if (nReactions > 1) {
    const k_f1 = arrheniusRate(params.reaction1_A, params.reaction1_Ea, T);
    const k_r1 = arrheniusRate(params.reaction1_A * 0.01, params.reaction1_Ea * 0.8, T);
    
    const r_f1 = k_f1 * numberDensity.O2;
    const r_r1 = k_r1 * numberDensity.O * numberDensity.O;
    
    dXidt[1] = r_f1 - r_r1;
  }

  // Reaction 2: H + O ↔ OH (recombination)
  // Forward: H + O → OH
  // Reverse: OH → H + O (dissociation)
  if (nReactions > 2) {
    const k_f2 = arrheniusRate(params.reaction2_A, params.reaction2_Ea, T);
    const k_r2 = arrheniusRate(params.reaction2_A * 0.1, params.reaction2_Ea * 0.9, T);
    
    const r_f2 = k_f2 * numberDensity.H * numberDensity.O;
    const r_r2 = k_r2 * (numberDensity.OH || 0);
    
    dXidt[2] = r_f2 - r_r2;
  }

  // Reaction 3: OH ↔ O + H (hydroxyl dissociation)
  if (nReactions > 3 && params.reaction3_A && params.reaction3_Ea) {
    const k_f3 = arrheniusRate(params.reaction3_A, params.reaction3_Ea, T);
    const k_r3 = arrheniusRate(params.reaction3_A * 0.1, params.reaction3_Ea * 0.9, T);
    
    const r_f3 = k_f3 * (numberDensity.OH || 0);
    const r_r3 = k_r3 * numberDensity.O * numberDensity.H;
    
    dXidt[3] = r_f3 - r_r3;
  }

  // Reaction 4: N2 ↔ 2N (nitrogen dissociation)
  if (nReactions > 4 && params.reaction4_A && params.reaction4_Ea) {
    const k_f4 = arrheniusRate(params.reaction4_A, params.reaction4_Ea, T);
    const k_r4 = arrheniusRate(params.reaction4_A * 0.01, params.reaction4_Ea * 0.8, T);
    
    const r_f4 = k_f4 * numberDensity.N2;
    const r_r4 = k_r4 * (numberDensity.N || 0) * (numberDensity.N || 0);
    
    dXidt[4] = r_f4 - r_r4;
  }

  // Reaction 5: H + OH → H2O (three-body recombination)
  if (nReactions > 5 && params.reaction5_A && params.reaction5_Ea) {
    const k_f5 = arrheniusRate(params.reaction5_A, params.reaction5_Ea, T);
    
    // Three-body rate: k * [H] * [OH] * [M_eff]
    // where [M_eff] is the effective third-body concentration
    let M_eff = 0;
    
    if (params.useDetailedThreeBody && params.thirdBodyEfficiencies) {
      // Detailed: weighted sum with species-specific efficiencies
      // M_eff = sum_i (alpha_i * [i])
      // where alpha_i is the efficiency factor for species i
      const efficiencies = params.thirdBodyEfficiencies;
      for (const [species, n] of Object.entries(numberDensity)) {
        const alpha = efficiencies[species] || 1.0; // Default efficiency = 1.0
        M_eff += alpha * n;
      }
    } else {
      // Simplified: use total density (all species equally efficient)
      M_eff = Object.values(numberDensity).reduce((sum, n) => sum + n, 0);
    }
    
    const r_f5 = k_f5 * numberDensity.H * (numberDensity.OH || 0) * M_eff;
    // Reverse is slow (H2O dissociation)
    const k_r5 = arrheniusRate(params.reaction5_A * 0.001, params.reaction5_Ea * 1.1, T);
    const r_r5 = k_r5 * numberDensity.H2O;
    
    dXidt[5] = r_f5 - r_r5;
  }

  return { dXidt };
}

/**
 * Compute species number density changes from reactions
 * 
 * This should be called in addition to computeSpeciesDerivatives from thermoChem
 * to add reaction terms to the species evolution
 */
export function computeSpeciesReactionRates(
  state: BubbleFullState,
  params: ReactionParams
): Record<string, number> {
  const { T } = state.gas;
  const { numberDensity } = state.species;
  const { xi } = state.reactions;

  const reactionDeriv = computeReactionDerivatives(state, params);
  
  const rates: Record<string, number> = {
    H2O: 0,
    O2: 0,
    N2: 0,
    Ar: 0,
    Xe: 0,
    H: 0,
    O: 0,
    OH: 0,
    N: 0,
  };

  // Reaction 0: H2O ↔ H + OH
  if (xi.length > 0) {
    const dxi0_dt = reactionDeriv.dXidt[0];
    rates.H2O -= dxi0_dt; // H2O consumed
    rates.H += dxi0_dt;    // H produced
    rates.OH += dxi0_dt;   // OH produced
  }

  // Reaction 1: O2 ↔ 2O
  if (xi.length > 1) {
    const dxi1_dt = reactionDeriv.dXidt[1];
    rates.O2 -= dxi1_dt;   // O2 consumed
    rates.O += 2.0 * dxi1_dt; // 2O produced
  }

  // Reaction 2: H + O ↔ OH
  if (xi.length > 2) {
    const dxi2_dt = reactionDeriv.dXidt[2];
    rates.H -= dxi2_dt;    // H consumed
    rates.O -= dxi2_dt;    // O consumed
    rates.OH += dxi2_dt;   // OH produced
  }

  // Reaction 3: OH ↔ O + H
  if (xi.length > 3 && reactionDeriv.dXidt.length > 3) {
    const dxi3_dt = reactionDeriv.dXidt[3];
    rates.OH -= dxi3_dt;   // OH consumed
    rates.O += dxi3_dt;    // O produced
    rates.H += dxi3_dt;    // H produced
  }

  // Reaction 4: N2 ↔ 2N
  if (xi.length > 4 && reactionDeriv.dXidt.length > 4) {
    const dxi4_dt = reactionDeriv.dXidt[4];
    rates.N2 -= dxi4_dt;   // N2 consumed
    rates.N += 2.0 * dxi4_dt; // 2N produced
  }

  // Reaction 5: H + OH → H2O (three-body)
  if (xi.length > 5 && reactionDeriv.dXidt.length > 5) {
    const dxi5_dt = reactionDeriv.dXidt[5];
    rates.H -= dxi5_dt;    // H consumed
    rates.OH -= dxi5_dt;   // OH consumed
    rates.H2O += dxi5_dt;  // H2O produced
  }

  return rates;
}
