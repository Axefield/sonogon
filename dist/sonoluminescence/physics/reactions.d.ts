import { BubbleFullState } from "../model/types";
export interface ReactionParams {
    reaction0_A: number;
    reaction0_Ea: number;
    reaction1_A: number;
    reaction1_Ea: number;
    reaction2_A: number;
    reaction2_Ea: number;
    reaction3_A?: number;
    reaction3_Ea?: number;
    reaction4_A?: number;
    reaction4_Ea?: number;
    reaction5_A?: number;
    reaction5_Ea?: number;
    usePressureDependentRates?: boolean;
    reaction0_k0?: number;
    reaction0_kInf?: number;
    reaction0_Fc?: number;
    useDetailedThreeBody?: boolean;
    thirdBodyEfficiencies?: Record<string, number>;
}
export interface ReactionDerivatives {
    dXidt: number[];
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
export declare function computeReactionDerivatives(state: BubbleFullState, params: ReactionParams): ReactionDerivatives;
/**
 * Compute species number density changes from reactions
 *
 * This should be called in addition to computeSpeciesDerivatives from thermoChem
 * to add reaction terms to the species evolution
 */
export declare function computeSpeciesReactionRates(state: BubbleFullState, params: ReactionParams): Record<string, number>;
//# sourceMappingURL=reactions.d.ts.map