import { BubbleFullState } from "../model/types";
import { SonoluminescenceParams } from "../model/sonoluminescenceModel";
/**
 * Create an equilibrium initial state
 *
 * Bubble at acoustic equilibrium: R = R_eq, Rdot = 0
 * Gas pressure balances acoustic pressure
 */
export declare function createEquilibriumState(params: SonoluminescenceParams, R_equilibrium: number): BubbleFullState;
/**
 * Create initial state at maximum expansion
 *
 * Bubble at maximum radius (typically 10-50x equilibrium)
 * Rdot = 0 (turning point)
 */
export declare function createExpandedState(params: SonoluminescenceParams, R_max: number): BubbleFullState;
/**
 * Create initial state just before collapse
 *
 * Bubble starting to collapse: R = R_min_start, Rdot < 0
 */
export declare function createCollapseState(params: SonoluminescenceParams, R_min_start: number, Rdot_collapse?: number): BubbleFullState;
/**
 * Create initial state from preset with custom radius
 */
export declare function createFromPreset(params: SonoluminescenceParams, R_initial: number, Rdot_initial?: number): BubbleFullState;
/**
 * Validate initial state for physical consistency
 */
export declare function validateInitialState(state: BubbleFullState): {
    valid: boolean;
    errors: string[];
};
//# sourceMappingURL=initialStates.d.ts.map