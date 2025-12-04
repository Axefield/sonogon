import { BubbleFullState } from "../model/types";
export interface EmCavityParams {
    modeFrequencies0: number[];
    pumpCoefficient: number;
    decayTime: number;
    couplingStrength: number;
    refractiveIndexBase: number;
    useModeCoupling?: boolean;
    modeCouplingMatrix?: number[][];
    useFrequencyDependentQ?: boolean;
    Q0?: number;
    QFrequencyDependence?: (omega: number) => number;
    useDetailedParametricPumping?: boolean;
    parametricCoupling?: number;
}
export interface EmDerivatives {
    dModes: {
        dRe: number;
        dIm: number;
    }[];
    dStoredEnergyDt: number;
}
/**
 * Compute EM mode + stored energy derivatives.
 *
 * Mode evolution:
 * - Parametric amplification from boundary motion (Rdot pumps modes)
 * - Mode frequency changes with R(t) and plasma frequency
 * - Damping from cavity losses
 *
 * Stored energy:
 * - Pump term: α * |Rdot| * |gradient_terms|
 * - Decay term: -E_em / τ
 */
export declare function computeEmDerivatives(state: BubbleFullState, params: EmCavityParams): EmDerivatives;
//# sourceMappingURL=emCavity.d.ts.map