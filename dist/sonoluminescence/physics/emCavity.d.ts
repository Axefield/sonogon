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
    thermoGamma?: number;
    useRefractiveIndexFrequencyShift?: boolean;
    useRefractiveIndexTimeDerivative?: boolean;
    useDynamicQ?: boolean;
    useRadiationBackreaction?: boolean;
    radiationBackreactionCoeff?: number;
    includeBremsstrahlungEmission?: boolean;
    includeRecombinationEmission?: boolean;
    useMagneticFieldEffects?: boolean;
    useFaradayRotation?: boolean;
    staticMagneticField?: number;
    magneticFieldCoupling?: number;
    magneticContributionFactor?: number;
}
export interface EmDerivatives {
    dModes: {
        dRe: number;
        dIm: number;
    }[];
    dStoredEnergyDt: number;
    magneticFieldAmplitude?: number[];
    faradayRotationAngle?: number[];
    magneticEnergy?: number;
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
export interface EmDerivativesWithState {
    dModes: {
        dRe: number;
        dIm: number;
    }[];
    dStoredEnergyDt: number;
    radiationBackreaction?: {
        fx: number;
        fy: number;
        fz: number;
    };
    dynamicQ?: number[];
}
/**
 * Compute EM mode + stored energy derivatives (CAVITY-QED FORMALISM)
 *
 * Enhanced with:
 * - Refractive index frequency shift: ω_k(R(t), n(r,t))
 * - Mode squeezing with ∂n/∂t: ȧk = f(ak, R, Ṙ, ∂n/∂t)
 * - Dynamic Q-factor & cavity losses
 * - Radiation backreaction
 * - Bremsstrahlung + recombination emission
 */
export declare function computeEmDerivatives(state: BubbleFullState, params: EmCavityParams & {
    thermoGamma?: number;
}, statePrev?: BubbleFullState, // Previous state for computing time derivatives
dt?: number): EmDerivatives | EmDerivativesWithState;
//# sourceMappingURL=emCavity.d.ts.map