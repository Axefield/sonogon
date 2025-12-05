import { BubbleFullState } from "../model/types";
import { computeSpatialGradients, GradientMetrics } from "./observables";
export interface EnergyBudget {
    acousticWorkIn: number;
    compressionWork: number;
    heatLoss: number;
    lightEmission: number;
    viscousDissipation: number;
    chemicalEnergy: number;
    emStoredEnergy: number;
    totalEnergy: number;
}
export interface StabilityMetrics {
    isStable: boolean;
    limitCycleDetected: boolean;
    runawayHeating: boolean;
    energyConserved: boolean;
    periodEstimate?: number;
}
export interface RefractiveIndexEstimate {
    n_effective: number;
    n_plasma: number;
    n_neutral: number;
    n_gradient: number;
}
export interface PlasmaDiagnostics {
    plasmaFrequency: number;
    debyeLength: number;
    opticalFrequency: number;
    plasmaCutoff: boolean;
    meanFreePath: number;
    collisional: boolean;
}
/**
 * Compute energy budget for the bubble system
 *
 * Tracks energy flows to detect conservation and identify
 * dominant energy pathways.
 */
export declare function computeEnergyBudget(state: BubbleFullState, statePrev: BubbleFullState, dt: number, acousticPower: number): EnergyBudget;
/**
 * Detect extreme gradient events
 *
 * Identifies when the system enters non-classical regimes
 * with extreme spatial or temporal gradients.
 */
export declare function detectExtremeGradients(statePrev: BubbleFullState, stateNext: BubbleFullState, dt: number): {
    detected: boolean;
    metrics: GradientMetrics;
    spatialGradients: ReturnType<typeof computeSpatialGradients>;
};
/**
 * Estimate refractive index from bubble state
 *
 * Combines contributions from:
 * - Neutral gas (density-dependent)
 * - Plasma (frequency-dependent, can be negative)
 * - Ionization effects
 */
export declare function estimateRefractiveIndex(state: BubbleFullState, opticalFrequency: number): RefractiveIndexEstimate;
/**
 * Compute plasma diagnostics
 *
 * Analyzes plasma properties and their relationship to
 * optical frequencies and bubble size.
 */
export declare function computePlasmaDiagnostics(state: BubbleFullState, opticalFrequency?: number): PlasmaDiagnostics;
/**
 * Compute stability metrics
 *
 * Detects limit cycles, runaway heating, and energy conservation
 * to assess system stability.
 */
export declare function computeStabilityMetrics(states: BubbleFullState[], times: number[]): StabilityMetrics;
/**
 * Compute mode squeezing metrics
 *
 * Analyzes EM mode amplitudes to detect parametric amplification
 * and mode squeezing effects. Reports pumped energy into modes vs released energy.
 *
 * This quantifies the negative-space behavior: energy pumped into squeezed states
 * during compression vs energy released as photons during decay.
 */
export declare function computeModeSqueezingMetrics(state: BubbleFullState, statePrev?: BubbleFullState, dt?: number): {
    modeAmplitudes: number[];
    totalModeEnergy: number;
    squeezingDetected: boolean;
    pumpEfficiency: number;
    pumpedEnergy: number;
    releasedEnergy: number;
    netEnergyChange: number;
};
//# sourceMappingURL=diagnostics.d.ts.map