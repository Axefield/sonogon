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
    atomicDisturbances?: ReturnType<typeof computeAtomicSubatomicDisturbances>;
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
 * Compute atomic/subatomic level disturbance metrics
 *
 * At extreme gradients during collapse, conditions approach regimes where:
 * - Nuclear densities are reached (~10¹⁸ kg/m³)
 * - Temperatures reach MeV scales (nuclear binding energies)
 * - Strong force effects may become relevant
 * - Atomic structure is probed at fundamental levels
 *
 * This function tracks conditions where exotic hadron-like states might form
 * or where quark-level physics might manifest, similar to the tetraquark
 * discoveries at CERN where extreme conditions reveal multi-quark bound states.
 *
 * The connection: extreme gradients in sonoluminescence create transient
 * conditions where atomic-level disturbances probe deeper into subatomic structure.
 *
 * EXPANDED VERSION: Now includes detailed nuclear physics, quark deconfinement
 * conditions, exotic hadron formation probabilities, comparison to accelerator conditions,
 * quark flavor tracking, quantum field effects, nuclear fusion conditions, momentum space
 * distributions, and Planck-scale physics.
 */
export declare function computeAtomicSubatomicDisturbances(state: BubbleFullState): {
    nuclearDensityReached: boolean;
    nuclearDensityRatio: number;
    massDensity: number;
    mevTemperatureReached: boolean;
    temperatureMeV: number;
    temperatureRatioToDeconfinement: number;
    strongForceRelevant: boolean;
    strongForceFieldEstimate: number;
    strongCouplingConstant: number;
    strongForceRange: number;
    quarkLevelConditions: boolean;
    approachingQuarkDeconfinement: boolean;
    deconfinementProximity: number;
    exoticStateFormationPossible: boolean;
    tetraquarkFormationProbability: number;
    pentaquarkFormationProbability: number;
    hexaquarkFormationProbability: number;
    hybridMesonFormationProbability: number;
    glueballFormationProbability: number;
    multiQuarkStateFormationRate: number;
    quarkFlavorRelevant: boolean;
    lightQuarkEnergy: number;
    charmQuarkThreshold: boolean;
    strangeQuarkThreshold: boolean;
    heavyQuarkProductionPossible: boolean;
    qcdVacuumEnergy: number;
    vacuumFluctuationsRelevant: boolean;
    casimirEffectEstimate: number;
    gluonFieldStrength: number;
    colorChargeDensity: number;
    fusionConditionsApproached: boolean;
    fusionCrossSection: number;
    gamowFactor: number;
    fusionRate: number;
    fermiEnergy: number;
    degeneracyPressure: number;
    quantumStatisticsRelevant: boolean;
    fermiMomentum: number;
    degeneracyParameter: number;
    approachingPlanckScale: boolean;
    planckDensityRatio: number;
    planckLengthRatio: number;
    planckTimeRatio: number;
    quantumGravityEffectsPossible: boolean;
    nuclearStructureAffected: boolean;
    nuclearBindingEnergyRatio: number;
    nuclearRadius: number;
    comparableToLHC: boolean;
    lhcEnergyRatio: number;
    acceleratorConditions: {
        comparableToLHC: boolean;
        comparableToRHIC: boolean;
        comparableToSPS: boolean;
    };
    qcdPhase: 'hadronic' | 'quark-gluon-plasma' | 'crossover' | 'unknown';
    qcdPhaseDiagramPosition: {
        temperatureRatio: number;
        densityRatio: number;
    };
    atomicStructureProbing: boolean;
    probingDepth: 'atomic' | 'nuclear' | 'quark' | 'planck' | 'none';
    nuclearTimeScale: number;
    strongInteractionTimeScale: number;
    fusionTimeScale: number;
    planckTimeScale: number;
};
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