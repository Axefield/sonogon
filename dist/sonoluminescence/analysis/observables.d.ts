import { BubbleFullState } from "../model/types";
export interface EmissionSnapshot {
    totalPower: number;
    blackbodyPower: number;
    bremsstrahlungPower: number;
    emDecayPower: number;
    spectrum?: SpectralDistribution;
}
export interface SpectralDistribution {
    wavelengths: number[];
    powerDensity: number[];
    blackbody: number[];
    bremsstrahlung: number[];
    emModes: number[];
}
export interface GradientMetrics {
    dRdt: number;
    dRdotdt: number;
    dPgdt: number;
    dTdt: number;
    dnedt: number;
    maxGradient: number;
    extremeGradient: boolean;
}
/**
 * Compute spectral emission distribution
 *
 * Returns power per unit wavelength across specified wavelength range
 */
export declare function computeSpectralDistribution(state: BubbleFullState, wavelengthMin?: number, // 100 nm (UV)
wavelengthMax?: number, // 1000 nm (near-IR)
nBins?: number): SpectralDistribution;
/**
 * Estimate total emission power from the bubble
 *
 * Combines:
 * - Blackbody radiation from gas temperature
 * - Bremsstrahlung from plasma (ne, Te)
 * - EM stored energy decay (photons from cavity modes)
 *
 * Optionally computes spectral distribution if computeSpectrum is true
 */
export declare function estimateEmission(state: BubbleFullState, computeSpectrum?: boolean): EmissionSnapshot;
/**
 * Compute gradient metrics between two states
 *
 * Calculates time derivatives and detects extreme gradients
 * that indicate non-classical behavior.
 */
export declare function computeGradients(statePrev: BubbleFullState, stateNext: BubbleFullState, dt: number): GradientMetrics;
/**
 * Compute spatial gradient estimates (simplified)
 *
 * Estimates spatial gradients like |∇P|, |∇n| using
 * bubble geometry and state variables.
 */
export declare function computeSpatialGradients(state: BubbleFullState): {
    pressureGradient: number;
    densityGradient: number;
    temperatureGradient: number;
};
//# sourceMappingURL=observables.d.ts.map