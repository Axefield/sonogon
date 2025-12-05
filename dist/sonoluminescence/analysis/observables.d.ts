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
 * EXPLICITLY COMBINES THREE COMPONENTS:
 *
 * 1. BLACKBODY-ISH T TERM:
 *    P_bb = σ * T⁴ * A
 *    Where σ is Stefan-Boltzmann constant, T is gas temperature, A is bubble surface area.
 *    This represents thermal radiation from the hot gas.
 *
 * 2. BREMSSTRAHLUNG FROM ne, Te:
 *    P_brems ~ ne² * sqrt(Te) * V
 *    Where ne is electron density, Te is electron temperature, V is bubble volume.
 *    This represents free-free radiation from electron-ion collisions in the plasma.
 *
 * 3. EM DECAY FROM E_em:
 *    P_em = E_em / τ
 *    Where E_em is stored EM energy in the negative-space squeezed state, τ is decay time (~1 ns).
 *    This represents photons emitted from the decay of the squeezed vacuum state created
 *    by parametric amplification during collapse. This is the "light from negative cavity state."
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