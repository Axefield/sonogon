import { AcousticState } from "../model/types";
export interface AcousticParams {
    Pa?: number;
    omega?: number;
    frequencies?: Array<{
        amplitude: number;
        frequency: number;
        phase: number;
    }>;
    waveform?: 'sinusoidal' | 'square' | 'sawtooth' | 'custom';
    customWaveform?: (phase: number) => number;
    enableGradients?: boolean;
    waveVector?: {
        x: number;
        y: number;
        z: number;
    };
    standingWave?: boolean;
    nodePosition?: {
        x: number;
        y: number;
        z: number;
    };
}
export interface AcousticDerivatives {
    dPhaseDt: number;
    Pacoustic: number;
    gradient?: {
        x: number;
        y: number;
        z: number;
    };
    laplacian?: number;
}
/**
 * Compute acoustic state with support for:
 * - Single frequency (backward compatible)
 * - Multi-frequency driving
 * - Non-sinusoidal waveforms
 * - Spatial gradients using actual bubble position
 */
export declare function computeAcousticState(t: number, acoustic: AcousticState, params: AcousticParams, bubblePosition?: {
    x: number;
    y: number;
    z: number;
}): AcousticDerivatives;
//# sourceMappingURL=acoustic.d.ts.map