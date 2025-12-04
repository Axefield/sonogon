"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeAcousticState = computeAcousticState;
/**
 * Compute acoustic state with support for:
 * - Single frequency (backward compatible)
 * - Multi-frequency driving
 * - Non-sinusoidal waveforms
 */
function computeAcousticState(t, acoustic, params) {
    // Determine primary frequency for phase evolution
    let primaryOmega = params.omega || 0;
    if (params.frequencies && params.frequencies.length > 0) {
        primaryOmega = params.frequencies[0].frequency;
    }
    const dPhaseDt = primaryOmega;
    // Compute acoustic pressure
    let Pacoustic = 0;
    if (params.frequencies && params.frequencies.length > 0) {
        // Multi-frequency mode
        for (const freq of params.frequencies) {
            const phase = freq.frequency * t + freq.phase;
            Pacoustic -= freq.amplitude * Math.sin(phase);
        }
    }
    else if (params.Pa !== undefined && params.omega !== undefined) {
        // Single frequency mode (backward compatible)
        Pacoustic = -params.Pa * Math.sin(acoustic.phase);
    }
    // Apply waveform transformation if specified
    if (params.waveform && params.waveform !== 'sinusoidal') {
        const normalizedPhase = (acoustic.phase % (2 * Math.PI)) / (2 * Math.PI);
        switch (params.waveform) {
            case 'square':
                Pacoustic = Pacoustic > 0 ? Math.abs(Pacoustic) : -Math.abs(Pacoustic);
                break;
            case 'sawtooth':
                Pacoustic *= 2 * (normalizedPhase - 0.5);
                break;
            case 'custom':
                if (params.customWaveform) {
                    Pacoustic = params.customWaveform(acoustic.phase);
                }
                break;
        }
    }
    // Compute gradients if enabled
    let gradient;
    let laplacian;
    if (params.enableGradients) {
        // For a plane wave: P = P₀ * sin(k·r - ωt)
        // Gradient: ∇P = P₀ * k * cos(k·r - ωt)
        // Laplacian: ∇²P = -P₀ * k² * sin(k·r - ωt)
        const k = params.waveVector || { x: 0, y: 0, z: 0 };
        const k_mag = Math.sqrt(k.x * k.x + k.y * k.y + k.z * k.z);
        if (k_mag > 0) {
            // Compute phase at bubble position (if translation enabled)
            // For now, assume bubble at origin
            const bubblePos = { x: 0, y: 0, z: 0 }; // TODO: get from state.translation
            const k_dot_r = k.x * bubblePos.x + k.y * bubblePos.y + k.z * bubblePos.z;
            const phase_spatial = k_dot_r - primaryOmega * t;
            // Gradient
            const cos_phase = Math.cos(phase_spatial);
            if (params.frequencies && params.frequencies.length > 0) {
                const P0 = params.frequencies[0].amplitude;
                gradient = {
                    x: P0 * k.x * cos_phase,
                    y: P0 * k.y * cos_phase,
                    z: P0 * k.z * cos_phase,
                };
                // Laplacian
                const sin_phase = Math.sin(phase_spatial);
                laplacian = -P0 * k_mag * k_mag * sin_phase;
            }
            else if (params.Pa !== undefined) {
                const cos_phase_local = Math.cos(acoustic.phase);
                gradient = {
                    x: params.Pa * k.x * cos_phase_local,
                    y: params.Pa * k.y * cos_phase_local,
                    z: params.Pa * k.z * cos_phase_local,
                };
                laplacian = -params.Pa * k_mag * k_mag * Math.sin(acoustic.phase);
            }
        }
        // Standing wave pattern (if specified)
        if (params.standingWave && params.nodePosition) {
            const node = params.nodePosition;
            // Standing wave: P = P₀ * sin(k·r) * sin(ωt)
            // Gradient: ∇P = P₀ * k * cos(k·r) * sin(ωt)
            // For simplicity, use simplified form
            if (gradient) {
                const standing_factor = Math.sin(primaryOmega * t);
                gradient.x *= standing_factor;
                gradient.y *= standing_factor;
                gradient.z *= standing_factor;
            }
        }
    }
    return { dPhaseDt, Pacoustic, gradient, laplacian };
}
//# sourceMappingURL=acoustic.js.map