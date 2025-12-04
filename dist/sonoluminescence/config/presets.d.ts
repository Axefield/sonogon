import { SonoluminescenceParams } from "../model/sonoluminescenceModel";
/**
 * Preset: Water with dissolved air (typical experimental setup)
 */
export declare function createWaterAirPreset(): SonoluminescenceParams;
/**
 * Preset: Argon bubble in water (common for sonoluminescence)
 */
export declare function createArgonBubblePreset(): SonoluminescenceParams;
/**
 * Preset: Xenon bubble (high light yield)
 */
export declare function createXenonBubblePreset(): SonoluminescenceParams;
/**
 * Preset: High-intensity driving (extreme conditions)
 */
export declare function createHighIntensityPreset(): SonoluminescenceParams;
/**
 * Validate parameter set for physical reasonableness
 */
export declare function validateParams(params: SonoluminescenceParams): {
    valid: boolean;
    errors: string[];
};
/**
 * Get default preset (water with air)
 */
export declare function getDefaultPreset(): SonoluminescenceParams;
//# sourceMappingURL=presets.d.ts.map