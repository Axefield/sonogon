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
 * Preset: Pistol Shrimp (Alpheidae)
 *
 * The pistol shrimp creates cavitation bubbles through rapid claw closure,
 * ejecting a high-speed water jet (30-60 m/s) that forms a cavitation bubble.
 * The bubble collapses with extreme violence, generating:
 * - Temperatures up to 8,000°F (4,427°C) - comparable to sun's surface
 * - Sound intensity up to 218 dB (louder than gunshot/jet engine)
 * - Sonoluminescence (light flash)
 * - Shockwaves that stun or kill prey
 *
 * Typical conditions:
 * - Initial bubble radius: 1-5 mm
 * - Jet velocity: 30-60 m/s
 * - Jet duration: ~0.1-1 ms
 * - Collapse velocity: can exceed sound speed (supersonic)
 */
export declare function createPistolShrimpPreset(): SonoluminescenceParams;
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