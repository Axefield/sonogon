import { BubbleFullState } from "../model/types";
export interface BubbleTranslationParams {
    rho: number;
    mu: number;
    enableTranslation?: boolean;
    dragCoeff?: number;
    acousticGradient?: {
        x: number;
        y: number;
        z: number;
    };
    acousticLaplacian?: number;
}
export interface BubbleTranslationDerivatives {
    dx_dt: number;
    dy_dt: number;
    dz_dt: number;
    dvx_dt: number;
    dvy_dt: number;
    dvz_dt: number;
}
/**
 * Compute bubble translation derivatives
 *
 * Position: dx/dt = vx, dy/dt = vy, dz/dt = vz
 *
 * Velocity: dv/dt = F_total / m
 *
 * Forces:
 * 1. Primary Bjerknes force: F = -V * ∇P_acoustic
 * 2. Drag force: F_drag = -6π*μ*R*v (Stokes)
 * 3. Buoyancy (if gravity included)
 */
export declare function computeBubbleTranslationDerivatives(state: BubbleFullState, params: BubbleTranslationParams): BubbleTranslationDerivatives;
/**
 * Compute secondary Bjerknes force between two bubbles
 *
 * Secondary Bjerknes force: F = -ρ * V₁ * V₂ * ∇²P / (4π * r²)
 *
 * This is the interaction force between bubbles.
 * For now, we implement the force on a single bubble from the field.
 */
export declare function computeSecondaryBjerknesForce(V: number, // Bubble volume [m³]
V_other: number, // Other bubble volume [m³]
r: number, // Distance between bubbles [m]
laplacianP: number, // ∇²P_acoustic [Pa/m²]
rho: number): number;
//# sourceMappingURL=bubbleTranslation.d.ts.map