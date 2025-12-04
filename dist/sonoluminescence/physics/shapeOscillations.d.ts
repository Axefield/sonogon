import { BubbleFullState, ShapeOscillationState } from "../model/types";
export interface ShapeOscillationParams {
    sigma: number;
    mu: number;
    rho: number;
    enableShapeOscillations?: boolean;
    omega2_natural?: number;
    omega4_natural?: number;
    damping2?: number;
    damping4?: number;
    enableNonlinearCoupling?: boolean;
    enableShapeRadialCoupling?: boolean;
    enableAcousticCoupling?: boolean;
    couplingCoefficient2?: number;
    couplingCoefficient4?: number;
    nonlinearCoefficient?: number;
}
export interface ShapeOscillationDerivatives {
    da2_dt: number;
    da2dot_dt: number;
    da4_dt: number;
    da4dot_dt: number;
}
/**
 * Compute shape oscillation derivatives
 *
 * Shape modes follow a damped harmonic oscillator:
 *
 * d²a_n/dt² + 2*ζ*ω_n*da_n/dt + ω_n²*a_n = F_n
 *
 * where:
 * - ω_n is the natural frequency
 * - ζ is the damping ratio
 * - F_n is the driving force (from radial motion, acoustic field, etc.)
 *
 * The driving force comes from:
 * 1. Radial motion coupling (Rdot affects shape)
 * 2. Acoustic field gradients
 * 3. Pressure variations
 */
export declare function computeShapeOscillationDerivatives(state: BubbleFullState, params: ShapeOscillationParams): ShapeOscillationDerivatives;
/**
 * Compute effective radius accounting for shape deformations
 *
 * For shape R(θ) = R₀ + a₂*P₂(cos θ) + a₄*P₄(cos θ)
 * The effective radius (volume-equivalent) is approximately:
 * R_eff ≈ R₀ * (1 + (a₂² + a₄²) / (2*R₀²))
 *
 * Higher-order correction includes mode interactions
 */
export declare function computeEffectiveRadius(R0: number, shape?: ShapeOscillationState, includeNonlinear?: boolean): number;
/**
 * Compute energy in shape oscillations
 *
 * E_shape = (1/2) * m * (a₂_dot² + a₄_dot²) + (1/2) * k * (a₂² + a₄²)
 * where m is effective mass and k is spring constant
 */
export declare function computeShapeEnergy(shape: ShapeOscillationState, R: number, rho: number): number;
//# sourceMappingURL=shapeOscillations.d.ts.map