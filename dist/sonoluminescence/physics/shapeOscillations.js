"use strict";
// physics/shapeOscillations.ts
// Shape oscillation dynamics using spherical harmonics
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeShapeOscillationDerivatives = computeShapeOscillationDerivatives;
exports.computeEffectiveRadius = computeEffectiveRadius;
exports.computeShapeEnergy = computeShapeEnergy;
/**
 * Compute natural frequency for shape mode n
 *
 * For a spherical bubble, the natural frequency is:
 * ω_n = sqrt((n-1)(n+1)(n+2) * sigma / (rho * R³))
 *
 * This comes from the restoring force due to surface tension.
 */
function computeNaturalFrequency(n, // Mode number (2 for quadrupole, 4 for hexadecapole)
R, // Bubble radius [m]
sigma, // Surface tension [N/m]
rho // Liquid density [kg/m³]
) {
    const R_safe = Math.max(R, 1e-10);
    const prefactor = (n - 1) * (n + 1) * (n + 2);
    return Math.sqrt((prefactor * sigma) / (rho * R_safe * R_safe * R_safe));
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
function computeShapeOscillationDerivatives(state, params) {
    if (!params.enableShapeOscillations || !state.shape) {
        return { da2_dt: 0, da2dot_dt: 0, da4_dt: 0, da4dot_dt: 0 };
    }
    const { R, Rdot } = state.hydro;
    const { a2, a2_dot, a4, a4_dot } = state.shape;
    const { sigma, mu, rho } = params;
    const R_safe = Math.max(R, 1e-10);
    // Natural frequencies
    const omega2 = params.omega2_natural || computeNaturalFrequency(2, R, sigma, rho);
    const omega4 = params.omega4_natural || computeNaturalFrequency(4, R, sigma, rho);
    // Damping (from viscosity)
    // Simplified: ζ ≈ 2*mu / (rho * R² * ω_n)
    const zeta2 = params.damping2 || (2.0 * mu) / (rho * R_safe * R_safe * omega2);
    const zeta4 = params.damping4 || (2.0 * mu) / (rho * R_safe * R_safe * omega4);
    // Driving forces
    let F2 = 0;
    let F4 = 0;
    // 1. Radial motion coupling: shape is driven by radial acceleration
    if (params.enableShapeRadialCoupling !== false) {
        // Linear coupling: F_n ~ Rdot² / R (centrifugal-like effect)
        const coupling2 = params.couplingCoefficient2 || 0.1;
        const coupling4 = params.couplingCoefficient4 || 0.05;
        const F2_radial = coupling2 * (Rdot * Rdot) / R_safe;
        const F4_radial = coupling4 * (Rdot * Rdot) / R_safe;
        // Nonlinear coupling: shape affects radial motion feedback
        // F_n ~ a_n * Rdot (parametric coupling)
        const F2_nonlinear = params.enableNonlinearCoupling
            ? (params.nonlinearCoefficient || 0.01) * a2 * Rdot / R_safe
            : 0;
        const F4_nonlinear = params.enableNonlinearCoupling
            ? (params.nonlinearCoefficient || 0.005) * a4 * Rdot / R_safe
            : 0;
        F2 += F2_radial + F2_nonlinear;
        F4 += F4_radial + F4_nonlinear;
    }
    // 2. Mode-mode interactions (nonlinear coupling)
    if (params.enableNonlinearCoupling) {
        // Mode 2 and Mode 4 can interact
        // F_2 += α * a₂ * a₄ (quadratic coupling)
        // F_4 += β * a₂² (mode 2 drives mode 4)
        const alpha = params.nonlinearCoefficient || 0.01;
        const F2_modeCoupling = alpha * a2 * a4 / R_safe;
        const F4_modeCoupling = alpha * a2 * a2 / R_safe;
        F2 += F2_modeCoupling;
        F4 += F4_modeCoupling;
    }
    // 3. Acoustic field coupling
    if (params.enableAcousticCoupling && state.acoustic) {
        // Shape modes can be driven by acoustic pressure gradients
        // For standing waves, pressure gradients drive shape
        // F_n ~ ∇P_acoustic (simplified)
        // In a full implementation, would use acoustic gradient from acoustic module
        const acousticCoupling = 1e-6; // Small coupling
        // Placeholder: would use actual gradient if available
        F2 += acousticCoupling * Math.sin(state.acoustic.phase);
        F4 += acousticCoupling * 0.5 * Math.sin(state.acoustic.phase);
    }
    // 4. Pressure coupling (from gas pressure variations)
    // Shape is affected by pressure asymmetry
    const { Pg } = state.gas;
    const pressureCoupling = 1e-9; // Small effect
    F2 += pressureCoupling * Pg / R_safe;
    F4 += pressureCoupling * 0.5 * Pg / R_safe;
    // Mode 2 (quadrupole) evolution
    // da2/dt = a2_dot
    // da2dot/dt = -2*ζ₂*ω₂*a2_dot - ω₂²*a2 + F2
    // 
    // With nonlinear corrections:
    // Additional terms from mode interactions and energy exchange
    let da2_dt = a2_dot;
    let da2dot_dt = -2.0 * zeta2 * omega2 * a2_dot - omega2 * omega2 * a2 + F2;
    // Nonlinear frequency shift (amplitude-dependent)
    if (params.enableNonlinearCoupling) {
        const nonlinearShift = (params.nonlinearCoefficient || 0.01) * (a2 * a2 + a4 * a4) / (R_safe * R_safe);
        da2dot_dt -= omega2 * omega2 * a2 * nonlinearShift;
    }
    // Mode 4 (hexadecapole) evolution
    // da4/dt = a4_dot
    // da4dot/dt = -2*ζ₄*ω₄*a4_dot - ω₄²*a4 + F4
    let da4_dt = a4_dot;
    let da4dot_dt = -2.0 * zeta4 * omega4 * a4_dot - omega4 * omega4 * a4 + F4;
    // Nonlinear frequency shift for mode 4
    if (params.enableNonlinearCoupling) {
        const nonlinearShift = (params.nonlinearCoefficient || 0.01) * (a2 * a2 + a4 * a4) / (R_safe * R_safe);
        da4dot_dt -= omega4 * omega4 * a4 * nonlinearShift;
    }
    return { da2_dt, da2dot_dt, da4_dt, da4dot_dt };
}
/**
 * Compute effective radius accounting for shape deformations
 *
 * For shape R(θ) = R₀ + a₂*P₂(cos θ) + a₄*P₄(cos θ)
 * The effective radius (volume-equivalent) is approximately:
 * R_eff ≈ R₀ * (1 + (a₂² + a₄²) / (2*R₀²))
 *
 * Higher-order correction includes mode interactions
 */
function computeEffectiveRadius(R0, shape, includeNonlinear) {
    if (!shape)
        return R0;
    const { a2, a4 } = shape;
    const R_safe = Math.max(R0, 1e-10);
    // Volume-conserving correction (linear)
    let correction = (a2 * a2 + a4 * a4) / (2.0 * R_safe * R_safe);
    // Nonlinear correction (mode interactions)
    if (includeNonlinear) {
        const nonlinearCorrection = (a2 * a4) / (4.0 * R_safe * R_safe); // Cross-term
        correction += nonlinearCorrection;
    }
    return R0 * (1.0 + correction);
}
/**
 * Compute energy in shape oscillations
 *
 * E_shape = (1/2) * m * (a₂_dot² + a₄_dot²) + (1/2) * k * (a₂² + a₄²)
 * where m is effective mass and k is spring constant
 */
function computeShapeEnergy(shape, R, rho) {
    const { a2, a2_dot, a4, a4_dot } = shape;
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const m_eff = rho * volume; // Effective mass
    // Kinetic energy
    const E_kin = 0.5 * m_eff * (a2_dot * a2_dot + a4_dot * a4_dot);
    // Potential energy (from surface tension)
    const sigma = 0.0728; // Surface tension [N/m]
    const k2 = m_eff * computeNaturalFrequency(2, R, sigma, rho) * computeNaturalFrequency(2, R, sigma, rho);
    const k4 = m_eff * computeNaturalFrequency(4, R, sigma, rho) * computeNaturalFrequency(4, R, sigma, rho);
    const E_pot = 0.5 * (k2 * a2 * a2 + k4 * a4 * a4);
    return E_kin + E_pot;
}
//# sourceMappingURL=shapeOscillations.js.map