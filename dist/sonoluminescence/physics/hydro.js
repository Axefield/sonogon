"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeHydroDerivatives = computeHydroDerivatives;
/**
 * Compute hydrodynamic derivatives using Rayleigh-Plesset or Keller-Miksis equation
 *
 * Rayleigh-Plesset (default):
 * R * Rddot + (3/2) * Rdot² = (1/rho) * (Pg - P0 - Pacoustic - Pv - 2*sigma/R - 4*mu*Rdot/R)
 *
 * Keller-Miksis (enhanced, accounts for liquid compressibility):
 * (1 - Rdot/c) * R * Rddot + (3/2) * (1 - Rdot/(3*c)) * Rdot²
 *   = (1/rho) * (1 + Rdot/c + R/c * d/dt) * (Pg - P0 - Pacoustic - Pv - 2*sigma/R - 4*mu*Rdot/R)
 *
 * The Keller-Miksis equation is more accurate for high-speed collapse (Rdot comparable to sound speed).
 */
function computeHydroDerivatives(state, params, Pacoustic, gamma // Adiabatic index from thermo params (optional for backward compatibility)
) {
    const { R, Rdot } = state.hydro;
    const { Pg } = state.gas;
    const { rho, mu, sigma, Pv, P0, c = 1482, useKellerMiksis = false, enableShapeRadialCoupling = false, shapeCouplingCoefficient = 0.1, } = params;
    // Avoid division by zero for very small R
    const R_safe = Math.max(R, 1e-10);
    // Shape-radial coupling: shape oscillations affect effective radius
    // Effective radius: R_eff = R * (1 + (a₂² + a₄²) / (2*R²))
    let R_effective = R_safe;
    if (enableShapeRadialCoupling && state.shape) {
        const { a2, a4 } = state.shape;
        const shapeCorrection = (a2 * a2 + a4 * a4) / (2.0 * R_safe * R_safe);
        R_effective = R_safe * (1.0 + shapeCorrection);
    }
    // dR/dt = Rdot (kinematics)
    // With shape coupling: radial velocity is modified by shape deformation rate
    let dRdt = Rdot;
    if (enableShapeRadialCoupling && state.shape) {
        const { a2_dot, a4_dot } = state.shape;
        // Shape deformation contributes to effective radial motion
        const shapeContribution = shapeCouplingCoefficient * (a2_dot + a4_dot) / R_safe;
        dRdt += shapeContribution;
    }
    // Pressure terms
    // Use effective radius if shape coupling is enabled
    const R_for_pressure = enableShapeRadialCoupling ? R_effective : R_safe;
    const pressureTerm = Pg - P0 - Pacoustic - Pv;
    const surfaceTensionTerm = 2 * sigma / R_for_pressure;
    const viscousTerm = 4 * mu * Rdot / R_for_pressure;
    // Additional pressure term from shape oscillations
    // Shape oscillations create additional surface area and pressure
    let shapePressureTerm = 0;
    if (enableShapeRadialCoupling && state.shape) {
        const { a2, a2_dot, a4, a4_dot } = state.shape;
        // Shape oscillations create restoring forces
        // Additional pressure ~ sigma * (a₂² + a₄²) / R³
        const shapePressure = shapeCouplingCoefficient * sigma * (a2 * a2 + a4 * a4) / (R_safe * R_safe * R_safe);
        shapePressureTerm = shapePressure;
    }
    const netPressure = pressureTerm - surfaceTensionTerm - viscousTerm - shapePressureTerm;
    if (useKellerMiksis && c > 0) {
        // Keller-Miksis equation (accounts for liquid compressibility)
        // More accurate for supersonic collapse (Rdot ~ c)
        // Mach number
        const M = Rdot / c;
        // Compressibility correction factors
        const factor1 = 1.0 - M; // (1 - Rdot/c)
        const factor2 = 1.0 - M / 3.0; // (1 - Rdot/(3*c))
        // Time derivative of pressure term (simplified: assume dPg/dt dominates)
        // dPg/dt from adiabatic compression: -gamma * Pg * (3*Rdot/R)
        // Use provided gamma or fallback to 1.4 (monatomic/diatomic average)
        const effectiveGamma = gamma || 1.4;
        const dPgdt = -effectiveGamma * Pg * (3.0 * Rdot / R_safe);
        const dPdt_term = (R_safe / c) * dPgdt;
        // Modified pressure term with compressibility
        const modifiedPressure = (1.0 + M + dPdt_term) * netPressure;
        // Keller-Miksis: (1 - Rdot/c) * R * Rddot + (3/2) * (1 - Rdot/(3*c)) * Rdot² = (1/rho) * modifiedPressure
        // Solving for Rddot:
        const dRdotdt = (1.0 / (factor1 * R_safe)) * ((1.0 / rho) * modifiedPressure - (3.0 / 2.0) * factor2 * Rdot * Rdot);
        return { dRdt, dRdotdt };
    }
    else {
        // Standard Rayleigh-Plesset equation
        const dRdotdt = (1.0 / R_safe) * ((1.0 / rho) * netPressure - (3.0 / 2.0) * Rdot * Rdot);
        return { dRdt, dRdotdt };
    }
}
//# sourceMappingURL=hydro.js.map