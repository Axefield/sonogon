"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeEmDerivatives = computeEmDerivatives;
const units_1 = require("../core/units");
/**
 * Compute mode frequency as function of bubble radius
 *
 * For a spherical cavity, mode frequencies scale with 1/R
 * ω_k(t) = ω_k0 * (R0 / R(t))
 *
 * Also affected by plasma frequency: ω_eff = sqrt(ω_k² - ω_p²)
 */
function computeModeFrequency(omega0, R, R0, omega_p) {
    // Geometric scaling: frequency ~ 1/R for cavity modes
    const omega_geometric = omega0 * (R0 / Math.max(R, 1e-10));
    // Plasma frequency effect: ω_eff = sqrt(ω_k² - ω_p²) if ω_k > ω_p
    // If ω_p > ω_k, mode is cut off
    if (omega_p >= omega_geometric) {
        return 0; // Mode cut off by plasma
    }
    return Math.sqrt(omega_geometric * omega_geometric - omega_p * omega_p);
}
/**
 * Compute EM mode + stored energy derivatives.
 *
 * Mode evolution:
 * - Parametric amplification from boundary motion (Rdot pumps modes)
 * - Mode frequency changes with R(t) and plasma frequency
 * - Damping from cavity losses
 *
 * Stored energy:
 * - Pump term: α * |Rdot| * |gradient_terms|
 * - Decay term: -E_em / τ
 */
function computeEmDerivatives(state, params) {
    const { R, Rdot } = state.hydro;
    const { Pg } = state.gas;
    const { ne } = state.plasma;
    const modes = state.em.modes;
    const nModes = modes.length;
    const R0 = params.modeFrequencies0.length > 0 ? 1e-6 : 1e-6; // Reference radius
    // Compute plasma frequency
    const omega_p = units_1.Calculations.plasmaFrequency(ne);
    // Compute gradient magnitudes for pumping
    // More detailed calculation of parametric coupling strength
    const Rdot_magnitude = Math.abs(Rdot);
    // Pressure gradient from radial motion
    // dP/dR ≈ -gamma * P * (3/R) for adiabatic compression
    const gamma = params.thermoGamma || 1.4;
    const Pg_gradient_radial = gamma * Math.abs(Pg) / Math.max(R, 1e-10);
    // Pressure gradient from time variation
    // dP/dt ≈ -gamma * P * (3*Rdot/R) for adiabatic
    const Pg_gradient_temporal = gamma * Math.abs(Pg) * Math.abs(Rdot) / Math.max(R, 1e-10);
    // Combined pressure gradient (geometric mean for coupling strength)
    const Pg_gradient = Math.sqrt(Pg_gradient_radial * Pg_gradient_temporal);
    // Acoustic pressure gradient (if available from acoustic module)
    // This would come from acoustic field gradients
    const acousticGradient = 0; // TODO: pass from acoustic module if available
    // Combined gradient metric for parametric pumping
    // Parametric coupling strength scales with gradient magnitude
    const gradientMagnitude = Rdot_magnitude + Pg_gradient * 1e-6 + acousticGradient * 1e-8;
    const dModes = Array.from({ length: nModes }, (_, k) => {
        if (k >= params.modeFrequencies0.length) {
            return { dRe: 0, dIm: 0 };
        }
        const omega0 = params.modeFrequencies0[k];
        const omega_k = computeModeFrequency(omega0, R, R0, omega_p);
        // If mode is cut off, decay to zero
        if (omega_k === 0) {
            const decayRate = 1e12; // fast decay when cut off
            return {
                dRe: -state.em.modes[k].re * decayRate,
                dIm: -state.em.modes[k].im * decayRate,
            };
        }
        const a_re = state.em.modes[k].re;
        const a_im = state.em.modes[k].im;
        // Mode evolution equations
        // da/dt = -i*ω*a + pump_term - damping_term
        // Pumping term: parametric amplification from boundary motion
        let pumpRe = 0;
        let pumpIm = 0;
        if (params.useDetailedParametricPumping) {
            // Detailed parametric pumping Hamiltonian
            // H = g(t) * (a†² + a²) where g(t) = g0 * Rdot/R
            // This gives parametric amplification when Rdot < 0 (compression)
            // 
            // NEGATIVE-SPACE BEHAVIOR:
            // The parametric Hamiltonian H = g(t) * (a†² + a²) creates a squeezed
            // vacuum state. In the Wigner function representation, this appears as
            // a "negative-space" region where the Wigner function becomes negative,
            // indicating non-classical behavior. The squeezing occurs in quadrature
            // space: one quadrature is squeezed (reduced uncertainty) while the
            // conjugate quadrature is anti-squeezed (increased uncertainty).
            //
            // During extreme compression (large |Rdot|), g(t) becomes large and
            // the modes are strongly squeezed, storing energy in the negative-space
            // state. This energy is then released as photons when the cavity decays.
            //
            // From Heisenberg equations: da/dt = -i*[a, H] = -i*g(t)*(a† + a)
            // In terms of real/imaginary parts:
            // dRe/dt = g(t) * Im
            // dIm/dt = -g(t) * Re
            const g0 = params.parametricCoupling || params.couplingStrength;
            const R_safe = Math.max(R, 1e-10);
            const g_t = g0 * (Rdot / R_safe); // Time-dependent coupling
            // Parametric pumping (quadrature-dependent)
            // This cross-coupling (Re ↔ Im) is the signature of parametric amplification
            // and creates the squeezed/negative-space state
            pumpRe = g_t * a_im;
            pumpIm = -g_t * a_re;
            // Additional squeezing term (higher order)
            // For strong coupling, add squeezing: H_squeeze = g_squeeze * (a†² + a²)
            // This enhances the negative-space behavior for very strong compression
            const g_squeeze = g0 * 0.1; // Squeezing is weaker
            const squeezeRe = g_squeeze * a_re; // Squeezing in phase
            const squeezeIm = g_squeeze * a_im;
            pumpRe += squeezeRe;
            pumpIm += squeezeIm;
        }
        else {
            // Simplified pumping (original)
            // Pump strength ~ couplingStrength * Rdot * gradient
            const pumpStrength = params.couplingStrength * gradientMagnitude;
            pumpRe = pumpStrength * a_im; // Cross-coupling for parametric amplification
            pumpIm = -pumpStrength * a_re;
        }
        // Damping: cavity losses (quality factor Q)
        let Q = params.Q0 || 1000; // Base quality factor
        if (params.useFrequencyDependentQ) {
            if (params.QFrequencyDependence) {
                Q = params.QFrequencyDependence(omega_k);
            }
            else {
                // Default frequency dependence: Q decreases with frequency
                // Q(ω) = Q0 * (ω0/ω)^alpha
                const omega_ref = params.modeFrequencies0[0] || 2 * Math.PI * 3e14;
                const alpha = 0.5; // Power law exponent
                Q = Q * Math.pow(omega_ref / Math.max(omega_k, omega_ref * 0.1), alpha);
            }
        }
        const damping = omega_k / (2.0 * Q);
        // Cross-mode coupling (if enabled)
        let couplingRe = 0;
        let couplingIm = 0;
        if (params.useModeCoupling && params.modeCouplingMatrix) {
            const couplingMatrix = params.modeCouplingMatrix;
            // Coupling term: sum_j g_kj * a_j
            for (let j = 0; j < nModes && j < couplingMatrix.length; j++) {
                if (j !== k && couplingMatrix[k] && couplingMatrix[k][j] !== undefined) {
                    const g_kj = couplingMatrix[k][j];
                    // Coupling adds to both real and imaginary parts
                    couplingRe += g_kj * state.em.modes[j].re;
                    couplingIm += g_kj * state.em.modes[j].im;
                }
            }
        }
        // Mode evolution: da/dt = -i*ω*a + pump - damping*a + coupling
        // Real part: dRe/dt = ω*Im + pumpRe - damping*Re + couplingRe
        // Imag part: dIm/dt = -ω*Re + pumpIm - damping*Im + couplingIm
        const dRe = omega_k * a_im + pumpRe - damping * a_re + couplingRe;
        const dIm = -omega_k * a_re + pumpIm - damping * a_im + couplingIm;
        return { dRe, dIm };
    });
    // Stored energy evolution
    // dE_em/dt = pump_term - decay_term
    // 
    // Pump term: α * |Rdot| * |gradient_terms|
    // This represents parametric amplification from boundary motion.
    // During extreme compression (large |Rdot|, large |dPg/dt|), the boundary
    // motion pumps energy into EM modes, creating a "negative-space" squeezed state.
    // The negative-space refers to the quantum squeezed vacuum state in the
    // Wigner function representation, where the uncertainty in one quadrature
    // is reduced below the vacuum level at the expense of increased uncertainty
    // in the conjugate quadrature.
    //
    // Decay term: -E_em / τ
    // The stored energy decays over a short timescale (τ ~ 1 ns), converting
    // to photons (light emission). This decay is what produces the sonoluminescence
    // flash: E_em rises during collapse, then rapidly decays as photons are emitted.
    const pumpTerm = params.pumpCoefficient * gradientMagnitude;
    const decayTerm = state.em.storedEnergy / Math.max(params.decayTime, 1e-12);
    const dStoredEnergyDt = pumpTerm - decayTerm;
    return { dModes, dStoredEnergyDt };
}
//# sourceMappingURL=emCavity.js.map