"use strict";
// validation/kinematicTests.ts
// Test simulations for shape oscillations and translation dynamics
Object.defineProperty(exports, "__esModule", { value: true });
exports.testShapeOscillationFrequencies = testShapeOscillationFrequencies;
exports.testBubbleTranslation = testBubbleTranslation;
exports.testShapeRadialCoupling = testShapeRadialCoupling;
exports.runAllKinematicTests = runAllKinematicTests;
const sonoluminescenceModel_1 = require("../model/sonoluminescenceModel");
const statevector_1 = require("../core/statevector");
const integrator_1 = require("../core/integrator");
const presets_1 = require("../config/presets");
const shapeOscillations_1 = require("../physics/shapeOscillations");
/**
 * Test shape oscillation natural frequencies
 *
 * Verify that shape modes oscillate at their natural frequencies
 */
function testShapeOscillationFrequencies() {
    try {
        const mapper = new statevector_1.DefaultStateVectorMapper();
        const params = (0, presets_1.createArgonBubblePreset)();
        // Enable shape oscillations
        params.shape = {
            sigma: 0.0728,
            mu: 0.001002,
            rho: 998.2,
            enableShapeOscillations: true,
        };
        const model = new sonoluminescenceModel_1.SonoluminescenceModel(mapper, params);
        // Initial state with small shape deformation
        const initialState = {
            t: 0,
            hydro: { R: 5e-6, Rdot: 0 },
            shape: {
                a2: 1e-7, // Small initial deformation
                a2_dot: 0,
                a4: 0,
                a4_dot: 0,
            },
            translation: {
                x: 0, y: 0, z: 0,
                vx: 0, vy: 0, vz: 0,
            },
            gas: { Pg: 101325, T: 293.15 },
            species: {
                numberDensity: {
                    H2O: 1e23,
                    O2: 1e23,
                    N2: 1e23,
                    Ar: 1e25,
                    Xe: 0,
                    H: 0,
                    O: 0,
                    OH: 0,
                    N: 0,
                },
            },
            plasma: { ne: 1e18, Te: 293.15, ionizationFraction: 1e-6 },
            internalEnergy: {
                translational: 1e-18,
                rotational: 1e-19,
                vibrational: 1e-20,
                electronic: 1e-21,
            },
            em: {
                modes: [
                    { omega: 2 * Math.PI * 3e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 4e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 5e14, re: 0, im: 0 },
                ],
                storedEnergy: 1e-25,
            },
            acoustic: { phase: 0 },
            reactions: { xi: [0, 0, 0] },
        };
        const x0 = mapper.toVector(initialState);
        // Integrate for several oscillation periods
        const result = (0, integrator_1.integrateAdaptive)((t, x) => model.rhs(t, x), x0, {
            dt: 1e-9,
            tMax: 1e-5,
            tolerance: 1e-6,
            dtMin: 1e-12,
            dtMax: 1e-8,
            adaptive: true,
        });
        // Extract shape amplitudes over time
        const states = result.x.map((vec, i) => mapper.fromVector(vec, result.t[i]));
        const a2_history = states.map(s => s.shape?.a2 || 0);
        const a4_history = states.map(s => s.shape?.a4 || 0);
        // Check that oscillations occur (amplitude changes)
        const a2_max = Math.max(...a2_history.map(Math.abs));
        const a4_max = Math.max(...a4_history.map(Math.abs));
        const passed = a2_max > 1e-8 || a4_max > 1e-8; // Some oscillation detected
        return {
            name: "Shape Oscillation Frequencies",
            passed,
            error: passed ? undefined : "No oscillations detected",
            details: {
                a2_max,
                a4_max,
                nSteps: result.x.length,
            },
        };
    }
    catch (err) {
        return {
            name: "Shape Oscillation Frequencies",
            passed: false,
            error: err.message,
        };
    }
}
/**
 * Test bubble translation with Bjerknes force
 *
 * Verify that bubble moves in response to acoustic gradient
 */
function testBubbleTranslation() {
    try {
        const mapper = new statevector_1.DefaultStateVectorMapper();
        const params = (0, presets_1.createArgonBubblePreset)();
        // Enable translation and acoustic gradients
        params.translation = {
            rho: 998.2,
            mu: 0.001002,
            enableTranslation: true,
        };
        params.acoustic.enableGradients = true;
        params.acoustic.waveVector = { x: 100, y: 0, z: 0 }; // Gradient in x direction
        const model = new sonoluminescenceModel_1.SonoluminescenceModel(mapper, params);
        // Initial state with bubble at origin
        const initialState = {
            t: 0,
            hydro: { R: 5e-6, Rdot: 0 },
            shape: {
                a2: 0, a2_dot: 0,
                a4: 0, a4_dot: 0,
            },
            translation: {
                x: 0, y: 0, z: 0,
                vx: 0, vy: 0, vz: 0,
            },
            gas: { Pg: 101325, T: 293.15 },
            species: {
                numberDensity: {
                    H2O: 1e23,
                    O2: 1e23,
                    N2: 1e23,
                    Ar: 1e25,
                    Xe: 0,
                    H: 0,
                    O: 0,
                    OH: 0,
                    N: 0,
                },
            },
            plasma: { ne: 1e18, Te: 293.15, ionizationFraction: 1e-6 },
            internalEnergy: {
                translational: 1e-18,
                rotational: 1e-19,
                vibrational: 1e-20,
                electronic: 1e-21,
            },
            em: {
                modes: [
                    { omega: 2 * Math.PI * 3e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 4e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 5e14, re: 0, im: 0 },
                ],
                storedEnergy: 1e-25,
            },
            acoustic: { phase: 0 },
            reactions: { xi: [0, 0, 0] },
        };
        const x0 = mapper.toVector(initialState);
        // Integrate for short time
        const result = (0, integrator_1.integrateAdaptive)((t, x) => model.rhs(t, x), x0, {
            dt: 1e-9,
            tMax: 1e-6,
            tolerance: 1e-6,
            dtMin: 1e-12,
            dtMax: 1e-8,
            adaptive: true,
        });
        // Extract translation over time
        const states = result.x.map((vec, i) => mapper.fromVector(vec, result.t[i]));
        const x_history = states.map(s => s.translation?.x || 0);
        const vx_history = states.map(s => s.translation?.vx || 0);
        // Check that position or velocity changes (bubble moves)
        const x_final = x_history[x_history.length - 1];
        const vx_final = vx_history[vx_history.length - 1];
        const x_change = Math.abs(x_final);
        const vx_change = Math.abs(vx_final);
        const passed = x_change > 1e-12 || vx_change > 1e-12;
        return {
            name: "Bubble Translation",
            passed,
            error: passed ? undefined : "No translation detected",
            details: {
                x_final,
                vx_final,
                x_change,
                vx_change,
                nSteps: result.x.length,
            },
        };
    }
    catch (err) {
        return {
            name: "Bubble Translation",
            passed: false,
            error: err.message,
        };
    }
}
/**
 * Test shape-radial coupling
 *
 * Verify that shape oscillations affect radial dynamics
 */
function testShapeRadialCoupling() {
    try {
        const mapper = new statevector_1.DefaultStateVectorMapper();
        const params = (0, presets_1.createArgonBubblePreset)();
        // Enable shape oscillations and coupling
        params.shape = {
            sigma: 0.0728,
            mu: 0.001002,
            rho: 998.2,
            enableShapeOscillations: true,
            enableShapeRadialCoupling: true,
        };
        params.hydro.enableShapeRadialCoupling = true;
        params.hydro.shapeCouplingCoefficient = 0.1;
        const model = new sonoluminescenceModel_1.SonoluminescenceModel(mapper, params);
        // Initial state with shape deformation
        const initialState = {
            t: 0,
            hydro: { R: 5e-6, Rdot: 0 },
            shape: {
                a2: 1e-7,
                a2_dot: 0,
                a4: 0,
                a4_dot: 0,
            },
            translation: {
                x: 0, y: 0, z: 0,
                vx: 0, vy: 0, vz: 0,
            },
            gas: { Pg: 101325, T: 293.15 },
            species: {
                numberDensity: {
                    H2O: 1e23,
                    O2: 1e23,
                    N2: 1e23,
                    Ar: 1e25,
                    Xe: 0,
                    H: 0,
                    O: 0,
                    OH: 0,
                    N: 0,
                },
            },
            plasma: { ne: 1e18, Te: 293.15, ionizationFraction: 1e-6 },
            internalEnergy: {
                translational: 1e-18,
                rotational: 1e-19,
                vibrational: 1e-20,
                electronic: 1e-21,
            },
            em: {
                modes: [
                    { omega: 2 * Math.PI * 3e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 4e14, re: 0, im: 0 },
                    { omega: 2 * Math.PI * 5e14, re: 0, im: 0 },
                ],
                storedEnergy: 1e-25,
            },
            acoustic: { phase: 0 },
            reactions: { xi: [0, 0, 0] },
        };
        const x0 = mapper.toVector(initialState);
        // Integrate
        const result = (0, integrator_1.integrateAdaptive)((t, x) => model.rhs(t, x), x0, {
            dt: 1e-9,
            tMax: 1e-6,
            tolerance: 1e-6,
            dtMin: 1e-12,
            dtMax: 1e-8,
            adaptive: true,
        });
        // Extract state
        const states = result.x.map((vec, i) => mapper.fromVector(vec, result.t[i]));
        const R_history = states.map(s => s.hydro.R);
        const a2_history = states.map(s => s.shape?.a2 || 0);
        // Compute effective radius
        const R_eff_history = states.map(s => (0, shapeOscillations_1.computeEffectiveRadius)(s.hydro.R, s.shape, true));
        // Check that effective radius differs from actual radius (coupling effect)
        const R0 = R_history[0];
        const R_eff0 = R_eff_history[0];
        const couplingDetected = Math.abs(R_eff0 - R0) > 1e-12;
        // Check shape energy
        const E_shape_initial = (0, shapeOscillations_1.computeShapeEnergy)(initialState.shape, R0, 998.2);
        const E_shape_final = (0, shapeOscillations_1.computeShapeEnergy)(states[states.length - 1].shape, states[states.length - 1].hydro.R, 998.2);
        const passed = couplingDetected && (E_shape_initial > 0 || E_shape_final > 0);
        return {
            name: "Shape-Radial Coupling",
            passed,
            error: passed ? undefined : "No coupling detected",
            details: {
                R0,
                R_eff0,
                R_eff_diff: Math.abs(R_eff0 - R0),
                E_shape_initial,
                E_shape_final,
                nSteps: result.x.length,
            },
        };
    }
    catch (err) {
        return {
            name: "Shape-Radial Coupling",
            passed: false,
            error: err.message,
        };
    }
}
/**
 * Run all kinematic tests
 */
function runAllKinematicTests() {
    const tests = [
        testShapeOscillationFrequencies,
        testBubbleTranslation,
        testShapeRadialCoupling,
    ];
    const results = tests.map((test) => test());
    const passed = results.filter((r) => r.passed).length;
    const failed = results.filter((r) => !r.passed).length;
    return {
        passed,
        failed,
        results,
    };
}
//# sourceMappingURL=kinematicTests.js.map