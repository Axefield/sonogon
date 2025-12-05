"use strict";
// simulation/runner.ts
// High-level simulation runner API
Object.defineProperty(exports, "__esModule", { value: true });
exports.runSimulation = runSimulation;
exports.quickSimulation = quickSimulation;
const statevector_1 = require("../core/statevector");
const integrator_1 = require("../core/integrator");
const observables_1 = require("../analysis/observables");
const diagnostics_1 = require("../analysis/diagnostics");
/**
 * Run a complete sonoluminescence simulation
 *
 * This is a high-level API that:
 * 1. Integrates the ODE system
 * 2. Converts results to structured states
 * 3. Performs requested analysis
 * 4. Returns comprehensive results
 */
function runSimulation(config) {
    const { model, initialState, mapper = new statevector_1.DefaultStateVectorMapper(), integratorOptions, events, analysis = {}, progressCallback, } = config;
    // Convert initial state to vector
    const x0 = mapper.toVector(initialState);
    // Create progress-aware RHS function
    let lastProgressTime = 0;
    const rhs = (t, x) => {
        if (progressCallback && t - lastProgressTime > integratorOptions.tMax * 0.01) {
            const progress = t / integratorOptions.tMax;
            progressCallback(progress, t);
            lastProgressTime = t;
        }
        return model.rhs(t, x);
    };
    // Integrate
    const timeSeries = (0, integrator_1.integrateAdaptive)(rhs, x0, integratorOptions, events);
    // Convert results to structured states
    const states = timeSeries.x.map((vec, i) => mapper.fromVector(vec, timeSeries.t[i]));
    // Perform analysis if requested
    const result = {
        timeSeries,
        states,
    };
    // Explicit time series logging (R, Pg, T, ne, Te, E_em, totalPower)
    if (config.logTimeSeries) {
        result.timeSeriesLog = states.map((state, i) => {
            const emission = (0, observables_1.estimateEmission)(state, false);
            const t = timeSeries.t[i];
            // Compute dPg/dt if we have previous state
            let dPg_dt;
            if (i > 0) {
                const dt = timeSeries.t[i] - timeSeries.t[i - 1];
                const prevState = states[i - 1];
                dPg_dt = Math.abs((state.gas.Pg - prevState.gas.Pg) / dt);
            }
            return {
                t,
                R: state.hydro.R,
                Pg: state.gas.Pg,
                T: state.gas.T,
                ne: state.plasma.ne,
                Te: state.plasma.Te,
                E_em: state.em.storedEnergy,
                totalPower: emission.totalPower,
                Rdot: state.hydro.Rdot,
                dPg_dt,
            };
        });
    }
    if (analysis.computeEmission || analysis.computeSpectrum) {
        result.analysis = result.analysis || {};
        result.analysis.emissions = states.map((state) => (0, observables_1.estimateEmission)(state, analysis.computeSpectrum || false));
    }
    if (analysis.computeGradients) {
        result.analysis = result.analysis || {};
        result.analysis.gradients = [];
        for (let i = 1; i < states.length; i++) {
            const dt = timeSeries.t[i] - timeSeries.t[i - 1];
            result.analysis.gradients.push((0, observables_1.computeGradients)(states[i - 1], states[i], dt));
        }
    }
    if (analysis.computeEnergyBudget) {
        result.analysis = result.analysis || {};
        result.analysis.energyBudgets = [];
        // Estimate acoustic power (simplified)
        const acousticPower = 1e3; // 1 kW (typical)
        for (let i = 1; i < states.length; i++) {
            const dt = timeSeries.t[i] - timeSeries.t[i - 1];
            result.analysis.energyBudgets.push((0, diagnostics_1.computeEnergyBudget)(states[i - 1], states[i], dt, acousticPower));
        }
    }
    if (analysis.computePlasmaDiagnostics) {
        result.analysis = result.analysis || {};
        result.analysis.plasmaDiagnostics = states.map((state) => (0, diagnostics_1.computePlasmaDiagnostics)(state));
    }
    if (analysis.computeStability) {
        result.analysis = result.analysis || {};
        result.analysis.stability = (0, diagnostics_1.computeStabilityMetrics)(states, timeSeries.t);
    }
    // Convert events if present
    if (timeSeries.events) {
        result.events = timeSeries.events.map((event) => ({
            time: event.time,
            state: mapper.fromVector(event.state, event.time),
            eventIndex: event.eventIndex,
        }));
    }
    return result;
}
/**
 * Quick simulation helper for common use case
 */
function quickSimulation(model, initialState, tMax, dt = 1e-9) {
    return runSimulation({
        model,
        initialState,
        integratorOptions: {
            dt,
            tMax,
            adaptive: true,
            tolerance: 1e-6,
        },
        analysis: {
            computeEmission: true,
            computeGradients: true,
        },
    });
}
//# sourceMappingURL=runner.js.map