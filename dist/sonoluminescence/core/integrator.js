"use strict";
// core/integrator.ts
// Advanced ODE integrators for sonoluminescence with extreme gradients
Object.defineProperty(exports, "__esModule", { value: true });
exports.integrateRK4 = integrateRK4;
exports.integrateDormandPrince = integrateDormandPrince;
exports.integrateAdaptive = integrateAdaptive;
exports.integrate = integrate;
/**
 * Simple fixed-step RK4 integrator
 * Use for quick simulations or when adaptive control is not needed
 */
function integrateRK4(rhs, x0, options) {
    const { dt, tMax } = options;
    const times = [];
    const states = [];
    let t = 0;
    let x = new Float64Array(x0);
    let steps = 0;
    while (t < tMax) {
        times.push(t);
        states.push(new Float64Array(x));
        const k1 = rhs(t, x);
        const k2 = rhs(t + dt / 2, addVec(x, scaleVec(k1, dt / 2)));
        const k3 = rhs(t + dt / 2, addVec(x, scaleVec(k2, dt / 2)));
        const k4 = rhs(t + dt, addVec(x, scaleVec(k3, dt)));
        const incr = combineRK4(k1, k2, k3, k4, dt);
        x = new Float64Array(addVec(x, incr));
        t += dt;
        steps++;
        // Ensure we don't overshoot tMax
        if (t + dt > tMax) {
            const finalDt = tMax - t;
            if (finalDt > 1e-15) {
                const k1_final = rhs(t, x);
                const k2_final = rhs(t + finalDt / 2, addVec(x, scaleVec(k1_final, finalDt / 2)));
                const k3_final = rhs(t + finalDt / 2, addVec(x, scaleVec(k2_final, finalDt / 2)));
                const k4_final = rhs(t + finalDt, addVec(x, scaleVec(k3_final, finalDt)));
                const incr_final = combineRK4(k1_final, k2_final, k3_final, k4_final, finalDt);
                x = new Float64Array(addVec(x, incr_final));
                t = tMax;
                times.push(t);
                states.push(new Float64Array(x));
            }
            break;
        }
    }
    return {
        t: times,
        x: states,
        stats: {
            steps,
            rejectedSteps: 0,
            minDt: dt,
            maxDt: dt,
        },
    };
}
/**
 * Dormand-Prince 5(4) embedded Runge-Kutta method with adaptive step size
 *
 * This is a high-order method (5th order) with 4th order error estimation.
 * Automatically adjusts step size based on local truncation error.
 *
 * Ideal for sonoluminescence simulations with extreme gradients near collapse.
 */
function integrateDormandPrince(rhs, x0, options) {
    const { dt: dt0, tMax, dtMin = 1e-15, dtMax = 1e-6, tolerance = 1e-6, maxSteps = 1e6, } = options;
    const times = [];
    const states = [];
    let t = 0;
    let x = new Float64Array(x0);
    let dt = dt0;
    let steps = 0;
    let rejectedSteps = 0;
    let minDt = dt;
    let maxDt = dt;
    // Dormand-Prince 5(4) coefficients
    const c = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
    const a = [
        [],
        [1 / 5],
        [3 / 40, 9 / 40],
        [44 / 45, -56 / 15, 32 / 9],
        [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
        [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
        [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
    ];
    const b5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0]; // 5th order
    const b4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40]; // 4th order
    const k = [];
    while (t < tMax && steps < maxSteps) {
        // Store state
        times.push(t);
        states.push(new Float64Array(x));
        // Compute stages
        k[0] = new Float64Array(rhs(t, x));
        for (let i = 1; i < 7; i++) {
            let sum = new Float64Array(x.length);
            for (let j = 0; j < i; j++) {
                sum = new Float64Array(addVec(sum, scaleVec(k[j], a[i][j])));
            }
            k[i] = new Float64Array(rhs(t + c[i] * dt, new Float64Array(addVec(x, scaleVec(sum, dt)))));
        }
        // Compute 5th order solution
        let x5 = new Float64Array(x.length);
        for (let i = 0; i < 7; i++) {
            x5 = new Float64Array(addVec(x5, scaleVec(k[i], b5[i])));
        }
        x5 = new Float64Array(addVec(x, scaleVec(x5, dt)));
        // Compute 4th order solution (for error estimation)
        let x4 = new Float64Array(x.length);
        for (let i = 0; i < 7; i++) {
            x4 = new Float64Array(addVec(x4, scaleVec(k[i], b4[i])));
        }
        x4 = new Float64Array(addVec(x, scaleVec(x4, dt)));
        // Estimate error
        const error = estimateError(x4, x5, tolerance);
        const errorNorm = Math.sqrt(error.reduce((sum, e, i) => {
            const scale = Math.max(Math.abs(x[i]), 1.0);
            return sum + (e / scale) * (e / scale);
        }, 0) / x.length);
        // Accept or reject step
        if (errorNorm <= 1.0 || dt <= dtMin) {
            // Accept step
            x = x5; // Use 5th order solution
            t += dt;
            steps++;
            // Update step size statistics
            minDt = Math.min(minDt, dt);
            maxDt = Math.max(maxDt, dt);
            // Adjust step size for next step
            if (errorNorm > 0) {
                const safety = 0.9;
                const factor = safety * Math.pow(1.0 / errorNorm, 0.2);
                dt = Math.max(dtMin, Math.min(dtMax, dt * factor));
            }
            else {
                dt = Math.min(dtMax, dt * 1.5); // Increase if error is zero
            }
        }
        else {
            // Reject step, reduce dt
            rejectedSteps++;
            const safety = 0.9;
            const factor = safety * Math.pow(1.0 / errorNorm, 0.25);
            dt = Math.max(dtMin, dt * factor);
        }
        // Ensure we don't overshoot tMax
        if (t + dt > tMax) {
            dt = tMax - t;
            if (dt < dtMin)
                break;
        }
    }
    return {
        t: times,
        x: states,
        stats: {
            steps,
            rejectedSteps,
            minDt,
            maxDt,
        },
    };
}
/**
 * Adaptive integrator with event detection
 *
 * Detects events (e.g., minimum radius, extreme gradients) and
 * records them for analysis.
 */
function integrateAdaptive(rhs, x0, options, events) {
    const useAdaptive = options.adaptive !== false;
    if (useAdaptive) {
        // Use Dormand-Prince for adaptive integration
        const result = integrateDormandPrince(rhs, x0, options);
        // Detect events if provided
        if (events && events.length > 0) {
            const detectedEvents = [];
            for (let i = 1; i < result.t.length; i++) {
                const t_prev = result.t[i - 1];
                const t_curr = result.t[i];
                const x_prev = result.x[i - 1];
                const x_curr = result.x[i];
                // Check each event function
                for (let j = 0; j < events.length; j++) {
                    let val_prev = events[j](t_prev, x_prev);
                    let val_curr = events[j](t_curr, x_curr);
                    // Sign change indicates event crossing
                    if (val_prev * val_curr < 0) {
                        // Simple bisection to find event time (could be improved)
                        let t_event = (t_prev + t_curr) / 2;
                        let x_event = interpolateState(x_prev, x_curr, t_prev, t_curr, t_event);
                        let val_event = events[j](t_event, x_event);
                        // Refine event location
                        let t_low = t_prev;
                        let t_high = t_curr;
                        let x_low = new Float64Array(x_prev);
                        let x_high = new Float64Array(x_curr);
                        let val_low = val_prev;
                        let val_high = val_curr;
                        for (let iter = 0; iter < 10 && Math.abs(val_event) > 1e-10; iter++) {
                            if (val_low * val_event < 0) {
                                t_high = t_event;
                                x_high = new Float64Array(x_event);
                                val_high = val_event;
                            }
                            else {
                                t_low = t_event;
                                x_low = new Float64Array(x_event);
                                val_low = val_event;
                            }
                            t_event = (t_low + t_high) / 2;
                            x_event = interpolateState(x_low, x_high, t_low, t_high, t_event);
                            val_event = events[j](t_event, x_event);
                        }
                        detectedEvents.push({
                            time: t_event,
                            state: x_event,
                            eventIndex: j,
                        });
                    }
                }
            }
            return {
                ...result,
                events: detectedEvents,
            };
        }
        return result;
    }
    else {
        // Use fixed-step RK4
        return integrateRK4(rhs, x0, options);
    }
}
// Helper functions
function addVec(a, b) {
    const out = new Float64Array(a.length);
    for (let i = 0; i < a.length; i++)
        out[i] = a[i] + b[i];
    return out;
}
function scaleVec(a, s) {
    const out = new Float64Array(a.length);
    for (let i = 0; i < a.length; i++)
        out[i] = a[i] * s;
    return out;
}
function combineRK4(k1, k2, k3, k4, dt) {
    const out = new Float64Array(k1.length);
    for (let i = 0; i < k1.length; i++) {
        out[i] = (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
    return out;
}
function estimateError(x4, x5, tolerance) {
    const error = new Float64Array(x4.length);
    for (let i = 0; i < x4.length; i++) {
        error[i] = Math.abs(x5[i] - x4[i]);
    }
    return error;
}
function interpolateState(x1, x2, t1, t2, t) {
    const alpha = (t - t1) / (t2 - t1);
    const out = new Float64Array(x1.length);
    for (let i = 0; i < x1.length; i++) {
        out[i] = x1[i] + alpha * (x2[i] - x1[i]);
    }
    return out;
}
/**
 * Convenience function: automatically chooses best integrator
 */
function integrate(rhs, x0, options, events) {
    return integrateAdaptive(rhs, x0, options, events);
}
//# sourceMappingURL=integrator.js.map