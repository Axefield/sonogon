export type RhsFunction = (t: number, x: Float64Array) => Float64Array;
export type EventFunction = (t: number, x: Float64Array) => number;
export interface IntegratorOptions {
    dt: number;
    tMax: number;
    dtMin?: number;
    dtMax?: number;
    tolerance?: number;
    adaptive?: boolean;
    maxSteps?: number;
}
export interface IntegrationResult {
    t: number[];
    x: Float64Array[];
    events?: Array<{
        time: number;
        state: Float64Array;
        eventIndex: number;
    }>;
    stats: {
        steps: number;
        rejectedSteps: number;
        minDt: number;
        maxDt: number;
    };
}
/**
 * Simple fixed-step RK4 integrator
 * Use for quick simulations or when adaptive control is not needed
 */
export declare function integrateRK4(rhs: RhsFunction, x0: Float64Array, options: IntegratorOptions): IntegrationResult;
/**
 * Dormand-Prince 5(4) embedded Runge-Kutta method with adaptive step size
 *
 * This is a high-order method (5th order) with 4th order error estimation.
 * Automatically adjusts step size based on local truncation error.
 *
 * Ideal for sonoluminescence simulations with extreme gradients near collapse.
 */
export declare function integrateDormandPrince(rhs: RhsFunction, x0: Float64Array, options: IntegratorOptions): IntegrationResult;
/**
 * Adaptive integrator with event detection
 *
 * Detects events (e.g., minimum radius, extreme gradients) and
 * records them for analysis.
 */
export declare function integrateAdaptive(rhs: RhsFunction, x0: Float64Array, options: IntegratorOptions, events?: EventFunction[]): IntegrationResult;
/**
 * Convenience function: automatically chooses best integrator
 */
export declare function integrate(rhs: RhsFunction, x0: Float64Array, options: IntegratorOptions, events?: EventFunction[]): IntegrationResult;
//# sourceMappingURL=integrator.d.ts.map