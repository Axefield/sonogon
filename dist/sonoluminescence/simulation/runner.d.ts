import { SonoluminescenceModel } from "../model/sonoluminescenceModel";
import { StateVectorMapper } from "../core/statevector";
import { BubbleFullState } from "../model/types";
import { IntegratorOptions, EventFunction, IntegrationResult } from "../core/integrator";
import { EmissionSnapshot, GradientMetrics } from "../analysis/observables";
import { EnergyBudget, PlasmaDiagnostics, StabilityMetrics } from "../analysis/diagnostics";
export interface SimulationConfig {
    model: SonoluminescenceModel;
    initialState: BubbleFullState;
    mapper?: StateVectorMapper;
    integratorOptions: IntegratorOptions;
    events?: EventFunction[];
    analysis?: {
        computeEmission?: boolean;
        computeGradients?: boolean;
        computeEnergyBudget?: boolean;
        computePlasmaDiagnostics?: boolean;
        computeStability?: boolean;
        computeSpectrum?: boolean;
    };
    progressCallback?: (progress: number, t: number) => void;
}
export interface SimulationResult {
    timeSeries: IntegrationResult;
    states: BubbleFullState[];
    analysis?: {
        emissions?: EmissionSnapshot[];
        gradients?: GradientMetrics[];
        energyBudgets?: EnergyBudget[];
        plasmaDiagnostics?: PlasmaDiagnostics[];
        stability?: StabilityMetrics;
    };
    events?: Array<{
        time: number;
        state: BubbleFullState;
        eventIndex: number;
    }>;
}
/**
 * Run a complete sonoluminescence simulation
 *
 * This is a high-level API that:
 * 1. Integrates the ODE system
 * 2. Converts results to structured states
 * 3. Performs requested analysis
 * 4. Returns comprehensive results
 */
export declare function runSimulation(config: SimulationConfig): SimulationResult;
/**
 * Quick simulation helper for common use case
 */
export declare function quickSimulation(model: SonoluminescenceModel, initialState: BubbleFullState, tMax: number, dt?: number): SimulationResult;
//# sourceMappingURL=runner.d.ts.map