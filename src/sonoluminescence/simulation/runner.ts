// simulation/runner.ts
// High-level simulation runner API

import { SonoluminescenceModel } from "../model/sonoluminescenceModel";
import { StateVectorMapper, DefaultStateVectorMapper } from "../core/statevector";
import { BubbleFullState } from "../model/types";
import {
  integrateAdaptive,
  IntegratorOptions,
  EventFunction,
  IntegrationResult,
} from "../core/integrator";
import {
  estimateEmission,
  computeGradients,
  EmissionSnapshot,
  GradientMetrics,
} from "../analysis/observables";
import {
  computeEnergyBudget,
  computePlasmaDiagnostics,
  computeStabilityMetrics,
  EnergyBudget,
  PlasmaDiagnostics,
  StabilityMetrics,
} from "../analysis/diagnostics";

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
  logTimeSeries?: boolean; // Enable explicit time series logging (R, Pg, T, ne, Te, E_em, totalPower)
  progressCallback?: (progress: number, t: number) => void;
}

export interface TimeSeriesLog {
  t: number;
  R: number;
  Pg: number;
  T: number;
  ne: number;
  Te: number;
  E_em: number;
  totalPower: number;
  Rdot: number;
  dPg_dt?: number; // Optional: pressure derivative magnitude
}

export interface SimulationResult {
  timeSeries: IntegrationResult;
  states: BubbleFullState[];
  timeSeriesLog?: TimeSeriesLog[]; // Explicit logging of key variables
  analysis?: {
    emissions?: EmissionSnapshot[];
    gradients?: GradientMetrics[];
    energyBudgets?: EnergyBudget[];
    plasmaDiagnostics?: PlasmaDiagnostics[];
    stability?: StabilityMetrics;
  };
  events?: Array<{ time: number; state: BubbleFullState; eventIndex: number }>;
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
export function runSimulation(config: SimulationConfig): SimulationResult {
  const {
    model,
    initialState,
    mapper = new DefaultStateVectorMapper(),
    integratorOptions,
    events,
    analysis = {},
    progressCallback,
  } = config;

  // Convert initial state to vector
  const x0 = mapper.toVector(initialState);

  // Create progress-aware RHS function
  let lastProgressTime = 0;
  const rhs = (t: number, x: Float64Array): Float64Array => {
    if (progressCallback && t - lastProgressTime > integratorOptions.tMax * 0.01) {
      const progress = t / integratorOptions.tMax;
      progressCallback(progress, t);
      lastProgressTime = t;
    }
    return model.rhs(t, x);
  };

  // Integrate
  const timeSeries = integrateAdaptive(rhs, x0, integratorOptions, events);

  // Convert results to structured states
  const states = timeSeries.x.map((vec, i) =>
    mapper.fromVector(vec, timeSeries.t[i])
  );

  // Perform analysis if requested
  const result: SimulationResult = {
    timeSeries,
    states,
  };

  // Explicit time series logging (R, Pg, T, ne, Te, E_em, totalPower)
  if (config.logTimeSeries) {
    result.timeSeriesLog = states.map((state, i) => {
      const emission = estimateEmission(state, false);
      const t = timeSeries.t[i];
      
      // Compute dPg/dt if we have previous state
      let dPg_dt: number | undefined;
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
    result.analysis.emissions = states.map((state) =>
      estimateEmission(state, analysis.computeSpectrum || false)
    );
  }

  if (analysis.computeGradients) {
    result.analysis = result.analysis || {};
    result.analysis.gradients = [];
    for (let i = 1; i < states.length; i++) {
      const dt = timeSeries.t[i] - timeSeries.t[i - 1];
      result.analysis.gradients.push(
        computeGradients(states[i - 1], states[i], dt)
      );
    }
  }

  if (analysis.computeEnergyBudget) {
    result.analysis = result.analysis || {};
    result.analysis.energyBudgets = [];
    // Estimate acoustic power (simplified)
    const acousticPower = 1e3; // 1 kW (typical)
    for (let i = 1; i < states.length; i++) {
      const dt = timeSeries.t[i] - timeSeries.t[i - 1];
      result.analysis.energyBudgets.push(
        computeEnergyBudget(states[i - 1], states[i], dt, acousticPower)
      );
    }
  }

  if (analysis.computePlasmaDiagnostics) {
    result.analysis = result.analysis || {};
    result.analysis.plasmaDiagnostics = states.map((state) =>
      computePlasmaDiagnostics(state)
    );
  }

  if (analysis.computeStability) {
    result.analysis = result.analysis || {};
    result.analysis.stability = computeStabilityMetrics(states, timeSeries.t);
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
export function quickSimulation(
  model: SonoluminescenceModel,
  initialState: BubbleFullState,
  tMax: number,
  dt: number = 1e-9
): SimulationResult {
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

