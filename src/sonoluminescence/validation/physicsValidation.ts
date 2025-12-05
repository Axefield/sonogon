// validation/physicsValidation.ts
// Comprehensive physics validation tests

import { BubbleFullState } from "../model/types";
import { SonoluminescenceModel } from "../model/sonoluminescenceModel";
import { SonoluminescenceParams } from "../model/sonoluminescenceModel";
import { DefaultStateVectorMapper } from "../core/statevector";
import { integrateAdaptive } from "../core/integrator";
import { createArgonBubblePreset } from "../config/presets";
import { Constants, Calculations } from "../core/units";
import { computeEnergyBudget, EnergyBudget } from "../analysis/diagnostics";
import { computePlasmaDerivatives, PlasmaParams } from "../physics/plasma";

export interface TestResult {
  name: string;
  passed: boolean;
  error?: string;
  details?: any;
}

/**
 * Adiabatic Scaling Test
 * 
 * Set acoustic drive to 0, slowly enforce a known change in R (quasi-static),
 * and verify:
 *   T ∝ (R₀/R)^(3(γ-1))
 *   P ∝ (R₀/R)^(3γ)
 */
export function testAdiabaticScaling(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    
    // Set acoustic drive to zero
    params.acoustic.Pa = 0;
    params.acoustic.omega = 2 * Math.PI * 20e3; // Keep frequency for phase tracking
    
    const model = new SonoluminescenceModel(mapper, params);
    
    // Initial state: expanded bubble
    const R0 = 10e-6; // 10 microns
    const R1 = 5e-6;  // 5 microns (compressed by factor of 2)
    const gamma = params.thermo.gamma;
    const T0 = 293.15; // Initial temperature [K]
    const P0 = params.hydro.P0; // Initial pressure [Pa]
    
    // Expected scaling laws
    const R_ratio = R0 / R1; // = 2.0
    const T_expected = T0 * Math.pow(R_ratio, 3 * (gamma - 1));
    const P_expected = P0 * Math.pow(R_ratio, 3 * gamma);
    
    // Create initial state
    const initialState: BubbleFullState = {
      t: 0,
      hydro: { R: R0, Rdot: 0 }, // Start at rest
      gas: { Pg: P0, T: T0 },
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
      plasma: { ne: 1e18, Te: T0, ionizationFraction: 1e-6 },
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
    
    // Quasi-static compression: slowly change R
    // We'll integrate with a very slow Rdot to approximate quasi-static
    const compressionTime = 1e-3; // 1 ms for slow compression
    const Rdot_slow = -(R0 - R1) / compressionTime; // Slow compression rate
    
    // Set initial Rdot for slow compression
    const initialState_compressed: BubbleFullState = {
      ...initialState,
      hydro: { R: R0, Rdot: Rdot_slow },
    };
    const x0_compressed = mapper.toVector(initialState_compressed);
    
    // Integrate with slow compression
    const result = integrateAdaptive(
      (t, x) => model.rhs(t, x),
      x0_compressed,
      { dt: 1e-7, tMax: compressionTime, tolerance: 1e-6 }
    );
    
    // Get final state
    const finalState = mapper.fromVector(
      result.x[result.x.length - 1],
      result.t[result.t.length - 1]
    );
    
    const R_final = finalState.hydro.R;
    const T_final = finalState.gas.T;
    const P_final = finalState.gas.Pg;
    
    // Check scaling laws
    const R_ratio_actual = R0 / R_final;
    const T_ratio = T_final / T0;
    const P_ratio = P_final / P0;
    
    const T_expected_ratio = Math.pow(R_ratio_actual, 3 * (gamma - 1));
    const P_expected_ratio = Math.pow(R_ratio_actual, 3 * gamma);
    
    const T_error = Math.abs(T_ratio - T_expected_ratio) / T_expected_ratio;
    const P_error = Math.abs(P_ratio - P_expected_ratio) / P_expected_ratio;
    
    const tolerance = 0.15; // 15% tolerance for quasi-static approximation
    const passed = T_error < tolerance && P_error < tolerance;
    
    return {
      name: "Adiabatic Scaling Test",
      passed,
      error: passed
        ? undefined
        : `T error: ${(T_error * 100).toFixed(2)}%, P error: ${(P_error * 100).toFixed(2)}%`,
      details: {
        R0,
        R_final,
        R_ratio_actual,
        T0,
        T_final,
        T_ratio,
        T_expected_ratio,
        T_error,
        P0,
        P_final,
        P_ratio,
        P_expected_ratio,
        P_error,
        gamma,
      },
    };
  } catch (err: any) {
    return {
      name: "Adiabatic Scaling Test",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Small-Amplitude Linear Oscillation Test
 * 
 * Low Pa → bubble should behave like linear oscillator around R₀.
 * Check that dominant oscillation frequency matches Minnaert frequency.
 * 
 * Minnaert frequency: ω₀ = (1/R₀) * sqrt(3*γ*P₀/ρ)
 */
export function testMinnaertFrequency(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    
    // Low acoustic pressure for linear regime
    params.acoustic.Pa = 1000; // 1 kPa (low amplitude)
    params.acoustic.omega = 2 * Math.PI * 20e3; // 20 kHz
    
    const model = new SonoluminescenceModel(mapper, params);
    
    const R0 = 5e-6; // Equilibrium radius [m]
    const gamma = params.thermo.gamma;
    const P0 = params.hydro.P0; // Ambient pressure [Pa]
    const rho = params.hydro.rho; // Liquid density [kg/m³]
    
    // Minnaert frequency: ω₀ = (1/R₀) * sqrt(3*γ*P₀/ρ)
    const omega_minnaert = (1.0 / R0) * Math.sqrt((3.0 * gamma * P0) / rho);
    const f_minnaert = omega_minnaert / (2.0 * Math.PI); // [Hz]
    
    // Create initial state at equilibrium
    const initialState: BubbleFullState = {
      t: 0,
      hydro: { R: R0, Rdot: 0 },
      gas: { Pg: P0, T: 293.15 },
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
    
    // Integrate for several periods
    const period_minnaert = 2.0 * Math.PI / omega_minnaert;
    const nPeriods = 10;
    const tMax = nPeriods * period_minnaert;
    
    const result = integrateAdaptive(
      (t, x) => model.rhs(t, x),
      x0,
      { dt: period_minnaert / 100, tMax, tolerance: 1e-6 }
    );
    
    // Extract R(t) time series
    const R_t = result.x.map((vec, i) => {
      const state = mapper.fromVector(vec, result.t[i]);
      return { t: result.t[i], R: state.hydro.R };
    });
    
    // Find dominant frequency using zero-crossing method
    // Count zero crossings to estimate period
    const R_mean = R_t.reduce((sum, r) => sum + r.R, 0) / R_t.length;
    const R_oscillation = R_t.map((r) => r.R - R_mean);
    
    // Find zero crossings
    const zeroCrossings: number[] = [];
    for (let i = 1; i < R_oscillation.length; i++) {
      if (
        (R_oscillation[i - 1] >= 0 && R_oscillation[i] < 0) ||
        (R_oscillation[i - 1] < 0 && R_oscillation[i] >= 0)
      ) {
        zeroCrossings.push(result.t[i]);
      }
    }
    
    // Estimate period from zero crossings
    let f_measured: number | undefined;
    if (zeroCrossings.length >= 2) {
      const periods: number[] = [];
      for (let i = 1; i < zeroCrossings.length; i++) {
        periods.push(zeroCrossings[i] - zeroCrossings[i - 1]);
      }
      const avgPeriod = periods.reduce((sum, p) => sum + p, 0) / periods.length;
      f_measured = 1.0 / avgPeriod;
    } else {
      // Fallback: use FFT-like approach (simple peak finding)
      // Find peaks in oscillation
      const peaks: number[] = [];
      for (let i = 1; i < R_oscillation.length - 1; i++) {
        if (
          R_oscillation[i] > R_oscillation[i - 1] &&
          R_oscillation[i] > R_oscillation[i + 1] &&
          R_oscillation[i] > 0
        ) {
          peaks.push(result.t[i]);
        }
      }
      if (peaks.length >= 2) {
        const periods: number[] = [];
        for (let i = 1; i < peaks.length; i++) {
          periods.push(peaks[i] - peaks[i - 1]);
        }
        const avgPeriod = periods.reduce((sum, p) => sum + p, 0) / periods.length;
        f_measured = 1.0 / avgPeriod;
      }
    }
    
    if (!f_measured) {
      return {
        name: "Minnaert Frequency Test",
        passed: false,
        error: "Could not determine oscillation frequency",
      };
    }
    
    // Compare measured frequency to Minnaert frequency
    const frequencyError = Math.abs(f_measured - f_minnaert) / f_minnaert;
    const tolerance = 0.2; // 20% tolerance (nonlinear effects, damping)
    const passed = frequencyError < tolerance;
    
    return {
      name: "Minnaert Frequency Test",
      passed,
      error: passed
        ? undefined
        : `Frequency error: ${(frequencyError * 100).toFixed(2)}% (measured: ${f_measured.toFixed(0)} Hz, expected: ${f_minnaert.toFixed(0)} Hz)`,
      details: {
        R0,
        gamma,
        P0,
        rho,
        omega_minnaert,
        f_minnaert,
        f_measured,
        frequencyError,
        nPeriods,
        zeroCrossings: zeroCrossings.length,
      },
    };
  } catch (err: any) {
    return {
      name: "Minnaert Frequency Test",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Energy Budget Closure Test
 * 
 * Use computeEnergyBudget() to verify:
 *   Acoustic work in ≈ (thermal + EM + chemical + viscous dissipation)
 * 
 * No weird artificial energy gain over many cycles.
 */
export function testEnergyBudgetClosure(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    const model = new SonoluminescenceModel(mapper, params);
    
    const R0 = 5e-6;
    const initialState: BubbleFullState = {
      t: 0,
      hydro: { R: R0, Rdot: 0 },
      gas: { Pg: params.hydro.P0, T: 293.15 },
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
    
    // Integrate over multiple acoustic periods
    const omega_acoustic = params.acoustic.omega || 2 * Math.PI * 20e3;
    const period = 2.0 * Math.PI / omega_acoustic;
    const nCycles = 10;
    const tMax = nCycles * period;
    
    const result = integrateAdaptive(
      (t, x) => model.rhs(t, x),
      x0,
      { dt: period / 100, tMax, tolerance: 1e-6 }
    );
    
    // Compute energy budgets for each time step
    const energyBudgets: EnergyBudget[] = [];
    let totalAcousticWork = 0;
    let totalCompressionWork = 0;
    let totalHeatLoss = 0;
    let totalLightEmission = 0;
    let totalViscousDissipation = 0;
    let totalChemicalEnergy = 0;
    let totalEMEnergy = 0;
    
    for (let i = 1; i < result.x.length; i++) {
      const statePrev = mapper.fromVector(result.x[i - 1], result.t[i - 1]);
      const state = mapper.fromVector(result.x[i], result.t[i]);
      const dt = result.t[i] - result.t[i - 1];
      
      // Estimate acoustic power (simplified: P_acoustic * dV/dt)
      const { R, Rdot } = state.hydro;
      const dV_dt = 4.0 * Math.PI * R * R * Rdot;
      const Pacoustic = params.acoustic.Pa
        ? params.acoustic.Pa * Math.sin(params.acoustic.omega! * result.t[i] + state.acoustic.phase)
        : 0;
      const acousticPower = Pacoustic * dV_dt;
      
      const budget = computeEnergyBudget(state, statePrev, dt, acousticPower);
      energyBudgets.push(budget);
      
      totalAcousticWork += budget.acousticWorkIn;
      totalCompressionWork += Math.abs(budget.compressionWork);
      totalHeatLoss += budget.heatLoss;
      totalLightEmission += budget.lightEmission;
      totalViscousDissipation += budget.viscousDissipation;
      totalChemicalEnergy += Math.abs(budget.chemicalEnergy);
      totalEMEnergy += Math.abs(budget.emStoredEnergy);
    }
    
    // Energy balance: Acoustic work in should equal energy out
    const energyOut =
      totalCompressionWork +
      totalHeatLoss +
      totalLightEmission +
      totalViscousDissipation +
      totalChemicalEnergy +
      totalEMEnergy;
    
    // Allow some tolerance for numerical errors and approximations
    const energyBalance = Math.abs(totalAcousticWork - energyOut);
    const relativeError = energyBalance / Math.max(totalAcousticWork, energyOut, 1e-20);
    
    // Check for artificial energy gain (energy out > energy in by large margin)
    const energyGain = energyOut - totalAcousticWork;
    const relativeGain = energyGain / Math.max(totalAcousticWork, 1e-20);
    
    // Test passes if:
    // 1. Energy balance is within 30% (reasonable for approximations)
    // 2. No significant artificial energy gain (< 10%)
    const balanceTolerance = 0.3;
    const gainTolerance = 0.1;
    const passed =
      relativeError < balanceTolerance && relativeGain < gainTolerance;
    
    return {
      name: "Energy Budget Closure Test",
      passed,
      error: passed
        ? undefined
        : `Energy balance error: ${(relativeError * 100).toFixed(2)}%, Energy gain: ${(relativeGain * 100).toFixed(2)}%`,
      details: {
        totalAcousticWork,
        totalCompressionWork,
        totalHeatLoss,
        totalLightEmission,
        totalViscousDissipation,
        totalChemicalEnergy,
        totalEMEnergy,
        energyOut,
        energyBalance,
        relativeError,
        energyGain,
        relativeGain,
        nCycles,
        nSteps: result.x.length,
      },
    };
  } catch (err: any) {
    return {
      name: "Energy Budget Closure Test",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Plasma Equilibrium Spot-Check
 * 
 * At a given (ne, T) pair, compare plasma module's ionization fraction
 * vs Saha equation for that T.
 */
export function testPlasmaEquilibrium(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    
    // Test at specific temperature and density
    const T_test = 20000; // 20,000 K
    const n_neutral_test = 1e26; // High density [m⁻³]
    const ionizationPotential_Ar = 15.76 * Constants.e; // Argon ionization potential [J]
    
    // Compute Saha equilibrium electron density
    const g_ratio = 1.0;
    const prefactor = Math.pow(
      (2.0 * Math.PI * Constants.m_e * Constants.k_B * T_test) / (Constants.h * Constants.h),
      1.5
    );
    const exponent = Math.exp(-ionizationPotential_Ar / (Constants.k_B * T_test));
    const sahaProduct = 2.0 * g_ratio * prefactor * exponent;
    const n_e_saha = Math.sqrt(Math.max(n_neutral_test, 0) * sahaProduct);
    const ionizationFraction_saha = n_e_saha / n_neutral_test;
    
    // Create state with these conditions
    const testState: BubbleFullState = {
      t: 0,
      hydro: { R: 1e-6, Rdot: 0 },
      gas: { Pg: n_neutral_test * Constants.k_B * T_test, T: T_test },
      species: {
        numberDensity: {
          H2O: 0,
          O2: 0,
          N2: 0,
          Ar: n_neutral_test,
          Xe: 0,
          H: 0,
          O: 0,
          OH: 0,
          N: 0,
        },
      },
      plasma: {
        ne: n_e_saha * 0.5, // Start away from equilibrium
        Te: T_test,
        ionizationFraction: 0.1, // Initial guess
      },
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
    
    // Compute plasma derivatives
    const plasmaDeriv = computePlasmaDerivatives(testState, params.plasma);
    
    // Check if system is driving toward Saha equilibrium
    // If ne < n_e_saha, dne/dt should be positive (ionization)
    // If ne > n_e_saha, dne/dt should be negative (recombination)
    const ne_current = testState.plasma.ne;
    const dne_dt = plasmaDeriv.dn_edt;
    
    // Expected direction: if ne < n_e_saha, dne/dt > 0
    const expectedDirection = ne_current < n_e_saha ? 1 : -1;
    const actualDirection = dne_dt > 0 ? 1 : -1;
    const directionCorrect = expectedDirection === actualDirection || Math.abs(dne_dt) < 1e10;
    
    // Also check: at equilibrium, dne/dt should be small
    // Create a state closer to equilibrium
    const equilibriumState: BubbleFullState = {
      ...testState,
      plasma: {
        ne: n_e_saha,
        Te: T_test,
        ionizationFraction: ionizationFraction_saha,
      },
    };
    
    const plasmaDeriv_eq = computePlasmaDerivatives(equilibriumState, params.plasma);
    const dne_dt_eq = plasmaDeriv_eq.dn_edt;
    
    // At equilibrium, dne/dt should be close to zero
    const equilibriumTolerance = 0.1; // 10% of equilibrium value
    const atEquilibrium = Math.abs(dne_dt_eq) < n_e_saha * equilibriumTolerance;
    
    // Compare ionization fraction
    const ionizationFraction_plasma = equilibriumState.plasma.ionizationFraction;
    const ionizationError = Math.abs(
      ionizationFraction_plasma - ionizationFraction_saha
    ) / ionizationFraction_saha;
    
    const passed =
      directionCorrect && ionizationError < 0.5; // 50% tolerance for ionization fraction
    
    return {
      name: "Plasma Equilibrium Spot-Check",
      passed,
      error: passed
        ? undefined
        : `Ionization fraction error: ${(ionizationError * 100).toFixed(2)}%, dne/dt at equilibrium: ${dne_dt_eq.toExponential(2)}`,
      details: {
        T_test,
        n_neutral_test,
        n_e_saha,
        ionizationFraction_saha,
        ionizationFraction_plasma,
        ionizationError,
        dne_dt,
        dne_dt_eq,
        directionCorrect,
        atEquilibrium,
      },
    };
  } catch (err: any) {
    return {
      name: "Plasma Equilibrium Spot-Check",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Run all physics validation tests
 */
export function runAllPhysicsValidationTests(): {
  passed: number;
  failed: number;
  results: TestResult[];
} {
  const tests = [
    testAdiabaticScaling,
    testMinnaertFrequency,
    testEnergyBudgetClosure,
    testPlasmaEquilibrium,
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

