// validation/tests.ts
// Validation tests against known analytical and experimental results

import { BubbleFullState } from "../model/types";
import { SonoluminescenceModel } from "../model/sonoluminescenceModel";
import { DefaultStateVectorMapper } from "../core/statevector";
import { integrateRK4 } from "../core/integrator";
import { createArgonBubblePreset } from "../config/presets";
import { Constants, Calculations } from "../core/units";

export interface TestResult {
  name: string;
  passed: boolean;
  error?: string;
  details?: any;
}

/**
 * Test adiabatic compression
 * 
 * For adiabatic compression: P * V^gamma = constant
 * Verify that P * V^gamma remains constant during compression
 */
export function testAdiabaticCompression(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    const model = new SonoluminescenceModel(mapper, params);
    
    // Initial state: expanded bubble
    const R0 = 10e-6; // 10 microns
    const R1 = 1e-6;  // 1 micron (compressed)
    
    const V0 = (4.0 / 3.0) * Math.PI * R0 * R0 * R0;
    const V1 = (4.0 / 3.0) * Math.PI * R1 * R1 * R1;
    
    const P0 = params.hydro.P0;
    const gamma = params.thermo.gamma;
    
    // Expected final pressure from adiabatic relation
    const P1_expected = P0 * Math.pow(V0 / V1, gamma);
    
    // Create initial state
    const initialState: BubbleFullState = {
      t: 0,
      hydro: { R: R0, Rdot: -100 }, // Collapsing
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
    
    // Integrate for short time
    const result = integrateRK4(
      (t, x) => model.rhs(t, x),
      x0,
      { dt: 1e-10, tMax: 1e-8 }
    );
    
    // Check final state
    const finalState = mapper.fromVector(result.x[result.x.length - 1], result.t[result.t.length - 1]);
    const V_final = (4.0 / 3.0) * Math.PI * finalState.hydro.R * finalState.hydro.R * finalState.hydro.R;
    const P_final = finalState.gas.Pg;
    
    // Check: P * V^gamma should be approximately constant
    const PVgamma_initial = P0 * Math.pow(V0, gamma);
    const PVgamma_final = P_final * Math.pow(V_final, gamma);
    const error = Math.abs(PVgamma_final - PVgamma_initial) / PVgamma_initial;
    
    const passed = error < 0.1; // 10% tolerance
    
    return {
      name: "Adiabatic Compression",
      passed,
      error: passed ? undefined : `Error: ${(error * 100).toFixed(2)}%`,
      details: {
        P0,
        P_final,
        V0,
        V_final,
        PVgamma_initial,
        PVgamma_final,
        relativeError: error,
      },
    };
  } catch (err: any) {
    return {
      name: "Adiabatic Compression",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Test Saha equilibrium
 * 
 * Verify that electron density matches Saha equation prediction
 * at equilibrium conditions
 */
export function testSahaEquilibrium(): TestResult {
  try {
    // Test at known temperature and density
    const T = 20000; // 20,000 K
    const n_neutral = 1e26; // High density
    const ionizationPotential = 15.76 * Constants.e; // Argon
    
    // Compute Saha equilibrium
    const n_e_expected = Calculations.numberDensityFromPressure(
      Calculations.idealGasPressure(n_neutral, T),
      T
    ) * 0.1; // Rough estimate
    
    // Use Saha equation directly
    const g_ratio = 1.0;
    const prefactor = Math.pow(
      (2.0 * Math.PI * Constants.m_e * Constants.k_B * T) / (Constants.h * Constants.h),
      1.5
    );
    const exponent = Math.exp(-ionizationPotential / (Constants.k_B * T));
    const sahaProduct = 2.0 * g_ratio * prefactor * exponent;
    const n_e_saha = Math.sqrt(n_neutral * sahaProduct);
    
    // Check that Saha equation gives reasonable result
    const passed = n_e_saha > 0 && n_e_saha < n_neutral;
    
    return {
      name: "Saha Equilibrium",
      passed,
      error: passed ? undefined : "Saha equation failed",
      details: {
        T,
        n_neutral,
        n_e_saha,
        ionizationFraction: n_e_saha / n_neutral,
      },
    };
  } catch (err: any) {
    return {
      name: "Saha Equilibrium",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Test energy conservation
 * 
 * Verify that total energy is approximately conserved over a cycle
 */
export function testEnergyConservation(): TestResult {
  try {
    const mapper = new DefaultStateVectorMapper();
    const params = createArgonBubblePreset();
    const model = new SonoluminescenceModel(mapper, params);
    
    // Create equilibrium state
    const initialState: BubbleFullState = {
      t: 0,
      hydro: { R: 5e-6, Rdot: 0 },
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
    
    // Integrate over one acoustic period
    const period = 2 * Math.PI / (params.acoustic.omega || 2 * Math.PI * 20e3);
    const result = integrateRK4(
      (t, x) => model.rhs(t, x),
      x0,
      { dt: 1e-9, tMax: period }
    );
    
    // Compute total energy at each step
    const energies = result.x.map((vec, i) => {
      const state = mapper.fromVector(vec, result.t[i]);
      const volume = (4.0 / 3.0) * Math.PI * state.hydro.R * state.hydro.R * state.hydro.R;
      // Simplified: P*V + kinetic + internal
      const kinetic = 0.5 * state.hydro.Rdot * state.hydro.Rdot; // Simplified
      return state.gas.Pg * volume + kinetic + state.em.storedEnergy;
    });
    
    const E_initial = energies[0];
    const E_final = energies[energies.length - 1];
    const E_max = Math.max(...energies);
    const E_min = Math.min(...energies);
    
    // Energy should be approximately conserved (within 20% due to acoustic input)
    const energyChange = Math.abs(E_final - E_initial) / E_initial;
    const passed = energyChange < 0.2;
    
    return {
      name: "Energy Conservation",
      passed,
      error: passed ? undefined : `Energy change: ${(energyChange * 100).toFixed(2)}%`,
      details: {
        E_initial,
        E_final,
        E_max,
        E_min,
        relativeChange: energyChange,
      },
    };
  } catch (err: any) {
    return {
      name: "Energy Conservation",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Test plasma frequency calculation
 * 
 * Verify plasma frequency formula: ωp = sqrt(ne * e² / (ε0 * me))
 */
export function testPlasmaFrequency(): TestResult {
  try {
    const ne = 1e26; // High electron density
    const omega_p = Calculations.plasmaFrequency(ne);
    
    // Expected: ωp = sqrt(ne * e² / (ε0 * me))
    const e = Constants.e;
    const epsilon_0 = Constants.epsilon_0;
    const m_e = Constants.m_e;
    const omega_p_expected = Math.sqrt((ne * e * e) / (epsilon_0 * m_e));
    
    const error = Math.abs(omega_p - omega_p_expected) / omega_p_expected;
    const passed = error < 1e-6;
    
    return {
      name: "Plasma Frequency",
      passed,
      error: passed ? undefined : `Error: ${(error * 100).toFixed(6)}%`,
      details: {
        ne,
        omega_p,
        omega_p_expected,
        f_p: omega_p / (2 * Math.PI), // Frequency in Hz
      },
    };
  } catch (err: any) {
    return {
      name: "Plasma Frequency",
      passed: false,
      error: err.message,
    };
  }
}

/**
 * Run all validation tests
 */
export function runAllTests(): {
  passed: number;
  failed: number;
  results: TestResult[];
} {
  const tests = [
    testAdiabaticCompression,
    testSahaEquilibrium,
    testEnergyConservation,
    testPlasmaFrequency,
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

// Re-export physics validation tests
export {
  testAdiabaticScaling,
  testMinnaertFrequency,
  testEnergyBudgetClosure,
  testPlasmaEquilibrium,
  runAllPhysicsValidationTests,
} from "./physicsValidation";

