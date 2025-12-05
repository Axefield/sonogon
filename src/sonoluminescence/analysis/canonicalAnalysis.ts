// analysis/canonicalAnalysis.ts
// Canonical "this is what Sonogon tells us about a collapsing argon bubble" analysis

import { BubbleFullState } from "../model/types";
import { SimulationResult } from "../simulation/runner";
import { estimateEmission } from "./observables";
import { computeEnergyBudget, detectExtremeGradients, computeModeSqueezingMetrics, EnergyBudget } from "./diagnostics";

export interface CanonicalAnalysisResult {
  // Peak values
  peakEmission: {
    power: number; // [W]
    time: number; // [s]
    index: number;
    breakdown: {
      blackbody: number;
      bremsstrahlung: number;
      emDecay: number;
    };
  };
  
  maxElectronTemperature: {
    Te: number; // [K]
    time: number; // [s]
    index: number;
  };
  
  maxE_em: {
    E_em: number; // [J]
    time: number; // [s]
    index: number;
  };
  
  // Energy budget (averaged over cycle or at peak)
  energyBudget: {
    acousticWorkIn: number;
    compressionWork: number;
    heatLoss: number;
    lightEmission: number;
    viscousDissipation: number;
    chemicalEnergy: number;
    emStoredEnergy: number;
    totalEnergy: number;
  };
  
  // Extreme gradients
  extremeGradients: {
    detected: boolean;
    maxGradient: number;
    time: number;
    index: number;
    metrics: {
      dR_dt: number;
      dPg_dt: number;
      dT_dt: number;
    };
  };
  
  // Mode squeezing metrics (at peak E_em)
  modeSqueezing: {
    pumpedEnergy: number;
    releasedEnergy: number;
    netEnergyChange: number;
    squeezingDetected: boolean;
    pumpEfficiency: number;
  };
}

/**
 * Perform canonical analysis on a simulation result
 * 
 * This is the "this is what Sonogon tells us about a collapsing argon bubble" script.
 * Computes:
 * - Peak emission (and breakdown by component)
 * - Max electron temperature
 * - Max E_em (stored energy in negative-space state)
 * - Energy budget
 * - Extreme gradients
 * - Mode squeezing metrics
 */
export function performCanonicalAnalysis(
  result: SimulationResult
): CanonicalAnalysisResult {
  const { states, timeSeries, analysis } = result;
  
  if (states.length === 0) {
    throw new Error("No states in simulation result");
  }
  
  // Find peak emission
  let peakEmissionPower = -Infinity;
  let peakEmissionIndex = 0;
  let peakEmissionBreakdown = { blackbody: 0, bremsstrahlung: 0, emDecay: 0 };
  
  for (let i = 0; i < states.length; i++) {
    const emission = estimateEmission(states[i], false);
    if (emission.totalPower > peakEmissionPower) {
      peakEmissionPower = emission.totalPower;
      peakEmissionIndex = i;
      peakEmissionBreakdown = {
        blackbody: emission.blackbodyPower,
        bremsstrahlung: emission.bremsstrahlungPower,
        emDecay: emission.emDecayPower,
      };
    }
  }
  
  // Find max electron temperature
  let maxTe = -Infinity;
  let maxTeIndex = 0;
  for (let i = 0; i < states.length; i++) {
    if (states[i].plasma.Te > maxTe) {
      maxTe = states[i].plasma.Te;
      maxTeIndex = i;
    }
  }
  
  // Find max E_em
  let maxE_em = -Infinity;
  let maxE_emIndex = 0;
  for (let i = 0; i < states.length; i++) {
    if (states[i].em.storedEnergy > maxE_em) {
      maxE_em = states[i].em.storedEnergy;
      maxE_emIndex = i;
    }
  }
  
  // Compute energy budget (at peak emission or averaged)
  let energyBudget: EnergyBudget;
  if (analysis?.energyBudgets && analysis.energyBudgets.length > 0) {
    // Use provided energy budgets (averaged)
    const budgets = analysis.energyBudgets;
    energyBudget = {
      acousticWorkIn: budgets.reduce((sum: number, b: EnergyBudget) => sum + b.acousticWorkIn, 0) / budgets.length,
      compressionWork: budgets.reduce((sum: number, b: EnergyBudget) => sum + Math.abs(b.compressionWork), 0) / budgets.length,
      heatLoss: budgets.reduce((sum: number, b: EnergyBudget) => sum + b.heatLoss, 0) / budgets.length,
      lightEmission: budgets.reduce((sum: number, b: EnergyBudget) => sum + b.lightEmission, 0) / budgets.length,
      viscousDissipation: budgets.reduce((sum: number, b: EnergyBudget) => sum + b.viscousDissipation, 0) / budgets.length,
      chemicalEnergy: budgets.reduce((sum: number, b: EnergyBudget) => sum + Math.abs(b.chemicalEnergy), 0) / budgets.length,
      emStoredEnergy: budgets.reduce((sum: number, b: EnergyBudget) => sum + Math.abs(b.emStoredEnergy), 0) / budgets.length,
      totalEnergy: budgets.reduce((sum: number, b: EnergyBudget) => sum + Math.abs(b.totalEnergy), 0) / budgets.length,
    };
  } else {
    // Compute energy budget at peak emission
    if (peakEmissionIndex > 0) {
      const dt = timeSeries.t[peakEmissionIndex] - timeSeries.t[peakEmissionIndex - 1];
      const acousticPower = 1e3; // Estimate: 1 kW (typical)
      energyBudget = computeEnergyBudget(
        states[peakEmissionIndex - 1],
        states[peakEmissionIndex],
        dt,
        acousticPower
      );
    } else {
      // Fallback: use first state
      energyBudget = {
        acousticWorkIn: 0,
        compressionWork: 0,
        heatLoss: 0,
        lightEmission: 0,
        viscousDissipation: 0,
        chemicalEnergy: 0,
        emStoredEnergy: 0,
        totalEnergy: 0,
      };
    }
  }
  
  // Detect extreme gradients
  let extremeGradients;
  if (peakEmissionIndex > 0) {
    const dt = timeSeries.t[peakEmissionIndex] - timeSeries.t[peakEmissionIndex - 1];
    const gradientResult = detectExtremeGradients(
      states[peakEmissionIndex - 1],
      states[peakEmissionIndex],
      dt
    );
    extremeGradients = {
      detected: gradientResult.detected,
      maxGradient: gradientResult.metrics.maxGradient,
      time: timeSeries.t[peakEmissionIndex],
      index: peakEmissionIndex,
      metrics: {
        dR_dt: gradientResult.metrics.dRdt,
        dPg_dt: gradientResult.metrics.dPgdt,
        dT_dt: gradientResult.metrics.dTdt,
      },
    };
  } else {
    extremeGradients = {
      detected: false,
      maxGradient: 0,
      time: timeSeries.t[0],
      index: 0,
      metrics: { dR_dt: 0, dPg_dt: 0, dT_dt: 0 },
    };
  }
  
  // Compute mode squeezing metrics at max E_em
  let modeSqueezing;
  if (maxE_emIndex > 0) {
    const dt = timeSeries.t[maxE_emIndex] - timeSeries.t[maxE_emIndex - 1];
    modeSqueezing = computeModeSqueezingMetrics(
      states[maxE_emIndex],
      states[maxE_emIndex - 1],
      dt
    );
  } else {
    modeSqueezing = computeModeSqueezingMetrics(states[maxE_emIndex]);
  }
  
  return {
    peakEmission: {
      power: peakEmissionPower,
      time: timeSeries.t[peakEmissionIndex],
      index: peakEmissionIndex,
      breakdown: peakEmissionBreakdown,
    },
    maxElectronTemperature: {
      Te: maxTe,
      time: timeSeries.t[maxTeIndex],
      index: maxTeIndex,
    },
    maxE_em: {
      E_em: maxE_em,
      time: timeSeries.t[maxE_emIndex],
      index: maxE_emIndex,
    },
    energyBudget,
    extremeGradients,
    modeSqueezing: {
      pumpedEnergy: modeSqueezing.pumpedEnergy,
      releasedEnergy: modeSqueezing.releasedEnergy,
      netEnergyChange: modeSqueezing.netEnergyChange,
      squeezingDetected: modeSqueezing.squeezingDetected,
      pumpEfficiency: modeSqueezing.pumpEfficiency,
    },
  };
}

/**
 * Print canonical analysis results in human-readable format
 */
export function printCanonicalAnalysis(analysis: CanonicalAnalysisResult): void {
  console.log("\n=== SONOGON CANONICAL ANALYSIS ===");
  console.log("This is what Sonogon tells us about a collapsing argon bubble:\n");
  
  console.log("PEAK EMISSION:");
  console.log(`  Total Power: ${analysis.peakEmission.power.toExponential(2)} W`);
  console.log(`  Time: ${analysis.peakEmission.time.toExponential(2)} s`);
  console.log(`  Breakdown:`);
  console.log(`    Blackbody (T term): ${analysis.peakEmission.breakdown.blackbody.toExponential(2)} W`);
  console.log(`    Bremsstrahlung (ne, Te): ${analysis.peakEmission.breakdown.bremsstrahlung.toExponential(2)} W`);
  console.log(`    EM Decay (E_em): ${analysis.peakEmission.breakdown.emDecay.toExponential(2)} W`);
  
  console.log("\nMAX ELECTRON TEMPERATURE:");
  console.log(`  Te: ${analysis.maxElectronTemperature.Te.toFixed(0)} K`);
  console.log(`  Time: ${analysis.maxElectronTemperature.time.toExponential(2)} s`);
  
  console.log("\nMAX E_em (STORED ENERGY IN NEGATIVE-SPACE STATE):");
  console.log(`  E_em: ${analysis.maxE_em.E_em.toExponential(2)} J`);
  console.log(`  Time: ${analysis.maxE_em.time.toExponential(2)} s`);
  
  console.log("\nENERGY BUDGET:");
  console.log(`  Acoustic Work In: ${analysis.energyBudget.acousticWorkIn.toExponential(2)} J`);
  console.log(`  Compression Work: ${analysis.energyBudget.compressionWork.toExponential(2)} J`);
  console.log(`  Heat Loss: ${analysis.energyBudget.heatLoss.toExponential(2)} J`);
  console.log(`  Light Emission: ${analysis.energyBudget.lightEmission.toExponential(2)} J`);
  console.log(`  Viscous Dissipation: ${analysis.energyBudget.viscousDissipation.toExponential(2)} J`);
  console.log(`  Chemical Energy: ${analysis.energyBudget.chemicalEnergy.toExponential(2)} J`);
  console.log(`  EM Stored Energy: ${analysis.energyBudget.emStoredEnergy.toExponential(2)} J`);
  console.log(`  Total Energy (balance): ${analysis.energyBudget.totalEnergy.toExponential(2)} J`);
  
  console.log("\nEXTREME GRADIENTS:");
  console.log(`  Detected: ${analysis.extremeGradients.detected}`);
  console.log(`  Max Gradient: ${analysis.extremeGradients.maxGradient.toExponential(2)}`);
  console.log(`  Time: ${analysis.extremeGradients.time.toExponential(2)} s`);
  console.log(`  Metrics:`);
  console.log(`    dR/dt: ${analysis.extremeGradients.metrics.dR_dt.toExponential(2)} m/s`);
  console.log(`    dPg/dt: ${analysis.extremeGradients.metrics.dPg_dt.toExponential(2)} Pa/s`);
  console.log(`    dT/dt: ${analysis.extremeGradients.metrics.dT_dt.toExponential(2)} K/s`);
  
  console.log("\nMODE SQUEEZING METRICS (NEGATIVE-SPACE BEHAVIOR):");
  console.log(`  Pumped Energy: ${analysis.modeSqueezing.pumpedEnergy.toExponential(2)} J`);
  console.log(`  Released Energy: ${analysis.modeSqueezing.releasedEnergy.toExponential(2)} J`);
  console.log(`  Net Energy Change: ${analysis.modeSqueezing.netEnergyChange.toExponential(2)} J`);
  console.log(`  Squeezing Detected: ${analysis.modeSqueezing.squeezingDetected}`);
  console.log(`  Pump Efficiency: ${(analysis.modeSqueezing.pumpEfficiency * 100).toFixed(2)}%`);
  
  console.log("\n=== END ANALYSIS ===\n");
}

