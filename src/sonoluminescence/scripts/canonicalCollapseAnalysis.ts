// scripts/canonicalCollapseAnalysis.ts
// Main script: "This is what Sonogon tells us about a collapsing argon bubble"

import { DefaultStateVectorMapper } from "../core/statevector";
import { SonoluminescenceModel } from "../model/sonoluminescenceModel";
import { createArgonBubblePreset } from "../config/presets";
import { createEquilibriumState } from "../analysis/initialStates";
import { runSimulation } from "../simulation/runner";
import { performCanonicalAnalysis, printCanonicalAnalysis } from "../analysis/canonicalAnalysis";
import { generateCollapseVisualizationData, printVisualizationInstructions } from "./collapseVisualization";

/**
 * Main canonical analysis script
 * 
 * This is the "this is what Sonogon tells us about a collapsing argon bubble" script.
 * 
 * Outputs:
 * - Peak emission (with breakdown by component)
 * - Max electron temperature
 * - Max E_em (stored energy in negative-space state)
 * - Energy budget
 * - Extreme gradients
 * - Mode squeezing metrics
 * 
 * Also generates visualization data for plotting R(t), T(t), ne(t), E_em(t), totalPower(t)
 */
export function runCanonicalCollapseAnalysis(): void {
  console.log("=== SONOGON CANONICAL COLLAPSE ANALYSIS ===\n");
  console.log("Running simulation of collapsing argon bubble...\n");
  
  // Setup
  const mapper = new DefaultStateVectorMapper();
  const params = createArgonBubblePreset();
  
  // Enable advanced features for accurate physics
  params.hydro.useKellerMiksis = true;
  params.thermo.useDetailedHeatTransfer = true;
  params.thermo.useIterativeVanDerWaals = true;
  params.plasma.useDetailedCollisions = true;
  params.em.useDetailedParametricPumping = true;
  params.em.useFrequencyDependentQ = true;
  params.em.Q0 = 1000;
  
  const model = new SonoluminescenceModel(mapper, params);
  const initialState = createEquilibriumState(params, 5e-6); // 5 micron bubble
  
  // Run simulation
  const result = runSimulation({
    model,
    initialState,
    integratorOptions: {
      adaptive: true,
      tolerance: 1e-6,
      dt: 1e-9,
      dtMin: 1e-12,
      dtMax: 1e-8,
      tMax: 1e-5, // 10 μs (multiple acoustic periods)
    },
    logTimeSeries: true, // Enable time series logging
    analysis: {
      computeEmission: true,
      computeGradients: true,
      computeEnergyBudget: true,
    },
  });
  
  console.log(`Simulation complete: ${result.states.length} time steps\n`);
  
  // Perform canonical analysis
  const analysis = performCanonicalAnalysis(result);
  printCanonicalAnalysis(analysis);
  
  // Generate visualization data
  if (result.timeSeriesLog) {
    const vizData = generateCollapseVisualizationData(result);
    
    // Save CSV (in real usage, write to file)
    console.log("Time series CSV data generated (ready for plotting)");
    
    // Print visualization instructions
    printVisualizationInstructions(vizData.csv, vizData.collapseCycle);
  }
}

// Export for use as main script
if (require.main === module) {
  runCanonicalCollapseAnalysis();
}

