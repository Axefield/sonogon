// scripts/collapseVisualization.ts
// Visualization script: "collapse emits/decays inside a negative space" made real

import { SimulationResult } from "../simulation/runner";
import { exportTimeSeriesToCSV, findCollapseCycle } from "../io/timeSeriesExport";

/**
 * Generate visualization data for collapse cycle
 * 
 * Plots R(t), T(t), ne(t), E_em(t), totalPower(t)
 * This plot is the "collapse emits/decays inside a negative space" made real.
 */
export function generateCollapseVisualizationData(
  result: SimulationResult
): {
  csv: string;
  collapseCycle?: {
    start: number;
    end: number;
    minRIndex: number;
    maxE_emIndex: number;
    maxPowerIndex: number;
  };
} {
  if (!result.timeSeriesLog || result.timeSeriesLog.length === 0) {
    throw new Error("Time series log not available. Run simulation with logTimeSeries: true");
  }
  
  // Export full time series to CSV
  const csv = exportTimeSeriesToCSV(result.timeSeriesLog);
  
  // Find collapse cycle
  const collapseCycle = findCollapseCycle(result.timeSeriesLog);
  
  return {
    csv,
    collapseCycle: collapseCycle ? {
      start: collapseCycle.collapseStart,
      end: collapseCycle.collapseEnd,
      minRIndex: collapseCycle.minRIndex,
      maxE_emIndex: collapseCycle.maxE_emIndex,
      maxPowerIndex: collapseCycle.maxPowerIndex,
    } : undefined,
  };
}

/**
 * Print visualization instructions
 */
export function printVisualizationInstructions(
  csv: string,
  collapseCycle?: { start: number; end: number; minRIndex: number; maxE_emIndex: number; maxPowerIndex: number }
): void {
  console.log("\n=== COLLAPSE VISUALIZATION ===");
  console.log("Plot: R(t), T(t), ne(t), E_em(t), totalPower(t)");
  console.log("This plot is the 'collapse emits/decays inside a negative space' made real.\n");
  
  console.log("CSV data exported. Use the following Python code to plot:\n");
  
  console.log(`
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load CSV
df = pd.read_csv('time_series.csv')

# Create figure with subplots
fig, axes = plt.subplots(5, 1, figsize=(12, 10), sharex=True)

# Plot R(t)
axes[0].plot(df['t [s]'], df['R [m]'] * 1e6, 'b-', linewidth=2)
axes[0].set_ylabel('R [μm]', fontsize=12)
axes[0].grid(True, alpha=0.3)
axes[0].set_title('Bubble Collapse: Negative-Space Emission', fontsize=14, fontweight='bold')

# Plot T(t)
axes[1].plot(df['t [s]'], df['T [K]'], 'r-', linewidth=2)
axes[1].set_ylabel('T [K]', fontsize=12)
axes[1].grid(True, alpha=0.3)

# Plot ne(t)
axes[2].plot(df['t [s]'], df['ne [m^-3]'], 'g-', linewidth=2)
axes[2].set_ylabel('ne [m⁻³]', fontsize=12)
axes[2].set_yscale('log')
axes[2].grid(True, alpha=0.3)

# Plot E_em(t) - THE NEGATIVE-SPACE STATE
axes[3].plot(df['t [s]'], df['E_em [J]'], 'm-', linewidth=2, label='E_em (stored energy)')
axes[3].set_ylabel('E_em [J]', fontsize=12)
axes[3].set_yscale('log')
axes[3].grid(True, alpha=0.3)
axes[3].legend()

# Plot totalPower(t) - THE PHOTON EMISSION
axes[4].plot(df['t [s]'], df['totalPower [W]'], 'orange', linewidth=2, label='totalPower (photon emission)')
axes[4].set_ylabel('totalPower [W]', fontsize=12)
axes[4].set_xlabel('Time [s]', fontsize=12)
axes[4].set_yscale('log')
axes[4].grid(True, alpha=0.3)
axes[4].legend()

plt.tight_layout()
plt.savefig('collapse_negative_space.png', dpi=300, bbox_inches='tight')
print("Plot saved as 'collapse_negative_space.png'")

# Highlight collapse cycle if available
${collapseCycle ? `
# Mark collapse cycle
collapse_start_idx = ${collapseCycle.start}
collapse_end_idx = ${collapseCycle.end}
minR_idx = ${collapseCycle.minRIndex}
maxE_em_idx = ${collapseCycle.maxE_emIndex}
maxPower_idx = ${collapseCycle.maxPowerIndex}

# Shade collapse region
for ax in axes:
    ax.axvspan(df['t [s]'].iloc[collapse_start_idx], 
               df['t [s]'].iloc[collapse_end_idx], 
               alpha=0.2, color='yellow', label='Collapse cycle')
    ax.axvline(df['t [s]'].iloc[minR_idx], color='blue', linestyle='--', alpha=0.5, label='Min R')
    ax.axvline(df['t [s]'].iloc[maxE_em_idx], color='magenta', linestyle='--', alpha=0.5, label='Max E_em')
    ax.axvline(df['t [s]'].iloc[maxPower_idx], color='orange', linestyle='--', alpha=0.5, label='Max Power')
    if ax == axes[0]:
        ax.legend(loc='upper right', fontsize=8)

plt.savefig('collapse_negative_space_annotated.png', dpi=300, bbox_inches='tight')
print("Annotated plot saved as 'collapse_negative_space_annotated.png'")
` : ''}
`);
  
  if (collapseCycle) {
    console.log("\nCollapse cycle detected:");
    console.log(`  Min R at index: ${collapseCycle.minRIndex}`);
    console.log(`  Max E_em at index: ${collapseCycle.maxE_emIndex}`);
    console.log(`  Max Power at index: ${collapseCycle.maxPowerIndex}`);
    console.log(`  Cycle window: indices ${collapseCycle.start} to ${collapseCycle.end}`);
  }
  
  console.log("\n=== END VISUALIZATION INSTRUCTIONS ===\n");
}

