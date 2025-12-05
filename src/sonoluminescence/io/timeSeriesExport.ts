// io/timeSeriesExport.ts
// Export time series logs to CSV for visualization

import { TimeSeriesLog } from "../simulation/runner";

/**
 * Export time series log to CSV format
 * 
 * Columns: t, R, Pg, T, ne, Te, E_em, totalPower, Rdot, dPg_dt
 */
export function exportTimeSeriesToCSV(
  timeSeriesLog: TimeSeriesLog[],
  filename?: string
): string {
  const headers = [
    "t [s]",
    "R [m]",
    "Pg [Pa]",
    "T [K]",
    "ne [m^-3]",
    "Te [K]",
    "E_em [J]",
    "totalPower [W]",
    "Rdot [m/s]",
    "dPg_dt [Pa/s]",
  ];
  
  const rows = timeSeriesLog.map((log) => [
    log.t.toExponential(6),
    log.R.toExponential(6),
    log.Pg.toExponential(6),
    log.T.toFixed(2),
    log.ne.toExponential(6),
    log.Te.toFixed(2),
    log.E_em.toExponential(6),
    log.totalPower.toExponential(6),
    log.Rdot.toExponential(6),
    log.dPg_dt?.toExponential(6) || "N/A",
  ]);
  
  const csv = [
    headers.join(","),
    ...rows.map((row) => row.join(",")),
  ].join("\n");
  
  return csv;
}

/**
 * Export time series log to JSON format
 */
export function exportTimeSeriesToJSON(
  timeSeriesLog: TimeSeriesLog[]
): string {
  return JSON.stringify(timeSeriesLog, null, 2);
}

/**
 * Find collapse cycle in time series
 * 
 * Identifies the time window where:
 * 1. R reaches minimum (collapse)
 * 2. E_em rises during extreme gradients
 * 3. E_em decays while totalPower spikes
 */
export function findCollapseCycle(
  timeSeriesLog: TimeSeriesLog[]
): {
  collapseStart: number;
  collapseEnd: number;
  minRIndex: number;
  maxE_emIndex: number;
  maxPowerIndex: number;
} | null {
  if (timeSeriesLog.length < 10) return null;
  
  // Find minimum R (collapse point)
  let minR = Infinity;
  let minRIndex = 0;
  for (let i = 0; i < timeSeriesLog.length; i++) {
    if (timeSeriesLog[i].R < minR) {
      minR = timeSeriesLog[i].R;
      minRIndex = i;
    }
  }
  
  // Find maximum E_em (stored energy peak)
  let maxE_em = -Infinity;
  let maxE_emIndex = 0;
  for (let i = 0; i < timeSeriesLog.length; i++) {
    if (timeSeriesLog[i].E_em > maxE_em) {
      maxE_em = timeSeriesLog[i].E_em;
      maxE_emIndex = i;
    }
  }
  
  // Find maximum totalPower (emission peak)
  let maxPower = -Infinity;
  let maxPowerIndex = 0;
  for (let i = 0; i < timeSeriesLog.length; i++) {
    if (timeSeriesLog[i].totalPower > maxPower) {
      maxPower = timeSeriesLog[i].totalPower;
      maxPowerIndex = i;
    }
  }
  
  // Define collapse window: from before min R to after max power
  const windowBefore = Math.max(0, minRIndex - 50);
  const windowAfter = Math.min(timeSeriesLog.length - 1, maxPowerIndex + 50);
  
  return {
    collapseStart: windowBefore,
    collapseEnd: windowAfter,
    minRIndex,
    maxE_emIndex,
    maxPowerIndex,
  };
}

