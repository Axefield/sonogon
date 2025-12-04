// io/export.ts
// Utilities for exporting simulation results

import { IntegrationResult } from "../core/integrator";
import { StateVectorMapper } from "../core/statevector";
import { DimensionId } from "../model/types";
import { BubbleFullState } from "../model/types";
import * as fs from "fs";

export interface ExportOptions {
  filename: string;
  variables?: DimensionId[]; // If not specified, export all
  includeTime?: boolean; // Include time column (default: true)
  delimiter?: string; // CSV delimiter (default: ',')
}

/**
 * Export time series data to CSV
 * 
 * Creates a CSV file with columns: time, variable1, variable2, ...
 */
export function exportToCSV(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  options: ExportOptions
): void {
  const {
    filename,
    variables,
    includeTime = true,
    delimiter = ",",
  } = options;

  const layout = mapper.layout;
  const allVariables = variables || Array.from({ length: layout.size }, (_, i) => {
    // Find DimensionId for index i (reverse lookup)
    for (let dim = 0; dim < 100; dim++) {
      try {
        if (layout.indexOf(dim as DimensionId) === i) {
          return dim as DimensionId;
        }
      } catch {
        // Not found, continue
      }
    }
    return undefined;
  }).filter((d): d is DimensionId => d !== undefined);

  // Create header
  const headers: string[] = [];
  if (includeTime) {
    headers.push("time");
  }
  allVariables.forEach((dim) => {
    headers.push(`dim_${dim}`);
  });

  // Create rows
  const rows: string[] = [headers.join(delimiter)];

  for (let i = 0; i < result.t.length; i++) {
    const row: string[] = [];
    if (includeTime) {
      row.push(result.t[i].toString());
    }
    allVariables.forEach((dim) => {
      const idx = layout.indexOf(dim);
      row.push(result.x[i][idx].toString());
    });
    rows.push(row.join(delimiter));
  }

  // Write to file
  fs.writeFileSync(filename, rows.join("\n"), "utf-8");
}

/**
 * Export complete state history to JSON
 * 
 * Exports full BubbleFullState objects for each time step
 */
export function exportToJSON(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  filename: string
): void {
  const states = result.x.map((vec, i) => mapper.fromVector(vec, result.t[i]));

  const data = {
    metadata: {
      nSteps: result.t.length,
      tMin: Math.min(...result.t),
      tMax: Math.max(...result.t),
      nDimensions: mapper.layout.size,
    },
    timeSeries: result.t,
    states: states,
    stats: result.stats,
  };

  fs.writeFileSync(filename, JSON.stringify(data, null, 2), "utf-8");
}

/**
 * Export selected variables to CSV with readable names
 */
export function exportSelectedVariables(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  options: ExportOptions & {
    variableNames?: Partial<Record<DimensionId, string>>; // Custom names for variables
  }
): void {
  const {
    filename,
    variables = [],
    variableNames = {},
    includeTime = true,
    delimiter = ",",
  } = options;

  if (variables.length === 0) {
    throw new Error("At least one variable must be specified");
  }

  // Create header with readable names
  const headers: string[] = [];
  if (includeTime) {
    headers.push("time");
  }
  variables.forEach((dim) => {
    const name = variableNames[dim] || `dim_${dim}`;
    headers.push(name);
  });

  // Create rows
  const rows: string[] = [headers.join(delimiter)];

  const layout = mapper.layout;
  for (let i = 0; i < result.t.length; i++) {
    const row: string[] = [];
    if (includeTime) {
      row.push(result.t[i].toString());
    }
    variables.forEach((dim) => {
      const idx = layout.indexOf(dim);
      row.push(result.x[i][idx].toString());
    });
    rows.push(row.join(delimiter));
  }

  fs.writeFileSync(filename, rows.join("\n"), "utf-8");
}

/**
 * Export common observables to CSV
 * 
 * Exports: time, R, Rdot, Pg, T, ne, Te, totalPower
 */
export function exportObservables(
  result: IntegrationResult,
  mapper: StateVectorMapper,
  filename: string
): void {
  const { estimateEmission } = require("../analysis/observables");

  const headers = [
    "time",
    "R",
    "Rdot",
    "Pg",
    "T",
    "ne",
    "Te",
    "totalPower",
    "blackbodyPower",
    "bremsstrahlungPower",
    "emDecayPower",
  ];

  const rows: string[] = [headers.join(",")];

  const layout = mapper.layout;
  for (let i = 0; i < result.t.length; i++) {
    const state = mapper.fromVector(result.x[i], result.t[i]);
    const emission = estimateEmission(state);

    const row = [
      result.t[i].toString(),
      state.hydro.R.toString(),
      state.hydro.Rdot.toString(),
      state.gas.Pg.toString(),
      state.gas.T.toString(),
      state.plasma.ne.toString(),
      state.plasma.Te.toString(),
      emission.totalPower.toString(),
      emission.blackbodyPower.toString(),
      emission.bremsstrahlungPower.toString(),
      emission.emDecayPower.toString(),
    ];

    rows.push(row.join(","));
  }

  fs.writeFileSync(filename, rows.join("\n"), "utf-8");
}

