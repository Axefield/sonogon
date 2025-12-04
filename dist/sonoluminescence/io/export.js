"use strict";
// io/export.ts
// Utilities for exporting simulation results
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || (function () {
    var ownKeys = function(o) {
        ownKeys = Object.getOwnPropertyNames || function (o) {
            var ar = [];
            for (var k in o) if (Object.prototype.hasOwnProperty.call(o, k)) ar[ar.length] = k;
            return ar;
        };
        return ownKeys(o);
    };
    return function (mod) {
        if (mod && mod.__esModule) return mod;
        var result = {};
        if (mod != null) for (var k = ownKeys(mod), i = 0; i < k.length; i++) if (k[i] !== "default") __createBinding(result, mod, k[i]);
        __setModuleDefault(result, mod);
        return result;
    };
})();
Object.defineProperty(exports, "__esModule", { value: true });
exports.exportToCSV = exportToCSV;
exports.exportToJSON = exportToJSON;
exports.exportSelectedVariables = exportSelectedVariables;
exports.exportObservables = exportObservables;
const fs = __importStar(require("fs"));
/**
 * Export time series data to CSV
 *
 * Creates a CSV file with columns: time, variable1, variable2, ...
 */
function exportToCSV(result, mapper, options) {
    const { filename, variables, includeTime = true, delimiter = ",", } = options;
    const layout = mapper.layout;
    const allVariables = variables || Array.from({ length: layout.size }, (_, i) => {
        // Find DimensionId for index i (reverse lookup)
        for (let dim = 0; dim < 100; dim++) {
            try {
                if (layout.indexOf(dim) === i) {
                    return dim;
                }
            }
            catch {
                // Not found, continue
            }
        }
        return undefined;
    }).filter((d) => d !== undefined);
    // Create header
    const headers = [];
    if (includeTime) {
        headers.push("time");
    }
    allVariables.forEach((dim) => {
        headers.push(`dim_${dim}`);
    });
    // Create rows
    const rows = [headers.join(delimiter)];
    for (let i = 0; i < result.t.length; i++) {
        const row = [];
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
function exportToJSON(result, mapper, filename) {
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
function exportSelectedVariables(result, mapper, options) {
    const { filename, variables = [], variableNames = {}, includeTime = true, delimiter = ",", } = options;
    if (variables.length === 0) {
        throw new Error("At least one variable must be specified");
    }
    // Create header with readable names
    const headers = [];
    if (includeTime) {
        headers.push("time");
    }
    variables.forEach((dim) => {
        const name = variableNames[dim] || `dim_${dim}`;
        headers.push(name);
    });
    // Create rows
    const rows = [headers.join(delimiter)];
    const layout = mapper.layout;
    for (let i = 0; i < result.t.length; i++) {
        const row = [];
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
function exportObservables(result, mapper, filename) {
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
    const rows = [headers.join(",")];
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
//# sourceMappingURL=export.js.map