import { IntegrationResult } from "../core/integrator";
import { StateVectorMapper } from "../core/statevector";
import { DimensionId } from "../model/types";
export interface ExportOptions {
    filename: string;
    variables?: DimensionId[];
    includeTime?: boolean;
    delimiter?: string;
}
/**
 * Export time series data to CSV
 *
 * Creates a CSV file with columns: time, variable1, variable2, ...
 */
export declare function exportToCSV(result: IntegrationResult, mapper: StateVectorMapper, options: ExportOptions): void;
/**
 * Export complete state history to JSON
 *
 * Exports full BubbleFullState objects for each time step
 */
export declare function exportToJSON(result: IntegrationResult, mapper: StateVectorMapper, filename: string): void;
/**
 * Export selected variables to CSV with readable names
 */
export declare function exportSelectedVariables(result: IntegrationResult, mapper: StateVectorMapper, options: ExportOptions & {
    variableNames?: Partial<Record<DimensionId, string>>;
}): void;
/**
 * Export common observables to CSV
 *
 * Exports: time, R, Rdot, Pg, T, ne, Te, totalPower
 */
export declare function exportObservables(result: IntegrationResult, mapper: StateVectorMapper, filename: string): void;
//# sourceMappingURL=export.d.ts.map