import { DimensionId, BubbleFullState } from "../model/types";
export interface StateVectorLayout {
    size: number;
    indexOf(dim: DimensionId): number;
}
export interface StateVectorMapper {
    layout: StateVectorLayout;
    toVector(state: BubbleFullState): Float64Array;
    fromVector(vec: Float64Array, t: number): BubbleFullState;
}
/**
 * Concrete implementation of StateVectorLayout
 * Maps each DimensionId to a unique index in the state vector
 */
export declare class DefaultStateVectorLayout implements StateVectorLayout {
    private readonly dimensionMap;
    readonly size: number;
    constructor();
    indexOf(dim: DimensionId): number;
}
/**
 * Concrete implementation of StateVectorMapper
 * Converts between BubbleFullState and Float64Array
 */
export declare class DefaultStateVectorMapper implements StateVectorMapper {
    readonly layout: StateVectorLayout;
    private readonly speciesOrder;
    constructor(layout?: StateVectorLayout);
    toVector(state: BubbleFullState): Float64Array;
    fromVector(vec: Float64Array, t: number): BubbleFullState;
}
//# sourceMappingURL=statevector.d.ts.map