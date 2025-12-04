import { StateVectorMapper } from "../core/statevector";
import { HydroParams } from "../physics/hydro";
import { ThermoParams } from "../physics/thermoChem";
import { PlasmaParams } from "../physics/plasma";
import { EmCavityParams } from "../physics/emCavity";
import { AcousticParams } from "../physics/acoustic";
import { ReactionParams } from "../physics/reactions";
export interface SonoluminescenceParams {
    hydro: HydroParams;
    thermo: ThermoParams;
    plasma: PlasmaParams;
    em: EmCavityParams;
    acoustic: AcousticParams;
    reactions: ReactionParams;
}
export declare class SonoluminescenceModel {
    private mapper;
    private params;
    constructor(mapper: StateVectorMapper, params: SonoluminescenceParams);
    /**
     * Compute dX/dt as a dense vector, given X as a vector.
     *
     * This is the right-hand side (RHS) function for the ODE system:
     * dX/dt = rhs(t, X)
     *
     * All physics modules compute their derivatives, which are then
     * mapped back into the state vector using the mapper.
     */
    rhs(t: number, x: Float64Array): Float64Array;
}
//# sourceMappingURL=sonoluminescenceModel.d.ts.map