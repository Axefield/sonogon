import { BubbleFullState, SpeciesId } from "../model/types";
export interface ThermoParams {
    gamma: number;
    R0: number;
    Pg0: number;
    T0: number;
    heatLossCoeff: number;
    useNonIdealGas?: boolean;
    vanDerWaals_a?: number;
    vanDerWaals_b?: number;
    useTemperatureDependentGamma?: boolean;
    gammaFunction?: (T: number) => number;
    tau_trans?: number;
    tau_rot?: number;
    tau_vib?: number;
    tau_elec?: number;
    useDetailedEnergyExchange?: boolean;
    useLandauTellerRelaxation?: boolean;
    useDetailedHeatTransfer?: boolean;
    thermalConductivity?: number;
    liquidThermalConductivity?: number;
    convectiveCoeff?: number;
    includeRadiation?: boolean;
    useDetailedHeatCapacity?: boolean;
    useSpeciesDiffusion?: boolean;
    diffusionCoefficients?: Record<string, number>;
    liquidConcentrations?: Record<string, number>;
    useIterativeVanDerWaals?: boolean;
    vdwIterationTolerance?: number;
    vdwMaxIterations?: number;
}
export interface ThermoDerivatives {
    dPgdt: number;
    dTdt: number;
}
export interface SpeciesDerivatives {
    dNdt: Record<SpeciesId, number>;
}
export interface InternalEnergyDerivatives {
    dE_trans_dt: number;
    dE_rot_dt: number;
    dE_vib_dt: number;
    dE_elec_dt: number;
}
/**
 * Compute thermodynamic derivatives (pressure and temperature)
 *
 * For adiabatic compression/expansion:
 * - dPg/dt from volume change: dPg/dt = -gamma * Pg * (3*Rdot/R)
 * - dT/dt from adiabatic relation: dT/dt = (gamma-1) * T * (3*Rdot/R) + heat_loss
 *
 * Enhanced with:
 * - Non-ideal gas effects (van der Waals)
 * - Temperature-dependent gamma
 */
export declare function computeThermoDerivatives(state: BubbleFullState, params: ThermoParams): ThermoDerivatives;
/**
 * Compute species number density derivatives
 *
 * Species evolve due to:
 * 1. Compression/expansion (conservation of mass)
 * 2. Chemical reactions (dissociation, recombination)
 * 3. Diffusion through bubble wall (if enabled)
 */
export declare function computeSpeciesDerivatives(state: BubbleFullState, params?: ThermoParams): SpeciesDerivatives;
/**
 * Compute internal energy partition derivatives
 *
 * Tracks energy in different modes:
 * - Translational: (3/2) * n * k_B * T (for monatomic)
 * - Rotational: n * k_B * T (for diatomic, per rotational mode)
 * - Vibrational: quantum harmonic oscillator energy
 * - Electronic: excited state populations
 *
 * Energy exchange between modes via relaxation processes
 * Enhanced with temperature-dependent relaxation times
 */
export declare function computeInternalEnergyDerivatives(state: BubbleFullState, params: ThermoParams): InternalEnergyDerivatives;
//# sourceMappingURL=thermoChem.d.ts.map