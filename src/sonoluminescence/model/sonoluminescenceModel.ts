// model/sonoluminescenceModel.ts
import { BubbleFullState, DimensionId } from "./types";
import { StateVectorMapper } from "../core/statevector";
import { HydroParams, computeHydroDerivatives } from "../physics/hydro";
import {
  ThermoParams,
  computeThermoDerivatives,
  computeSpeciesDerivatives,
  computeInternalEnergyDerivatives,
} from "../physics/thermoChem";
import { PlasmaParams, computePlasmaDerivatives } from "../physics/plasma";
import { EmCavityParams, computeEmDerivatives } from "../physics/emCavity";
import { AcousticParams, computeAcousticState } from "../physics/acoustic";
import {
  ReactionParams,
  computeReactionDerivatives,
  computeSpeciesReactionRates,
} from "../physics/reactions";

export interface SonoluminescenceParams {
  hydro: HydroParams;
  thermo: ThermoParams;
  plasma: PlasmaParams;
  em: EmCavityParams;
  acoustic: AcousticParams;
  reactions: ReactionParams;
}

export class SonoluminescenceModel {
  constructor(
    private mapper: StateVectorMapper,
    private params: SonoluminescenceParams
  ) {}

  /**
   * Compute dX/dt as a dense vector, given X as a vector.
   * 
   * This is the right-hand side (RHS) function for the ODE system:
   * dX/dt = rhs(t, X)
   * 
   * All physics modules compute their derivatives, which are then
   * mapped back into the state vector using the mapper.
   */
  rhs(t: number, x: Float64Array): Float64Array {
    const state: BubbleFullState = this.mapper.fromVector(x, t);
    const dxdt = new Float64Array(x.length);
    const idx = (dim: DimensionId) => this.mapper.layout.indexOf(dim);

    // Acoustic
    const { dPhaseDt, Pacoustic } = computeAcousticState(
      t,
      state.acoustic,
      this.params.acoustic
    );
    dxdt[idx(DimensionId.AcousticPhase)] = dPhaseDt;

    // Hydro (pass gamma from thermo params for Keller-Miksis)
    const gamma = this.params.thermo.gamma;
    const hydroDeriv = computeHydroDerivatives(
      state,
      this.params.hydro,
      Pacoustic,
      gamma
    );
    dxdt[idx(DimensionId.Radius)] = hydroDeriv.dRdt;
    dxdt[idx(DimensionId.RadiusVelocity)] = hydroDeriv.dRdotdt;

    // Thermo
    const thermoDeriv = computeThermoDerivatives(state, this.params.thermo);
    dxdt[idx(DimensionId.GasPressure)] = thermoDeriv.dPgdt;
    dxdt[idx(DimensionId.GasTemperature)] = thermoDeriv.dTdt;

    // Species (compression + reactions + diffusion)
    const speciesDeriv = computeSpeciesDerivatives(state, this.params.thermo);
    const speciesReactionRates = computeSpeciesReactionRates(
      state,
      this.params.reactions
    );
    dxdt[idx(DimensionId.Species_H2O)] =
      speciesDeriv.dNdt.H2O + (speciesReactionRates.H2O || 0);
    dxdt[idx(DimensionId.Species_O2)] =
      speciesDeriv.dNdt.O2 + (speciesReactionRates.O2 || 0);
    dxdt[idx(DimensionId.Species_N2)] =
      speciesDeriv.dNdt.N2 + (speciesReactionRates.N2 || 0);
    dxdt[idx(DimensionId.Species_Ar)] =
      speciesDeriv.dNdt.Ar + (speciesReactionRates.Ar || 0);
    dxdt[idx(DimensionId.Species_Xe)] =
      speciesDeriv.dNdt.Xe + (speciesReactionRates.Xe || 0);
    dxdt[idx(DimensionId.Species_H)] =
      speciesDeriv.dNdt.H + (speciesReactionRates.H || 0);
    dxdt[idx(DimensionId.Species_O)] =
      speciesDeriv.dNdt.O + (speciesReactionRates.O || 0);
    dxdt[idx(DimensionId.Species_OH)] =
      speciesDeriv.dNdt.OH + (speciesReactionRates.OH || 0);
    dxdt[idx(DimensionId.Species_N)] =
      speciesDeriv.dNdt.N + (speciesReactionRates.N || 0);

    // Plasma
    const plasmaDeriv = computePlasmaDerivatives(state, this.params.plasma);
    dxdt[idx(DimensionId.ElectronDensity)] = plasmaDeriv.dn_edt;
    dxdt[idx(DimensionId.ElectronTemperature)] = plasmaDeriv.dTedT;
    dxdt[idx(DimensionId.IonizationFraction)] = plasmaDeriv.dIonFracDt;

    // Internal energy partitions
    const internalEnergyDeriv = computeInternalEnergyDerivatives(
      state,
      this.params.thermo
    );
    dxdt[idx(DimensionId.InternalEnergy_Translational)] =
      internalEnergyDeriv.dE_trans_dt;
    dxdt[idx(DimensionId.InternalEnergy_Rotational)] =
      internalEnergyDeriv.dE_rot_dt;
    dxdt[idx(DimensionId.InternalEnergy_Vibrational)] =
      internalEnergyDeriv.dE_vib_dt;
    dxdt[idx(DimensionId.InternalEnergy_Electronic)] =
      internalEnergyDeriv.dE_elec_dt;

    // EM cavity modes
    const emDeriv = computeEmDerivatives(state, this.params.em);
    if (emDeriv.dModes.length >= 1) {
      dxdt[idx(DimensionId.EmMode0_Re)] = emDeriv.dModes[0].dRe;
      dxdt[idx(DimensionId.EmMode0_Im)] = emDeriv.dModes[0].dIm;
    }
    if (emDeriv.dModes.length >= 2) {
      dxdt[idx(DimensionId.EmMode1_Re)] = emDeriv.dModes[1].dRe;
      dxdt[idx(DimensionId.EmMode1_Im)] = emDeriv.dModes[1].dIm;
    }
    if (emDeriv.dModes.length >= 3) {
      dxdt[idx(DimensionId.EmMode2_Re)] = emDeriv.dModes[2].dRe;
      dxdt[idx(DimensionId.EmMode2_Im)] = emDeriv.dModes[2].dIm;
    }
    dxdt[idx(DimensionId.EmStoredEnergy)] = emDeriv.dStoredEnergyDt;

    // Reactions
    const reactionDeriv = computeReactionDerivatives(
      state,
      this.params.reactions
    );
    if (reactionDeriv.dXidt.length >= 1) {
      dxdt[idx(DimensionId.ReactionProgress_0)] = reactionDeriv.dXidt[0];
    }
    if (reactionDeriv.dXidt.length >= 2) {
      dxdt[idx(DimensionId.ReactionProgress_1)] = reactionDeriv.dXidt[1];
    }
    if (reactionDeriv.dXidt.length >= 3) {
      dxdt[idx(DimensionId.ReactionProgress_2)] = reactionDeriv.dXidt[2];
    }

    return dxdt;
  }
}
