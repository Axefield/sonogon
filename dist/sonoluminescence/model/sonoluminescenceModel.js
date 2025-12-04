"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.SonoluminescenceModel = void 0;
// model/sonoluminescenceModel.ts
const types_1 = require("./types");
const hydro_1 = require("../physics/hydro");
const thermoChem_1 = require("../physics/thermoChem");
const plasma_1 = require("../physics/plasma");
const emCavity_1 = require("../physics/emCavity");
const acoustic_1 = require("../physics/acoustic");
const reactions_1 = require("../physics/reactions");
class SonoluminescenceModel {
    constructor(mapper, params) {
        this.mapper = mapper;
        this.params = params;
    }
    /**
     * Compute dX/dt as a dense vector, given X as a vector.
     *
     * This is the right-hand side (RHS) function for the ODE system:
     * dX/dt = rhs(t, X)
     *
     * All physics modules compute their derivatives, which are then
     * mapped back into the state vector using the mapper.
     */
    rhs(t, x) {
        const state = this.mapper.fromVector(x, t);
        const dxdt = new Float64Array(x.length);
        const idx = (dim) => this.mapper.layout.indexOf(dim);
        // Acoustic
        const { dPhaseDt, Pacoustic } = (0, acoustic_1.computeAcousticState)(t, state.acoustic, this.params.acoustic);
        dxdt[idx(types_1.DimensionId.AcousticPhase)] = dPhaseDt;
        // Hydro (pass gamma from thermo params for Keller-Miksis)
        const gamma = this.params.thermo.gamma;
        const hydroDeriv = (0, hydro_1.computeHydroDerivatives)(state, this.params.hydro, Pacoustic, gamma);
        dxdt[idx(types_1.DimensionId.Radius)] = hydroDeriv.dRdt;
        dxdt[idx(types_1.DimensionId.RadiusVelocity)] = hydroDeriv.dRdotdt;
        // Thermo
        const thermoDeriv = (0, thermoChem_1.computeThermoDerivatives)(state, this.params.thermo);
        dxdt[idx(types_1.DimensionId.GasPressure)] = thermoDeriv.dPgdt;
        dxdt[idx(types_1.DimensionId.GasTemperature)] = thermoDeriv.dTdt;
        // Species (compression + reactions + diffusion)
        const speciesDeriv = (0, thermoChem_1.computeSpeciesDerivatives)(state, this.params.thermo);
        const speciesReactionRates = (0, reactions_1.computeSpeciesReactionRates)(state, this.params.reactions);
        dxdt[idx(types_1.DimensionId.Species_H2O)] =
            speciesDeriv.dNdt.H2O + (speciesReactionRates.H2O || 0);
        dxdt[idx(types_1.DimensionId.Species_O2)] =
            speciesDeriv.dNdt.O2 + (speciesReactionRates.O2 || 0);
        dxdt[idx(types_1.DimensionId.Species_N2)] =
            speciesDeriv.dNdt.N2 + (speciesReactionRates.N2 || 0);
        dxdt[idx(types_1.DimensionId.Species_Ar)] =
            speciesDeriv.dNdt.Ar + (speciesReactionRates.Ar || 0);
        dxdt[idx(types_1.DimensionId.Species_Xe)] =
            speciesDeriv.dNdt.Xe + (speciesReactionRates.Xe || 0);
        dxdt[idx(types_1.DimensionId.Species_H)] =
            speciesDeriv.dNdt.H + (speciesReactionRates.H || 0);
        dxdt[idx(types_1.DimensionId.Species_O)] =
            speciesDeriv.dNdt.O + (speciesReactionRates.O || 0);
        dxdt[idx(types_1.DimensionId.Species_OH)] =
            speciesDeriv.dNdt.OH + (speciesReactionRates.OH || 0);
        dxdt[idx(types_1.DimensionId.Species_N)] =
            speciesDeriv.dNdt.N + (speciesReactionRates.N || 0);
        // Plasma
        const plasmaDeriv = (0, plasma_1.computePlasmaDerivatives)(state, this.params.plasma);
        dxdt[idx(types_1.DimensionId.ElectronDensity)] = plasmaDeriv.dn_edt;
        dxdt[idx(types_1.DimensionId.ElectronTemperature)] = plasmaDeriv.dTedT;
        dxdt[idx(types_1.DimensionId.IonizationFraction)] = plasmaDeriv.dIonFracDt;
        // Internal energy partitions
        const internalEnergyDeriv = (0, thermoChem_1.computeInternalEnergyDerivatives)(state, this.params.thermo);
        dxdt[idx(types_1.DimensionId.InternalEnergy_Translational)] =
            internalEnergyDeriv.dE_trans_dt;
        dxdt[idx(types_1.DimensionId.InternalEnergy_Rotational)] =
            internalEnergyDeriv.dE_rot_dt;
        dxdt[idx(types_1.DimensionId.InternalEnergy_Vibrational)] =
            internalEnergyDeriv.dE_vib_dt;
        dxdt[idx(types_1.DimensionId.InternalEnergy_Electronic)] =
            internalEnergyDeriv.dE_elec_dt;
        // EM cavity modes
        const emDeriv = (0, emCavity_1.computeEmDerivatives)(state, this.params.em);
        if (emDeriv.dModes.length >= 1) {
            dxdt[idx(types_1.DimensionId.EmMode0_Re)] = emDeriv.dModes[0].dRe;
            dxdt[idx(types_1.DimensionId.EmMode0_Im)] = emDeriv.dModes[0].dIm;
        }
        if (emDeriv.dModes.length >= 2) {
            dxdt[idx(types_1.DimensionId.EmMode1_Re)] = emDeriv.dModes[1].dRe;
            dxdt[idx(types_1.DimensionId.EmMode1_Im)] = emDeriv.dModes[1].dIm;
        }
        if (emDeriv.dModes.length >= 3) {
            dxdt[idx(types_1.DimensionId.EmMode2_Re)] = emDeriv.dModes[2].dRe;
            dxdt[idx(types_1.DimensionId.EmMode2_Im)] = emDeriv.dModes[2].dIm;
        }
        dxdt[idx(types_1.DimensionId.EmStoredEnergy)] = emDeriv.dStoredEnergyDt;
        // Reactions
        const reactionDeriv = (0, reactions_1.computeReactionDerivatives)(state, this.params.reactions);
        if (reactionDeriv.dXidt.length >= 1) {
            dxdt[idx(types_1.DimensionId.ReactionProgress_0)] = reactionDeriv.dXidt[0];
        }
        if (reactionDeriv.dXidt.length >= 2) {
            dxdt[idx(types_1.DimensionId.ReactionProgress_1)] = reactionDeriv.dXidt[1];
        }
        if (reactionDeriv.dXidt.length >= 3) {
            dxdt[idx(types_1.DimensionId.ReactionProgress_2)] = reactionDeriv.dXidt[2];
        }
        return dxdt;
    }
}
exports.SonoluminescenceModel = SonoluminescenceModel;
//# sourceMappingURL=sonoluminescenceModel.js.map