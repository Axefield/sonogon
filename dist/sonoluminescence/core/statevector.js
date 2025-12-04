"use strict";
// statevector.ts
Object.defineProperty(exports, "__esModule", { value: true });
exports.DefaultStateVectorMapper = exports.DefaultStateVectorLayout = void 0;
const types_1 = require("../model/types");
/**
 * Concrete implementation of StateVectorLayout
 * Maps each DimensionId to a unique index in the state vector
 */
class DefaultStateVectorLayout {
    constructor() {
        this.dimensionMap = new Map();
        let index = 0;
        // Hydrodynamic / geometric (2)
        this.dimensionMap.set(types_1.DimensionId.Radius, index++);
        this.dimensionMap.set(types_1.DimensionId.RadiusVelocity, index++);
        // Shape oscillations (4) - optional
        this.dimensionMap.set(types_1.DimensionId.ShapeMode2_Amplitude, index++);
        this.dimensionMap.set(types_1.DimensionId.ShapeMode2_Velocity, index++);
        this.dimensionMap.set(types_1.DimensionId.ShapeMode4_Amplitude, index++);
        this.dimensionMap.set(types_1.DimensionId.ShapeMode4_Velocity, index++);
        // Bubble translation (6) - optional
        this.dimensionMap.set(types_1.DimensionId.BubblePosition_X, index++);
        this.dimensionMap.set(types_1.DimensionId.BubblePosition_Y, index++);
        this.dimensionMap.set(types_1.DimensionId.BubblePosition_Z, index++);
        this.dimensionMap.set(types_1.DimensionId.BubbleVelocity_X, index++);
        this.dimensionMap.set(types_1.DimensionId.BubbleVelocity_Y, index++);
        this.dimensionMap.set(types_1.DimensionId.BubbleVelocity_Z, index++);
        // Gas macro state (2)
        this.dimensionMap.set(types_1.DimensionId.GasPressure, index++);
        this.dimensionMap.set(types_1.DimensionId.GasTemperature, index++);
        // Species number densities (9)
        this.dimensionMap.set(types_1.DimensionId.Species_H2O, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_O2, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_N2, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_Ar, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_Xe, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_H, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_O, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_OH, index++);
        this.dimensionMap.set(types_1.DimensionId.Species_N, index++);
        // Plasma state (3)
        this.dimensionMap.set(types_1.DimensionId.ElectronDensity, index++);
        this.dimensionMap.set(types_1.DimensionId.ElectronTemperature, index++);
        this.dimensionMap.set(types_1.DimensionId.IonizationFraction, index++);
        // Internal energy partitions (4)
        this.dimensionMap.set(types_1.DimensionId.InternalEnergy_Translational, index++);
        this.dimensionMap.set(types_1.DimensionId.InternalEnergy_Rotational, index++);
        this.dimensionMap.set(types_1.DimensionId.InternalEnergy_Vibrational, index++);
        this.dimensionMap.set(types_1.DimensionId.InternalEnergy_Electronic, index++);
        // EM cavity modes (6 for 3 modes)
        this.dimensionMap.set(types_1.DimensionId.EmMode0_Re, index++);
        this.dimensionMap.set(types_1.DimensionId.EmMode0_Im, index++);
        this.dimensionMap.set(types_1.DimensionId.EmMode1_Re, index++);
        this.dimensionMap.set(types_1.DimensionId.EmMode1_Im, index++);
        this.dimensionMap.set(types_1.DimensionId.EmMode2_Re, index++);
        this.dimensionMap.set(types_1.DimensionId.EmMode2_Im, index++);
        // EM stored energy (1)
        this.dimensionMap.set(types_1.DimensionId.EmStoredEnergy, index++);
        // Acoustic phase (1)
        this.dimensionMap.set(types_1.DimensionId.AcousticPhase, index++);
        // Reaction progress variables (3)
        this.dimensionMap.set(types_1.DimensionId.ReactionProgress_0, index++);
        this.dimensionMap.set(types_1.DimensionId.ReactionProgress_1, index++);
        this.dimensionMap.set(types_1.DimensionId.ReactionProgress_2, index++);
        this.size = index;
    }
    indexOf(dim) {
        const idx = this.dimensionMap.get(dim);
        if (idx === undefined) {
            throw new Error(`DimensionId ${dim} not found in layout`);
        }
        return idx;
    }
}
exports.DefaultStateVectorLayout = DefaultStateVectorLayout;
/**
 * Concrete implementation of StateVectorMapper
 * Converts between BubbleFullState and Float64Array
 */
class DefaultStateVectorMapper {
    constructor(layout) {
        this.speciesOrder = ["H2O", "O2", "N2", "Ar", "Xe", "H", "O"];
        this.layout = layout || new DefaultStateVectorLayout();
    }
    toVector(state) {
        const vec = new Float64Array(this.layout.size);
        const idx = (dim) => this.layout.indexOf(dim);
        // Hydrodynamic
        vec[idx(types_1.DimensionId.Radius)] = state.hydro.R;
        vec[idx(types_1.DimensionId.RadiusVelocity)] = state.hydro.Rdot;
        // Shape oscillations
        if (state.shape) {
            vec[idx(types_1.DimensionId.ShapeMode2_Amplitude)] = state.shape.a2;
            vec[idx(types_1.DimensionId.ShapeMode2_Velocity)] = state.shape.a2_dot;
            vec[idx(types_1.DimensionId.ShapeMode4_Amplitude)] = state.shape.a4;
            vec[idx(types_1.DimensionId.ShapeMode4_Velocity)] = state.shape.a4_dot;
        }
        else {
            // Default to zero if not provided
            vec[idx(types_1.DimensionId.ShapeMode2_Amplitude)] = 0;
            vec[idx(types_1.DimensionId.ShapeMode2_Velocity)] = 0;
            vec[idx(types_1.DimensionId.ShapeMode4_Amplitude)] = 0;
            vec[idx(types_1.DimensionId.ShapeMode4_Velocity)] = 0;
        }
        // Bubble translation
        if (state.translation) {
            vec[idx(types_1.DimensionId.BubblePosition_X)] = state.translation.x;
            vec[idx(types_1.DimensionId.BubblePosition_Y)] = state.translation.y;
            vec[idx(types_1.DimensionId.BubblePosition_Z)] = state.translation.z;
            vec[idx(types_1.DimensionId.BubbleVelocity_X)] = state.translation.vx;
            vec[idx(types_1.DimensionId.BubbleVelocity_Y)] = state.translation.vy;
            vec[idx(types_1.DimensionId.BubbleVelocity_Z)] = state.translation.vz;
        }
        else {
            // Default to zero if not provided
            vec[idx(types_1.DimensionId.BubblePosition_X)] = 0;
            vec[idx(types_1.DimensionId.BubblePosition_Y)] = 0;
            vec[idx(types_1.DimensionId.BubblePosition_Z)] = 0;
            vec[idx(types_1.DimensionId.BubbleVelocity_X)] = 0;
            vec[idx(types_1.DimensionId.BubbleVelocity_Y)] = 0;
            vec[idx(types_1.DimensionId.BubbleVelocity_Z)] = 0;
        }
        // Gas macro state
        vec[idx(types_1.DimensionId.GasPressure)] = state.gas.Pg;
        vec[idx(types_1.DimensionId.GasTemperature)] = state.gas.T;
        // Species number densities
        vec[idx(types_1.DimensionId.Species_H2O)] = state.species.numberDensity.H2O;
        vec[idx(types_1.DimensionId.Species_O2)] = state.species.numberDensity.O2;
        vec[idx(types_1.DimensionId.Species_N2)] = state.species.numberDensity.N2;
        vec[idx(types_1.DimensionId.Species_Ar)] = state.species.numberDensity.Ar;
        vec[idx(types_1.DimensionId.Species_Xe)] = state.species.numberDensity.Xe;
        vec[idx(types_1.DimensionId.Species_H)] = state.species.numberDensity.H;
        vec[idx(types_1.DimensionId.Species_O)] = state.species.numberDensity.O;
        vec[idx(types_1.DimensionId.Species_OH)] = state.species.numberDensity.OH || 0;
        vec[idx(types_1.DimensionId.Species_N)] = state.species.numberDensity.N || 0;
        // Plasma state
        vec[idx(types_1.DimensionId.ElectronDensity)] = state.plasma.ne;
        vec[idx(types_1.DimensionId.ElectronTemperature)] = state.plasma.Te;
        vec[idx(types_1.DimensionId.IonizationFraction)] = state.plasma.ionizationFraction;
        // Internal energy partitions
        vec[idx(types_1.DimensionId.InternalEnergy_Translational)] = state.internalEnergy.translational;
        vec[idx(types_1.DimensionId.InternalEnergy_Rotational)] = state.internalEnergy.rotational;
        vec[idx(types_1.DimensionId.InternalEnergy_Vibrational)] = state.internalEnergy.vibrational;
        vec[idx(types_1.DimensionId.InternalEnergy_Electronic)] = state.internalEnergy.electronic;
        // EM cavity modes
        const modes = state.em.modes;
        if (modes.length >= 1) {
            vec[idx(types_1.DimensionId.EmMode0_Re)] = modes[0].re;
            vec[idx(types_1.DimensionId.EmMode0_Im)] = modes[0].im;
        }
        if (modes.length >= 2) {
            vec[idx(types_1.DimensionId.EmMode1_Re)] = modes[1].re;
            vec[idx(types_1.DimensionId.EmMode1_Im)] = modes[1].im;
        }
        if (modes.length >= 3) {
            vec[idx(types_1.DimensionId.EmMode2_Re)] = modes[2].re;
            vec[idx(types_1.DimensionId.EmMode2_Im)] = modes[2].im;
        }
        // EM stored energy
        vec[idx(types_1.DimensionId.EmStoredEnergy)] = state.em.storedEnergy;
        // Acoustic phase
        vec[idx(types_1.DimensionId.AcousticPhase)] = state.acoustic.phase;
        // Reaction progress variables
        const xi = state.reactions.xi;
        if (xi.length >= 1)
            vec[idx(types_1.DimensionId.ReactionProgress_0)] = xi[0];
        if (xi.length >= 2)
            vec[idx(types_1.DimensionId.ReactionProgress_1)] = xi[1];
        if (xi.length >= 3)
            vec[idx(types_1.DimensionId.ReactionProgress_2)] = xi[2];
        return vec;
    }
    fromVector(vec, t) {
        if (vec.length !== this.layout.size) {
            throw new Error(`Vector length ${vec.length} does not match layout size ${this.layout.size}`);
        }
        const idx = (dim) => this.layout.indexOf(dim);
        // Extract state for EM frequency computation
        const R = vec[idx(types_1.DimensionId.Radius)];
        const ne = vec[idx(types_1.DimensionId.ElectronDensity)];
        // Compute EM mode frequencies from R(t) and plasma frequency
        // This is a simplified computation - full computation happens in emCavity module
        // Mode frequency scales as: ω_k(t) ≈ ω_k0 · (R0/R(t))
        // For now, use placeholder frequencies that will be updated by physics module
        // In a full implementation, we'd need access to EmCavityParams here
        const R0 = 1e-6; // Reference radius [m]
        const omega0_base = 2 * Math.PI * 3e14; // ~500 nm base frequency
        const omega1_base = 2 * Math.PI * 4e14; // ~400 nm
        const omega2_base = 2 * Math.PI * 5e14; // ~300 nm
        // Plasma frequency
        const e = 1.602176634e-19; // Electron charge [C]
        const epsilon_0 = 8.8541878128e-12; // Vacuum permittivity [F/m]
        const m_e = 9.1093837015e-31; // Electron mass [kg]
        const omega_p = Math.sqrt((ne * e * e) / (epsilon_0 * m_e));
        // Compute mode frequencies with plasma cutoff
        const computeModeFreq = (omega_base, R, omega_p) => {
            const omega_geometric = omega_base * (R0 / Math.max(R, 1e-10));
            if (omega_p >= omega_geometric)
                return 0; // Mode cut off
            return Math.sqrt(omega_geometric * omega_geometric - omega_p * omega_p);
        };
        const omega0 = computeModeFreq(omega0_base, R, omega_p);
        const omega1 = computeModeFreq(omega1_base, R, omega_p);
        const omega2 = computeModeFreq(omega2_base, R, omega_p);
        // Extract shape oscillations
        const shape = {
            a2: vec[idx(types_1.DimensionId.ShapeMode2_Amplitude)],
            a2_dot: vec[idx(types_1.DimensionId.ShapeMode2_Velocity)],
            a4: vec[idx(types_1.DimensionId.ShapeMode4_Amplitude)],
            a4_dot: vec[idx(types_1.DimensionId.ShapeMode4_Velocity)],
        };
        // Extract bubble translation
        const translation = {
            x: vec[idx(types_1.DimensionId.BubblePosition_X)],
            y: vec[idx(types_1.DimensionId.BubblePosition_Y)],
            z: vec[idx(types_1.DimensionId.BubblePosition_Z)],
            vx: vec[idx(types_1.DimensionId.BubbleVelocity_X)],
            vy: vec[idx(types_1.DimensionId.BubbleVelocity_Y)],
            vz: vec[idx(types_1.DimensionId.BubbleVelocity_Z)],
        };
        return {
            t,
            hydro: {
                R: vec[idx(types_1.DimensionId.Radius)],
                Rdot: vec[idx(types_1.DimensionId.RadiusVelocity)],
            },
            shape,
            translation,
            gas: {
                Pg: vec[idx(types_1.DimensionId.GasPressure)],
                T: vec[idx(types_1.DimensionId.GasTemperature)],
            },
            species: {
                numberDensity: {
                    H2O: vec[idx(types_1.DimensionId.Species_H2O)],
                    O2: vec[idx(types_1.DimensionId.Species_O2)],
                    N2: vec[idx(types_1.DimensionId.Species_N2)],
                    Ar: vec[idx(types_1.DimensionId.Species_Ar)],
                    Xe: vec[idx(types_1.DimensionId.Species_Xe)],
                    H: vec[idx(types_1.DimensionId.Species_H)],
                    O: vec[idx(types_1.DimensionId.Species_O)],
                    OH: vec[idx(types_1.DimensionId.Species_OH)],
                    N: vec[idx(types_1.DimensionId.Species_N)],
                },
            },
            plasma: {
                ne: vec[idx(types_1.DimensionId.ElectronDensity)],
                Te: vec[idx(types_1.DimensionId.ElectronTemperature)],
                ionizationFraction: vec[idx(types_1.DimensionId.IonizationFraction)],
            },
            internalEnergy: {
                translational: vec[idx(types_1.DimensionId.InternalEnergy_Translational)],
                rotational: vec[idx(types_1.DimensionId.InternalEnergy_Rotational)],
                vibrational: vec[idx(types_1.DimensionId.InternalEnergy_Vibrational)],
                electronic: vec[idx(types_1.DimensionId.InternalEnergy_Electronic)],
            },
            em: {
                modes: [
                    {
                        omega: omega0,
                        re: vec[idx(types_1.DimensionId.EmMode0_Re)],
                        im: vec[idx(types_1.DimensionId.EmMode0_Im)],
                    },
                    {
                        omega: omega1,
                        re: vec[idx(types_1.DimensionId.EmMode1_Re)],
                        im: vec[idx(types_1.DimensionId.EmMode1_Im)],
                    },
                    {
                        omega: omega2,
                        re: vec[idx(types_1.DimensionId.EmMode2_Re)],
                        im: vec[idx(types_1.DimensionId.EmMode2_Im)],
                    },
                ],
                storedEnergy: vec[idx(types_1.DimensionId.EmStoredEnergy)],
            },
            acoustic: {
                phase: vec[idx(types_1.DimensionId.AcousticPhase)],
            },
            reactions: {
                xi: [
                    vec[idx(types_1.DimensionId.ReactionProgress_0)],
                    vec[idx(types_1.DimensionId.ReactionProgress_1)],
                    vec[idx(types_1.DimensionId.ReactionProgress_2)],
                ],
            },
        };
    }
}
exports.DefaultStateVectorMapper = DefaultStateVectorMapper;
//# sourceMappingURL=statevector.js.map