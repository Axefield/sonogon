// statevector.ts

import { DimensionId, BubbleFullState, SpeciesId } from "../model/types";

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
export class DefaultStateVectorLayout implements StateVectorLayout {
  private readonly dimensionMap: Map<DimensionId, number>;
  public readonly size: number;

  constructor() {
    this.dimensionMap = new Map();
    let index = 0;

    // Hydrodynamic / geometric (2)
    this.dimensionMap.set(DimensionId.Radius, index++);
    this.dimensionMap.set(DimensionId.RadiusVelocity, index++);
    
    // Shape oscillations (4) - optional
    this.dimensionMap.set(DimensionId.ShapeMode2_Amplitude, index++);
    this.dimensionMap.set(DimensionId.ShapeMode2_Velocity, index++);
    this.dimensionMap.set(DimensionId.ShapeMode4_Amplitude, index++);
    this.dimensionMap.set(DimensionId.ShapeMode4_Velocity, index++);
    
    // Bubble translation (6) - optional
    this.dimensionMap.set(DimensionId.BubblePosition_X, index++);
    this.dimensionMap.set(DimensionId.BubblePosition_Y, index++);
    this.dimensionMap.set(DimensionId.BubblePosition_Z, index++);
    this.dimensionMap.set(DimensionId.BubbleVelocity_X, index++);
    this.dimensionMap.set(DimensionId.BubbleVelocity_Y, index++);
    this.dimensionMap.set(DimensionId.BubbleVelocity_Z, index++);

    // Gas macro state (2)
    this.dimensionMap.set(DimensionId.GasPressure, index++);
    this.dimensionMap.set(DimensionId.GasTemperature, index++);

    // Species number densities (9)
    this.dimensionMap.set(DimensionId.Species_H2O, index++);
    this.dimensionMap.set(DimensionId.Species_O2, index++);
    this.dimensionMap.set(DimensionId.Species_N2, index++);
    this.dimensionMap.set(DimensionId.Species_Ar, index++);
    this.dimensionMap.set(DimensionId.Species_Xe, index++);
    this.dimensionMap.set(DimensionId.Species_H, index++);
    this.dimensionMap.set(DimensionId.Species_O, index++);
    this.dimensionMap.set(DimensionId.Species_OH, index++);
    this.dimensionMap.set(DimensionId.Species_N, index++);

    // Plasma state (3)
    this.dimensionMap.set(DimensionId.ElectronDensity, index++);
    this.dimensionMap.set(DimensionId.ElectronTemperature, index++);
    this.dimensionMap.set(DimensionId.IonizationFraction, index++);

    // Internal energy partitions (4)
    this.dimensionMap.set(DimensionId.InternalEnergy_Translational, index++);
    this.dimensionMap.set(DimensionId.InternalEnergy_Rotational, index++);
    this.dimensionMap.set(DimensionId.InternalEnergy_Vibrational, index++);
    this.dimensionMap.set(DimensionId.InternalEnergy_Electronic, index++);

    // EM cavity modes (6 for 3 modes)
    this.dimensionMap.set(DimensionId.EmMode0_Re, index++);
    this.dimensionMap.set(DimensionId.EmMode0_Im, index++);
    this.dimensionMap.set(DimensionId.EmMode1_Re, index++);
    this.dimensionMap.set(DimensionId.EmMode1_Im, index++);
    this.dimensionMap.set(DimensionId.EmMode2_Re, index++);
    this.dimensionMap.set(DimensionId.EmMode2_Im, index++);

    // EM stored energy (1)
    this.dimensionMap.set(DimensionId.EmStoredEnergy, index++);

    // Acoustic phase (1)
    this.dimensionMap.set(DimensionId.AcousticPhase, index++);

    // Reaction progress variables (3)
    this.dimensionMap.set(DimensionId.ReactionProgress_0, index++);
    this.dimensionMap.set(DimensionId.ReactionProgress_1, index++);
    this.dimensionMap.set(DimensionId.ReactionProgress_2, index++);

    this.size = index;
  }

  indexOf(dim: DimensionId): number {
    const idx = this.dimensionMap.get(dim);
    if (idx === undefined) {
      throw new Error(`DimensionId ${dim} not found in layout`);
    }
    return idx;
  }
}

/**
 * Concrete implementation of StateVectorMapper
 * Converts between BubbleFullState and Float64Array
 */
export class DefaultStateVectorMapper implements StateVectorMapper {
  public readonly layout: StateVectorLayout;
  private readonly speciesOrder: SpeciesId[] = ["H2O", "O2", "N2", "Ar", "Xe", "H", "O"];

  constructor(layout?: StateVectorLayout) {
    this.layout = layout || new DefaultStateVectorLayout();
  }

  toVector(state: BubbleFullState): Float64Array {
    const vec = new Float64Array(this.layout.size);
    const idx = (dim: DimensionId) => this.layout.indexOf(dim);

    // Hydrodynamic
    vec[idx(DimensionId.Radius)] = state.hydro.R;
    vec[idx(DimensionId.RadiusVelocity)] = state.hydro.Rdot;

    // Gas macro state
    vec[idx(DimensionId.GasPressure)] = state.gas.Pg;
    vec[idx(DimensionId.GasTemperature)] = state.gas.T;

    // Species number densities
    vec[idx(DimensionId.Species_H2O)] = state.species.numberDensity.H2O;
    vec[idx(DimensionId.Species_O2)] = state.species.numberDensity.O2;
    vec[idx(DimensionId.Species_N2)] = state.species.numberDensity.N2;
    vec[idx(DimensionId.Species_Ar)] = state.species.numberDensity.Ar;
    vec[idx(DimensionId.Species_Xe)] = state.species.numberDensity.Xe;
    vec[idx(DimensionId.Species_H)] = state.species.numberDensity.H;
    vec[idx(DimensionId.Species_O)] = state.species.numberDensity.O;
    vec[idx(DimensionId.Species_OH)] = state.species.numberDensity.OH || 0;
    vec[idx(DimensionId.Species_N)] = state.species.numberDensity.N || 0;

    // Plasma state
    vec[idx(DimensionId.ElectronDensity)] = state.plasma.ne;
    vec[idx(DimensionId.ElectronTemperature)] = state.plasma.Te;
    vec[idx(DimensionId.IonizationFraction)] = state.plasma.ionizationFraction;

    // Internal energy partitions
    vec[idx(DimensionId.InternalEnergy_Translational)] = state.internalEnergy.translational;
    vec[idx(DimensionId.InternalEnergy_Rotational)] = state.internalEnergy.rotational;
    vec[idx(DimensionId.InternalEnergy_Vibrational)] = state.internalEnergy.vibrational;
    vec[idx(DimensionId.InternalEnergy_Electronic)] = state.internalEnergy.electronic;

    // EM cavity modes
    const modes = state.em.modes;
    if (modes.length >= 1) {
      vec[idx(DimensionId.EmMode0_Re)] = modes[0].re;
      vec[idx(DimensionId.EmMode0_Im)] = modes[0].im;
    }
    if (modes.length >= 2) {
      vec[idx(DimensionId.EmMode1_Re)] = modes[1].re;
      vec[idx(DimensionId.EmMode1_Im)] = modes[1].im;
    }
    if (modes.length >= 3) {
      vec[idx(DimensionId.EmMode2_Re)] = modes[2].re;
      vec[idx(DimensionId.EmMode2_Im)] = modes[2].im;
    }

    // EM stored energy
    vec[idx(DimensionId.EmStoredEnergy)] = state.em.storedEnergy;

    // Acoustic phase
    vec[idx(DimensionId.AcousticPhase)] = state.acoustic.phase;

    // Reaction progress variables
    const xi = state.reactions.xi;
    if (xi.length >= 1) vec[idx(DimensionId.ReactionProgress_0)] = xi[0];
    if (xi.length >= 2) vec[idx(DimensionId.ReactionProgress_1)] = xi[1];
    if (xi.length >= 3) vec[idx(DimensionId.ReactionProgress_2)] = xi[2];

    return vec;
  }

  fromVector(vec: Float64Array, t: number): BubbleFullState {
    if (vec.length !== this.layout.size) {
      throw new Error(
        `Vector length ${vec.length} does not match layout size ${this.layout.size}`
      );
    }

    const idx = (dim: DimensionId) => this.layout.indexOf(dim);

    // Extract state for EM frequency computation
    const R = vec[idx(DimensionId.Radius)];
    const ne = vec[idx(DimensionId.ElectronDensity)];
    
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
    const computeModeFreq = (omega_base: number, R: number, omega_p: number): number => {
      const omega_geometric = omega_base * (R0 / Math.max(R, 1e-10));
      if (omega_p >= omega_geometric) return 0; // Mode cut off
      return Math.sqrt(omega_geometric * omega_geometric - omega_p * omega_p);
    };
    
    const omega0 = computeModeFreq(omega0_base, R, omega_p);
    const omega1 = computeModeFreq(omega1_base, R, omega_p);
    const omega2 = computeModeFreq(omega2_base, R, omega_p);

    return {
      t,
      hydro: {
        R: vec[idx(DimensionId.Radius)],
        Rdot: vec[idx(DimensionId.RadiusVelocity)],
      },
      gas: {
        Pg: vec[idx(DimensionId.GasPressure)],
        T: vec[idx(DimensionId.GasTemperature)],
      },
      species: {
        numberDensity: {
          H2O: vec[idx(DimensionId.Species_H2O)],
          O2: vec[idx(DimensionId.Species_O2)],
          N2: vec[idx(DimensionId.Species_N2)],
          Ar: vec[idx(DimensionId.Species_Ar)],
          Xe: vec[idx(DimensionId.Species_Xe)],
          H: vec[idx(DimensionId.Species_H)],
          O: vec[idx(DimensionId.Species_O)],
          OH: vec[idx(DimensionId.Species_OH)],
          N: vec[idx(DimensionId.Species_N)],
        },
      },
      plasma: {
        ne: vec[idx(DimensionId.ElectronDensity)],
        Te: vec[idx(DimensionId.ElectronTemperature)],
        ionizationFraction: vec[idx(DimensionId.IonizationFraction)],
      },
      internalEnergy: {
        translational: vec[idx(DimensionId.InternalEnergy_Translational)],
        rotational: vec[idx(DimensionId.InternalEnergy_Rotational)],
        vibrational: vec[idx(DimensionId.InternalEnergy_Vibrational)],
        electronic: vec[idx(DimensionId.InternalEnergy_Electronic)],
      },
      em: {
        modes: [
          {
            omega: omega0,
            re: vec[idx(DimensionId.EmMode0_Re)],
            im: vec[idx(DimensionId.EmMode0_Im)],
          },
          {
            omega: omega1,
            re: vec[idx(DimensionId.EmMode1_Re)],
            im: vec[idx(DimensionId.EmMode1_Im)],
          },
          {
            omega: omega2,
            re: vec[idx(DimensionId.EmMode2_Re)],
            im: vec[idx(DimensionId.EmMode2_Im)],
          },
        ],
        storedEnergy: vec[idx(DimensionId.EmStoredEnergy)],
      },
      acoustic: {
        phase: vec[idx(DimensionId.AcousticPhase)],
      },
      reactions: {
        xi: [
          vec[idx(DimensionId.ReactionProgress_0)],
          vec[idx(DimensionId.ReactionProgress_1)],
          vec[idx(DimensionId.ReactionProgress_2)],
        ],
      },
    };
  }
}