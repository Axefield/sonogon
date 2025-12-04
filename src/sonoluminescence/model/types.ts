// model/types.ts
export enum DimensionId {
    // Hydrodynamic / geometric
    Radius,            // R(t)
    RadiusVelocity,    // Rdot(t)
    
    // Shape oscillations (spherical harmonics)
    ShapeMode2_Amplitude,      // a₂ (quadrupole)
    ShapeMode2_Velocity,       // a₂_dot
    ShapeMode4_Amplitude,      // a₄ (hexadecapole)
    ShapeMode4_Velocity,       // a₄_dot
    
    // Bubble translation
    BubblePosition_X,          // x position
    BubblePosition_Y,          // y position
    BubblePosition_Z,          // z position
    BubbleVelocity_X,          // vx
    BubbleVelocity_Y,          // vy
    BubbleVelocity_Z,          // vz
  
    // Gas macro state
    GasPressure,       // Pg
    GasTemperature,    // Tg
  
    // Species number densities (example set)
    Species_H2O,
    Species_O2,
    Species_N2,
    Species_Ar,
    Species_Xe,
    Species_H,
    Species_O,
    Species_OH,  // Hydroxyl radical
    Species_N,   // Atomic nitrogen
  
    // Plasma state
    ElectronDensity,
    ElectronTemperature,
    IonizationFraction,  // average ionization degree or similar
  
    // Internal energy partitions
    InternalEnergy_Translational,
    InternalEnergy_Rotational,
    InternalEnergy_Vibrational,
    InternalEnergy_Electronic,
  
    // EM cavity modes (example = 3 modes, each Re/Im)
    EmMode0_Re,
    EmMode0_Im,
    EmMode1_Re,
    EmMode1_Im,
    EmMode2_Re,
    EmMode2_Im,
  
    // EM “stored energy” scalar (for your mist/negative-space pump)
    EmStoredEnergy,
  
    // Acoustic phase tracking
    AcousticPhase,
  
    // Optional: chemical progress variables (reaction extents)
    ReactionProgress_0,
    ReactionProgress_1,
    ReactionProgress_2,
  
    // Add more as needed…

    
  }
  export interface HydroState {
    R: number;
    Rdot: number;
  }
  
  export interface ShapeOscillationState {
    a2: number;      // Quadrupole mode amplitude [m]
    a2_dot: number;  // Quadrupole mode velocity [m/s]
    a4: number;      // Hexadecapole mode amplitude [m]
    a4_dot: number;  // Hexadecapole mode velocity [m/s]
  }
  
  export interface BubbleTranslationState {
    x: number;       // Position [m]
    y: number;
    z: number;
    vx: number;      // Velocity [m/s]
    vy: number;
    vz: number;
  }
  
  export interface GasMacroState {
    Pg: number;
    T: number;
  }
  
  export type SpeciesId =
    | "H2O"
    | "O2"
    | "N2"
    | "Ar"
    | "Xe"
    | "H"
    | "O"
    | "OH"  // Hydroxyl radical
    | "N";   // Atomic nitrogen
  
  export interface SpeciesState {
    // number density [1/m³] per species
    numberDensity: Record<SpeciesId, number>;
  }
  
  export interface PlasmaState {
    ne: number;      // electron density
    Te: number;      // electron temperature
    ionizationFraction: number; // e.g. <Z>/Z_max
  }
  
  export interface InternalEnergyState {
    translational: number;
    rotational: number;
    vibrational: number;
    electronic: number;
  }
  
  export interface EmModeState {
    omega: number;   // instantaneous mode frequency [rad/s]
    re: number;      // Re(a_k)
    im: number;      // Im(a_k)
  }
  
  export interface EmFieldState {
    modes: EmModeState[]; // small set of cavity modes
    storedEnergy: number; // your “negative space” E_em
  }
  
  export interface AcousticState {
    phase: number;   // φ(t) = ω t mod 2π
  }
  
  export interface ReactionProgressState {
    xi: number[];    // reaction extent variables
  }
  
  export interface BubbleFullState {
    t: number;
    hydro: HydroState;
    shape?: ShapeOscillationState;  // Optional: shape oscillations
    translation?: BubbleTranslationState;  // Optional: bubble translation
    gas: GasMacroState;
    species: SpeciesState;
    plasma: PlasmaState;
    internalEnergy: InternalEnergyState;
    em: EmFieldState;
    acoustic: AcousticState;
    reactions: ReactionProgressState;
  }
  