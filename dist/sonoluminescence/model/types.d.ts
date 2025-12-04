export declare enum DimensionId {
    Radius = 0,// R(t)
    RadiusVelocity = 1,// Rdot(t)
    ShapeMode2_Amplitude = 2,// a₂ (quadrupole)
    ShapeMode2_Velocity = 3,// a₂_dot
    ShapeMode4_Amplitude = 4,// a₄ (hexadecapole)
    ShapeMode4_Velocity = 5,// a₄_dot
    BubblePosition_X = 6,// x position
    BubblePosition_Y = 7,// y position
    BubblePosition_Z = 8,// z position
    BubbleVelocity_X = 9,// vx
    BubbleVelocity_Y = 10,// vy
    BubbleVelocity_Z = 11,// vz
    GasPressure = 12,// Pg
    GasTemperature = 13,// Tg
    Species_H2O = 14,
    Species_O2 = 15,
    Species_N2 = 16,
    Species_Ar = 17,
    Species_Xe = 18,
    Species_H = 19,
    Species_O = 20,
    Species_OH = 21,// Hydroxyl radical
    Species_N = 22,// Atomic nitrogen
    ElectronDensity = 23,
    ElectronTemperature = 24,
    IonizationFraction = 25,// average ionization degree or similar
    InternalEnergy_Translational = 26,
    InternalEnergy_Rotational = 27,
    InternalEnergy_Vibrational = 28,
    InternalEnergy_Electronic = 29,
    EmMode0_Re = 30,
    EmMode0_Im = 31,
    EmMode1_Re = 32,
    EmMode1_Im = 33,
    EmMode2_Re = 34,
    EmMode2_Im = 35,
    EmStoredEnergy = 36,
    AcousticPhase = 37,
    ReactionProgress_0 = 38,
    ReactionProgress_1 = 39,
    ReactionProgress_2 = 40
}
export interface HydroState {
    R: number;
    Rdot: number;
}
export interface ShapeOscillationState {
    a2: number;
    a2_dot: number;
    a4: number;
    a4_dot: number;
}
export interface BubbleTranslationState {
    x: number;
    y: number;
    z: number;
    vx: number;
    vy: number;
    vz: number;
}
export interface GasMacroState {
    Pg: number;
    T: number;
}
export type SpeciesId = "H2O" | "O2" | "N2" | "Ar" | "Xe" | "H" | "O" | "OH" | "N";
export interface SpeciesState {
    numberDensity: Record<SpeciesId, number>;
}
export interface PlasmaState {
    ne: number;
    Te: number;
    ionizationFraction: number;
}
export interface InternalEnergyState {
    translational: number;
    rotational: number;
    vibrational: number;
    electronic: number;
}
export interface EmModeState {
    omega: number;
    re: number;
    im: number;
}
export interface EmFieldState {
    modes: EmModeState[];
    storedEnergy: number;
}
export interface AcousticState {
    phase: number;
}
export interface ReactionProgressState {
    xi: number[];
}
export interface BubbleFullState {
    t: number;
    hydro: HydroState;
    shape?: ShapeOscillationState;
    translation?: BubbleTranslationState;
    gas: GasMacroState;
    species: SpeciesState;
    plasma: PlasmaState;
    internalEnergy: InternalEnergyState;
    em: EmFieldState;
    acoustic: AcousticState;
    reactions: ReactionProgressState;
}
//# sourceMappingURL=types.d.ts.map