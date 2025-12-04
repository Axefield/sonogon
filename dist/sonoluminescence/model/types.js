"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.DimensionId = void 0;
// model/types.ts
var DimensionId;
(function (DimensionId) {
    // Hydrodynamic / geometric
    DimensionId[DimensionId["Radius"] = 0] = "Radius";
    DimensionId[DimensionId["RadiusVelocity"] = 1] = "RadiusVelocity";
    // Shape oscillations (spherical harmonics)
    DimensionId[DimensionId["ShapeMode2_Amplitude"] = 2] = "ShapeMode2_Amplitude";
    DimensionId[DimensionId["ShapeMode2_Velocity"] = 3] = "ShapeMode2_Velocity";
    DimensionId[DimensionId["ShapeMode4_Amplitude"] = 4] = "ShapeMode4_Amplitude";
    DimensionId[DimensionId["ShapeMode4_Velocity"] = 5] = "ShapeMode4_Velocity";
    // Bubble translation
    DimensionId[DimensionId["BubblePosition_X"] = 6] = "BubblePosition_X";
    DimensionId[DimensionId["BubblePosition_Y"] = 7] = "BubblePosition_Y";
    DimensionId[DimensionId["BubblePosition_Z"] = 8] = "BubblePosition_Z";
    DimensionId[DimensionId["BubbleVelocity_X"] = 9] = "BubbleVelocity_X";
    DimensionId[DimensionId["BubbleVelocity_Y"] = 10] = "BubbleVelocity_Y";
    DimensionId[DimensionId["BubbleVelocity_Z"] = 11] = "BubbleVelocity_Z";
    // Gas macro state
    DimensionId[DimensionId["GasPressure"] = 12] = "GasPressure";
    DimensionId[DimensionId["GasTemperature"] = 13] = "GasTemperature";
    // Species number densities (example set)
    DimensionId[DimensionId["Species_H2O"] = 14] = "Species_H2O";
    DimensionId[DimensionId["Species_O2"] = 15] = "Species_O2";
    DimensionId[DimensionId["Species_N2"] = 16] = "Species_N2";
    DimensionId[DimensionId["Species_Ar"] = 17] = "Species_Ar";
    DimensionId[DimensionId["Species_Xe"] = 18] = "Species_Xe";
    DimensionId[DimensionId["Species_H"] = 19] = "Species_H";
    DimensionId[DimensionId["Species_O"] = 20] = "Species_O";
    DimensionId[DimensionId["Species_OH"] = 21] = "Species_OH";
    DimensionId[DimensionId["Species_N"] = 22] = "Species_N";
    // Plasma state
    DimensionId[DimensionId["ElectronDensity"] = 23] = "ElectronDensity";
    DimensionId[DimensionId["ElectronTemperature"] = 24] = "ElectronTemperature";
    DimensionId[DimensionId["IonizationFraction"] = 25] = "IonizationFraction";
    // Internal energy partitions
    DimensionId[DimensionId["InternalEnergy_Translational"] = 26] = "InternalEnergy_Translational";
    DimensionId[DimensionId["InternalEnergy_Rotational"] = 27] = "InternalEnergy_Rotational";
    DimensionId[DimensionId["InternalEnergy_Vibrational"] = 28] = "InternalEnergy_Vibrational";
    DimensionId[DimensionId["InternalEnergy_Electronic"] = 29] = "InternalEnergy_Electronic";
    // EM cavity modes (example = 3 modes, each Re/Im)
    DimensionId[DimensionId["EmMode0_Re"] = 30] = "EmMode0_Re";
    DimensionId[DimensionId["EmMode0_Im"] = 31] = "EmMode0_Im";
    DimensionId[DimensionId["EmMode1_Re"] = 32] = "EmMode1_Re";
    DimensionId[DimensionId["EmMode1_Im"] = 33] = "EmMode1_Im";
    DimensionId[DimensionId["EmMode2_Re"] = 34] = "EmMode2_Re";
    DimensionId[DimensionId["EmMode2_Im"] = 35] = "EmMode2_Im";
    // EM “stored energy” scalar (for your mist/negative-space pump)
    DimensionId[DimensionId["EmStoredEnergy"] = 36] = "EmStoredEnergy";
    // Acoustic phase tracking
    DimensionId[DimensionId["AcousticPhase"] = 37] = "AcousticPhase";
    // Optional: chemical progress variables (reaction extents)
    DimensionId[DimensionId["ReactionProgress_0"] = 38] = "ReactionProgress_0";
    DimensionId[DimensionId["ReactionProgress_1"] = 39] = "ReactionProgress_1";
    DimensionId[DimensionId["ReactionProgress_2"] = 40] = "ReactionProgress_2";
    // Add more as needed…
})(DimensionId || (exports.DimensionId = DimensionId = {}));
//# sourceMappingURL=types.js.map