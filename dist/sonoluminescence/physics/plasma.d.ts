import { BubbleFullState } from "../model/types";
export interface PlasmaParams {
    ionizationPotential_Ar: number;
    ionizationPotential_Xe: number;
    ionizationPotential_H: number;
    ionizationPotential_O: number;
    recombinationCoeff: number;
    electronCollisionFreq: number;
    useDetailedCollisions?: boolean;
    collisionModel?: 'constant' | 'temperature-dependent' | 'density-dependent';
    includeRadiationLoss?: boolean;
    radiationCoeff?: number;
    useTemperatureDependentCrossSections?: boolean;
    ionizationCrossSectionRef?: number;
}
export interface PlasmaDerivatives {
    dn_edt: number;
    dTedT: number;
    dIonFracDt: number;
}
/**
 * Compute plasma derivatives
 *
 * Electron density evolution:
 * - Ionization (temperature-dependent, via Saha equation)
 * - Recombination (proportional to n_e²)
 *
 * Electron temperature:
 * - Heating from compression work
 * - Cooling from collisions and radiation
 *
 * Ionization fraction:
 * - Average degree of ionization across all species
 */
export declare function computePlasmaDerivatives(state: BubbleFullState, params: PlasmaParams): PlasmaDerivatives;
//# sourceMappingURL=plasma.d.ts.map