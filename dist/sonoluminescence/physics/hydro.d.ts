import { BubbleFullState } from "../model/types";
export interface HydroParams {
    rho: number;
    mu: number;
    sigma: number;
    Pv: number;
    P0: number;
    c: number;
    useKellerMiksis?: boolean;
}
export interface HydroDerivatives {
    dRdt: number;
    dRdotdt: number;
}
/**
 * Compute hydrodynamic derivatives using Rayleigh-Plesset or Keller-Miksis equation
 *
 * Rayleigh-Plesset (default):
 * R * Rddot + (3/2) * Rdot² = (1/rho) * (Pg - P0 - Pacoustic - Pv - 2*sigma/R - 4*mu*Rdot/R)
 *
 * Keller-Miksis (enhanced, accounts for liquid compressibility):
 * (1 - Rdot/c) * R * Rddot + (3/2) * (1 - Rdot/(3*c)) * Rdot²
 *   = (1/rho) * (1 + Rdot/c + R/c * d/dt) * (Pg - P0 - Pacoustic - Pv - 2*sigma/R - 4*mu*Rdot/R)
 *
 * The Keller-Miksis equation is more accurate for high-speed collapse (Rdot comparable to sound speed).
 */
export declare function computeHydroDerivatives(state: BubbleFullState, params: HydroParams, Pacoustic: number, gamma?: number): HydroDerivatives;
//# sourceMappingURL=hydro.d.ts.map