// physics/hydro.ts
import { BubbleFullState } from "../model/types";

export interface HydroParams {
  rho: number;   // fluid density [kg/m³]
  mu: number;    // viscosity [Pa·s]
  sigma: number; // surface tension [N/m]
  Pv: number;    // vapor pressure [Pa]
  P0: number;    // ambient pressure [Pa]
  c: number;     // speed of sound in liquid [m/s] (for Keller-Miksis)
  useKellerMiksis?: boolean; // Use Keller-Miksis instead of Rayleigh-Plesset
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
export function computeHydroDerivatives(
  state: BubbleFullState,
  params: HydroParams,
  Pacoustic: number,
  gamma?: number // Adiabatic index from thermo params (optional for backward compatibility)
): HydroDerivatives {
  const { R, Rdot } = state.hydro;
  const { Pg } = state.gas;
  const { rho, mu, sigma, Pv, P0, c = 1482, useKellerMiksis = false } = params;

  // Avoid division by zero for very small R
  const R_safe = Math.max(R, 1e-10);

  // dR/dt = Rdot (kinematics)
  const dRdt = Rdot;

  // Pressure terms
  const pressureTerm = Pg - P0 - Pacoustic - Pv;
  const surfaceTensionTerm = 2 * sigma / R_safe;
  const viscousTerm = 4 * mu * Rdot / R_safe;
  const netPressure = pressureTerm - surfaceTensionTerm - viscousTerm;

  if (useKellerMiksis && c > 0) {
    // Keller-Miksis equation (accounts for liquid compressibility)
    // More accurate for supersonic collapse (Rdot ~ c)
    
    // Mach number
    const M = Rdot / c;
    
    // Compressibility correction factors
    const factor1 = 1.0 - M; // (1 - Rdot/c)
    const factor2 = 1.0 - M / 3.0; // (1 - Rdot/(3*c))
    
    // Time derivative of pressure term (simplified: assume dPg/dt dominates)
    // dPg/dt from adiabatic compression: -gamma * Pg * (3*Rdot/R)
    // Use provided gamma or fallback to 1.4 (monatomic/diatomic average)
    const effectiveGamma = gamma || 1.4;
    const dPgdt = -effectiveGamma * Pg * (3.0 * Rdot / R_safe);
    const dPdt_term = (R_safe / c) * dPgdt;
    
    // Modified pressure term with compressibility
    const modifiedPressure = (1.0 + M + dPdt_term) * netPressure;
    
    // Keller-Miksis: (1 - Rdot/c) * R * Rddot + (3/2) * (1 - Rdot/(3*c)) * Rdot² = (1/rho) * modifiedPressure
    // Solving for Rddot:
    const dRdotdt = (1.0 / (factor1 * R_safe)) * (
      (1.0 / rho) * modifiedPressure - (3.0 / 2.0) * factor2 * Rdot * Rdot
    );
    
    return { dRdt, dRdotdt };
  } else {
    // Standard Rayleigh-Plesset equation
    const dRdotdt = (1.0 / R_safe) * (
      (1.0 / rho) * netPressure - (3.0 / 2.0) * Rdot * Rdot
    );
    
    return { dRdt, dRdotdt };
  }
}
