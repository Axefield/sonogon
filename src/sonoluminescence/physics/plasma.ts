// physics/plasma.ts
import { BubbleFullState } from "../model/types";
import { Constants, Calculations } from "../core/units";

export interface PlasmaParams {
  // Ionization potentials [J] for different species
  ionizationPotential_Ar: number; // Argon: 15.76 eV = 2.525e-18 J
  ionizationPotential_Xe: number; // Xenon: 12.13 eV = 1.944e-18 J
  ionizationPotential_H: number;   // Hydrogen: 13.6 eV = 2.179e-18 J
  ionizationPotential_O: number;   // Oxygen: 13.62 eV = 2.183e-18 J
  
  // Recombination coefficient [m³/s]
  recombinationCoeff: number; // ~1e-13 m³/s for typical plasmas
  
  // Electron collision frequency [1/s]
  electronCollisionFreq: number; // ~1e12-1e14 1/s for dense plasmas
  
  // Enhanced collision models
  useDetailedCollisions?: boolean; // Use temperature-dependent collision rates
  collisionModel?: 'constant' | 'temperature-dependent' | 'density-dependent';
  
  // Radiation transport
  includeRadiationLoss?: boolean; // Include radiative energy loss
  radiationCoeff?: number; // Radiation coefficient [W/(m³·K⁴)]
  
  // Collision cross sections
  useTemperatureDependentCrossSections?: boolean; // Temperature-dependent σ(E)
  ionizationCrossSectionRef?: number; // Reference ionization cross section [m²]
}

export interface PlasmaDerivatives {
  dn_edt: number;
  dTedT: number;
  dIonFracDt: number;
}

/**
 * Saha equation for ionization equilibrium
 * 
 * For a species with ionization potential I:
 * n_e * n_i / n_0 = (2 * g_i / g_0) * (2*π*m_e*k_B*T / h²)^(3/2) * exp(-I / (k_B*T))
 * 
 * Where:
 * - n_e: electron density
 * - n_i: ion density
 * - n_0: neutral density
 * - g_i, g_0: statistical weights
 * - I: ionization potential
 */
function computeSahaEquilibrium(
  T: number,
  n_neutral: number,
  ionizationPotential: number
): number {
  // Statistical weights (simplified: g_i/g_0 ≈ 1 for ground states)
  const g_ratio = 1.0;
  
  // Prefactor: (2*π*m_e*k_B*T / h²)^(3/2)
  const prefactor = Math.pow(
    (2.0 * Math.PI * Constants.m_e * Constants.k_B * T) / (Constants.h * Constants.h),
    1.5
  );
  
  // Exponential: exp(-I / (k_B*T))
  const exponent = Math.exp(-ionizationPotential / (Constants.k_B * T));
  
  // Saha equation gives: n_e * n_i / n_0
  // For simplicity, assume n_i ≈ n_e (single ionization)
  // Then: n_e² / n_0 = prefactor * exp(-I/(k_B*T))
  // So: n_e = sqrt(n_0 * prefactor * exp(-I/(k_B*T)))
  const sahaProduct = 2.0 * g_ratio * prefactor * exponent;
  
  // Solve for equilibrium electron density
  // n_e² = n_0 * sahaProduct
  // But we need to account for total density conservation
  // Simplified: use n_neutral as reference
  const n_e_eq = Math.sqrt(Math.max(n_neutral, 0) * sahaProduct);
  
  return n_e_eq;
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
export function computePlasmaDerivatives(
  state: BubbleFullState,
  params: PlasmaParams
): PlasmaDerivatives {
  const { R, Rdot } = state.hydro;
  const { T } = state.gas;
  const { ne, Te } = state.plasma;
  const R_safe = Math.max(R, 1e-10);
  
  // Volume change rate
  const volumeChangeRate = 3.0 * Rdot / R_safe;

  // Compute equilibrium electron density from Saha equation
  // Use dominant species (Ar or Xe typically)
  const n_Ar = state.species.numberDensity.Ar;
  const n_Xe = state.species.numberDensity.Xe;
  const n_H = state.species.numberDensity.H;
  const n_O = state.species.numberDensity.O;
  
  // Weighted average ionization potential (simplified)
  const totalNeutral = n_Ar + n_Xe + n_H + n_O;
  let weightedIonizationPotential = 0;
  if (totalNeutral > 0) {
    weightedIonizationPotential = (
      n_Ar * params.ionizationPotential_Ar +
      n_Xe * params.ionizationPotential_Xe +
      n_H * params.ionizationPotential_H +
      n_O * params.ionizationPotential_O
    ) / totalNeutral;
  }
  
  // Equilibrium electron density from Saha equation
  const n_e_equilibrium = totalNeutral > 0
    ? computeSahaEquilibrium(T, totalNeutral, weightedIonizationPotential)
    : 0;

  // Ionization and recombination rates
  let ionizationRate = 0;
  let recombinationRate = 0;
  
  if (params.useDetailedCollisions) {
    // Detailed rate equations based on collision cross sections
    
    // Ionization rate: electron impact ionization
    // Rate = n_e * n_neutral * <σ_ion * v>
    // where <σ_ion * v> is the rate coefficient
    const v_th_e = Math.sqrt((Constants.k_B * Te) / Constants.m_e); // Electron thermal velocity
    const v_th_i = Math.sqrt((Constants.k_B * T) / (Constants.m_p * 40)); // Ion thermal velocity (Ar)
    
    // Ionization cross section (simplified: σ_ion ~ (E - I) / I for E > I)
    // Rate coefficient: k_ion = <σ_ion * v> ≈ σ_0 * v_th * exp(-I / (k_B * Te))
    const sigma_0 = 1e-20; // Reference cross section [m²]
    const ionizationRateCoeff = sigma_0 * v_th_e * Math.exp(-weightedIonizationPotential / (Constants.k_B * Te));
    ionizationRate = ne * totalNeutral * ionizationRateCoeff;
    
    // Recombination rates (multiple channels)
    // 1. Radiative recombination: e + A⁺ → A + photon
    const alpha_rad = 2.6e-19 * Math.pow(Te / 10000, -0.5); // Approximate [m³/s]
    const recombination_radiative = alpha_rad * ne * (totalNeutral * 0.1); // Assume 10% ions
    
    // 2. Three-body recombination: e + e + A⁺ → e + A
    const alpha_3body = 1e-32 * Math.pow(Te / 10000, -4.5); // Approximate [m⁶/s]
    const n_ions = totalNeutral * 0.1; // Simplified: assume 10% ionization
    const recombination_3body = alpha_3body * ne * ne * n_ions;
    
    // 3. Dielectronic recombination (for high Z): e + A⁺ → A**
    const alpha_diel = 1e-20 * Math.pow(Te / 10000, -1.5); // Approximate [m³/s]
    const recombination_diel = alpha_diel * ne * n_ions;
    
    recombinationRate = recombination_radiative + recombination_3body + recombination_diel;
  } else {
    // Simplified rates (original)
    // Ionization rate: drive toward Saha equilibrium
    const tau_ionization = 1e-9; // characteristic ionization time [s]
    ionizationRate = Math.max(0, (n_e_equilibrium - ne) / tau_ionization);
    
    // Recombination rate: ~ recombinationCoeff * n_e²
    recombinationRate = params.recombinationCoeff * ne * ne;
  }
  
  // Net electron density change
  // dN/dt from compression: -ne * volumeChangeRate
  // Plus ionization minus recombination
  const dn_edt = 
    -ne * volumeChangeRate + // compression/expansion
    ionizationRate - // ionization
    recombinationRate; // recombination

  // Electron temperature evolution
  // Heating from compression work (similar to gas temperature)
  const adiabaticHeating = (5.0 / 3.0 - 1.0) * Te * volumeChangeRate; // gamma_e ≈ 5/3 for electrons
  
  // Collisional cooling: energy loss to heavy particles
  let collisionalCooling = 0;
  if (params.useDetailedCollisions && params.collisionModel === 'temperature-dependent') {
    // Temperature-dependent collision frequency with cross sections
    const n_neutral = totalNeutral;
    const v_th = Math.sqrt((Constants.k_B * Te) / Constants.m_e);
    
    // Temperature-dependent collision cross section
    // σ(E) typically decreases with energy for elastic, increases for inelastic
    // Simplified: σ_coll ~ σ_0 * (E_th/E)^alpha where E = k_B*Te
    let sigma_coll = 1e-19; // Default [m²]
    if (params.useTemperatureDependentCrossSections) {
      const E_th = 1.0 * Constants.e; // Threshold energy ~1 eV
      const E_e = Constants.k_B * Te;
      const alpha = 0.5; // Power law exponent
      const sigma_0 = params.ionizationCrossSectionRef || 1e-19;
      sigma_coll = sigma_0 * Math.pow(E_th / Math.max(E_e, E_th * 0.1), alpha);
    }
    
    const nu_coll = n_neutral * sigma_coll * v_th;
    // Energy exchange: ~(2*m_e/m_i) * nu_coll * (Te - T) per collision
    const mass_ratio = Constants.m_e / (Constants.m_p * 40); // m_e / m_Ar
    collisionalCooling = -(2.0 * mass_ratio) * nu_coll * (Te - T);
  } else if (params.useDetailedCollisions && params.collisionModel === 'density-dependent') {
    // Density-dependent: ν_coll ~ n_e * n_i * <σv>
    const n_total = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0);
    const n_ions = n_total * 0.1; // Simplified
    const sigma_v = 1e-14; // Approximate <σv> [m³/s]
    const nu_coll = ne * n_ions * sigma_v;
    const mass_ratio = Constants.m_e / (Constants.m_p * 40);
    collisionalCooling = -(2.0 * mass_ratio) * nu_coll * (Te - T);
  } else {
    // Constant collision frequency (original)
    collisionalCooling = -params.electronCollisionFreq * (Te - T) * 0.1;
  }
  
  // Bremsstrahlung cooling (energy radiated away)
  const bremsstrahlungPower = Calculations.bremsstrahlungPowerDensity(ne, Te);
  const bremsstrahlungCooling = -bremsstrahlungPower / (ne * Constants.k_B * 1e6); // Simplified normalization
  
  // Radiation loss (if enabled)
  let radiationCooling = 0;
  if (params.includeRadiationLoss) {
    const radiationCoeff = params.radiationCoeff || 1e-30; // [W/(m³·K⁴)]
    const radiationPower = radiationCoeff * Te * Te * Te * Te; // T⁴ dependence
    radiationCooling = -radiationPower / (ne * Constants.k_B * 1e6);
  }
  
  const dTedT = adiabaticHeating + collisionalCooling + bremsstrahlungCooling + radiationCooling;

  // Ionization fraction evolution
  // Average ionization fraction: <Z> = ne / n_total
  const n_total = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0);
  const ionizationFraction_eq = n_total > 0 ? ne / n_total : 0;
  const currentIonizationFraction = state.plasma.ionizationFraction;
  
  // Relaxation toward equilibrium
  const tau_ionFrac = 1e-9;
  const dIonFracDt = (ionizationFraction_eq - currentIonizationFraction) / tau_ionFrac;

  return { dn_edt, dTedT, dIonFracDt };
}
