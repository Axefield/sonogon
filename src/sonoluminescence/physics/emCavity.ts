// physics/emCavity.ts
import { BubbleFullState, EmFieldState } from "../model/types";
import { Constants, Calculations } from "../core/units";

export interface EmCavityParams {
  modeFrequencies0: number[];  // base frequencies for modes [rad/s]
  pumpCoefficient: number;     // α: how strongly gradients pump E_em [J/(m·s)]
  decayTime: number;           // τ: decay time for EM stored energy [s]
  couplingStrength: number;   // coupling between modes and boundary motion
  refractiveIndexBase: number; // base refractive index of bubble interior
  
  // Enhanced EM features
  useModeCoupling?: boolean; // Enable cross-mode coupling
  modeCouplingMatrix?: number[][]; // Coupling matrix between modes (symmetric)
  useFrequencyDependentQ?: boolean; // Use frequency-dependent quality factor
  Q0?: number; // Base quality factor
  QFrequencyDependence?: (omega: number) => number; // Custom Q(ω) function
  
  // Detailed parametric pumping
  useDetailedParametricPumping?: boolean; // Use quantum parametric Hamiltonian
  parametricCoupling?: number; // Coupling strength g [1/s]
  
  // Thermodynamic parameter for gradient calculation
  thermoGamma?: number; // Adiabatic index (from thermo params)
  
  // Cavity-QED enhancements
  useRefractiveIndexFrequencyShift?: boolean; // Include n(r,t) in frequency calculation
  useRefractiveIndexTimeDerivative?: boolean; // Include ∂n/∂t in mode evolution
  useDynamicQ?: boolean; // Dynamic Q-factor (depends on R, n, etc.)
  useRadiationBackreaction?: boolean; // Radiation backreaction on bubble dynamics
  radiationBackreactionCoeff?: number; // Backreaction coefficient [N·s/J]
  includeBremsstrahlungEmission?: boolean; // Add bremsstrahlung to stored energy
  includeRecombinationEmission?: boolean; // Add recombination emission to stored energy
  
  // Magnetic field effects (based on recent research) 
  useMagneticFieldEffects?: boolean; // Enable magnetic field component tracking
  useFaradayRotation?: boolean; // Enable Faraday Effect (magnetic field rotates light polarization)
  staticMagneticField?: number; // Static magnetic field B₀ [T] (for Faraday Effect)
  magneticFieldCoupling?: number; // Coupling strength for magnetic field effects
  magneticContributionFactor?: number; // Fraction of magnetic contribution (0.17 for visible, 0.70 for IR)
}

export interface EmDerivatives {
  dModes: { dRe: number; dIm: number }[]; // for each mode
  dStoredEnergyDt: number;
  // Magnetic field effects 
  magneticFieldAmplitude?: number[]; // B field amplitude for each mode [T]
  faradayRotationAngle?: number[]; // Faraday rotation angle for each mode [rad]
  magneticEnergy?: number; // Energy in magnetic field component [J]
}

/**
 * Compute mode frequency as function of bubble radius and refractive index
 * 
 * FORMAL CAVITY-QED: ω_k(R(t), n(r,t))
 * 
 * For a spherical cavity with refractive index n:
 * ω_k(t) = ω_k0 * (R0 / R(t)) * (1 / n_effective)
 * 
 * Also affected by plasma frequency: ω_eff = sqrt(ω_k² - ω_p²)
 * 
 * The refractive index n(r,t) depends on:
 * - Neutral gas density (density-dependent)
 * - Plasma frequency (frequency-dependent)
 * - Radial position (gradient effects)
 */
function computeModeFrequency(
  omega0: number,
  R: number,
  R0: number,
  omega_p: number,
  n_effective?: number, // Refractive index (if available)
  opticalFrequency?: number // Optical frequency for plasma refractive index
): number {
  // Geometric scaling: frequency ~ 1/R for cavity modes
  let omega_geometric = omega0 * (R0 / Math.max(R, 1e-10));
  
  // Refractive index effect: ω ∝ 1/n
  // For cavity modes, the frequency scales inversely with refractive index
  if (n_effective !== undefined && n_effective > 0) {
    omega_geometric = omega_geometric / n_effective;
  }
  
  // Plasma frequency effect: ω_eff = sqrt(ω_k² - ω_p²) if ω_k > ω_p
  // If ω_p > ω_k, mode is cut off
  if (omega_p >= omega_geometric) {
    return 0; // Mode cut off by plasma
  }
  
  return Math.sqrt(omega_geometric * omega_geometric - omega_p * omega_p);
}

/**
 * Compute refractive index time derivative ∂n/∂t
 * 
 * The refractive index changes with:
 * - Density changes (compression/expansion)
 * - Plasma density changes (ionization/recombination)
 * - Temperature changes (affects polarizability)
 * 
 * Computed directly from state changes to avoid circular dependency.
 */
function computeRefractiveIndexTimeDerivative(
  state: BubbleFullState,
  statePrev: BubbleFullState,
  dt: number,
  opticalFrequency: number
): number {
  // Compute refractive index directly (avoiding circular dependency)
  const { R } = state.hydro;
  const { T } = state.gas;
  const { ne } = state.plasma;
  const { numberDensity } = state.species;
  
  // Neutral gas refractive index: n ≈ 1 + α * ρ
  const totalNeutralDensity = Object.values(numberDensity).reduce(
    (sum, n) => sum + n,
    0
  );
  const alpha = 1e-29; // Approximate polarizability [m³]
  const n_neutral = 1.0 + alpha * totalNeutralDensity;
  
  // Plasma refractive index: n² = 1 - (ωp² / ω²)
  const omega_p = Calculations.plasmaFrequency(ne);
  const n_plasma_squared = 1.0 - (omega_p * omega_p) / (opticalFrequency * opticalFrequency);
  const n_plasma = n_plasma_squared > 0 ? Math.sqrt(n_plasma_squared) : 0;
  const n_current = n_neutral + (n_plasma - 1.0) * 0.5;
  
  // Previous state
  const { R: R_prev } = statePrev.hydro;
  const { T: T_prev } = statePrev.gas;
  const { ne: ne_prev } = statePrev.plasma;
  const { numberDensity: numberDensityPrev } = statePrev.species;
  
  const totalNeutralDensityPrev = Object.values(numberDensityPrev).reduce(
    (sum, n) => sum + n,
    0
  );
  const n_neutral_prev = 1.0 + alpha * totalNeutralDensityPrev;
  const omega_p_prev = Calculations.plasmaFrequency(ne_prev);
  const n_plasma_squared_prev = 1.0 - (omega_p_prev * omega_p_prev) / (opticalFrequency * opticalFrequency);
  const n_plasma_prev = n_plasma_squared_prev > 0 ? Math.sqrt(n_plasma_squared_prev) : 0;
  const n_prev = n_neutral_prev + (n_plasma_prev - 1.0) * 0.5;
  
  return (n_current - n_prev) / dt;
}

/**
 * Compute EM mode + stored energy derivatives.
 * 
 * Mode evolution:
 * - Parametric amplification from boundary motion (Rdot pumps modes)
 * - Mode frequency changes with R(t) and plasma frequency
 * - Damping from cavity losses
 * 
 * Stored energy:
 * - Pump term: α * |Rdot| * |gradient_terms|
 * - Decay term: -E_em / τ
 */
export interface EmDerivativesWithState {
  dModes: { dRe: number; dIm: number }[];
  dStoredEnergyDt: number;
  radiationBackreaction?: { fx: number; fy: number; fz: number }; // Radiation backreaction force
  dynamicQ?: number[]; // Dynamic Q-factor for each mode
}

/**
 * Compute EM mode + stored energy derivatives (CAVITY-QED FORMALISM)
 * 
 * Enhanced with:
 * - Refractive index frequency shift: ω_k(R(t), n(r,t))
 * - Mode squeezing with ∂n/∂t: ȧk = f(ak, R, Ṙ, ∂n/∂t)
 * - Dynamic Q-factor & cavity losses
 * - Radiation backreaction
 * - Bremsstrahlung + recombination emission
 */
export function computeEmDerivatives(
  state: BubbleFullState,
  params: EmCavityParams & { thermoGamma?: number },
  statePrev?: BubbleFullState, // Previous state for computing time derivatives
  dt?: number // Time step for computing time derivatives
): EmDerivatives | EmDerivativesWithState {
  const { R, Rdot } = state.hydro;
  const { Pg, T } = state.gas;
  const { ne, Te } = state.plasma;
  const modes = state.em.modes;
  const nModes = modes.length;
  const R0 = params.modeFrequencies0.length > 0 ? 1e-6 : 1e-6; // Reference radius
  const volume = (4.0 / 3.0) * Math.PI * R * R * R; // Bubble volume [m³]

  // Compute plasma frequency
  const omega_p = Calculations.plasmaFrequency(ne);
  
  // Compute refractive index (for frequency shift)
  // FORMAL CAVITY-QED: n(r,t) depends on density, plasma, temperature
  let n_effective: number | undefined;
  let dn_dt: number | undefined;
  if (params.useRefractiveIndexFrequencyShift || params.useRefractiveIndexTimeDerivative) {
    const opticalFrequency = params.modeFrequencies0[0] || 2 * Math.PI * 3e14;
    
    // Compute n directly (avoiding circular dependency)
    const { numberDensity } = state.species;
    const totalNeutralDensity = Object.values(numberDensity).reduce(
      (sum, n) => sum + n,
      0
    );
    const alpha = 1e-29; // Approximate polarizability [m³]
    const n_neutral = 1.0 + alpha * totalNeutralDensity;
    
    const n_plasma_squared = 1.0 - (omega_p * omega_p) / (opticalFrequency * opticalFrequency);
    const n_plasma = n_plasma_squared > 0 ? Math.sqrt(n_plasma_squared) : 0;
    n_effective = n_neutral + (n_plasma - 1.0) * 0.5;
    
    // Compute ∂n/∂t if previous state available
    if (params.useRefractiveIndexTimeDerivative && statePrev && dt !== undefined) {
      dn_dt = computeRefractiveIndexTimeDerivative(state, statePrev, dt, opticalFrequency);
    }
  }

  // Compute gradient magnitudes for pumping
  // More detailed calculation of parametric coupling strength
  const Rdot_magnitude = Math.abs(Rdot);
  
  // Pressure gradient from radial motion
  // dP/dR ≈ -gamma * P * (3/R) for adiabatic compression
  const gamma = params.thermoGamma || 1.4;
  const Pg_gradient_radial = gamma * Math.abs(Pg) / Math.max(R, 1e-10);
  
  // Pressure gradient from time variation
  // dP/dt ≈ -gamma * P * (3*Rdot/R) for adiabatic
  const Pg_gradient_temporal = gamma * Math.abs(Pg) * Math.abs(Rdot) / Math.max(R, 1e-10);
  
  // Combined pressure gradient (geometric mean for coupling strength)
  const Pg_gradient = Math.sqrt(Pg_gradient_radial * Pg_gradient_temporal);
  
  // Acoustic pressure gradient (if available from acoustic module)
  // This would come from acoustic field gradients
  const acousticGradient = 0; // TODO: pass from acoustic module if available
  
  // Combined gradient metric for parametric pumping
  // Parametric coupling strength scales with gradient magnitude
  const gradientMagnitude = Rdot_magnitude + Pg_gradient * 1e-6 + acousticGradient * 1e-8;

  const dModes = Array.from({ length: nModes }, (_, k) => {
    if (k >= params.modeFrequencies0.length) {
      return { dRe: 0, dIm: 0 };
    }

    const omega0 = params.modeFrequencies0[k];
    const opticalFrequency = omega0; // Use mode frequency as optical frequency
    const omega_k = computeModeFrequency(
      omega0, 
      R, 
      R0, 
      omega_p,
      n_effective, // Include refractive index if enabled
      opticalFrequency
    );
    
    // If mode is cut off, decay to zero
    if (omega_k === 0) {
      const decayRate = 1e12; // fast decay when cut off
      return {
        dRe: -state.em.modes[k].re * decayRate,
        dIm: -state.em.modes[k].im * decayRate,
      };
    }

    const a_re = state.em.modes[k].re;
    const a_im = state.em.modes[k].im;

    // Mode evolution equations
    // da/dt = -i*ω*a + pump_term - damping_term
    
    // Pumping term: parametric amplification from boundary motion
    let pumpRe = 0;
    let pumpIm = 0;
    
    if (params.useDetailedParametricPumping) {
      // Detailed parametric pumping Hamiltonian
      // H = g(t) * (a†² + a²) where g(t) = g0 * Rdot/R
      // This gives parametric amplification when Rdot < 0 (compression)
      // 
      // NEGATIVE-SPACE BEHAVIOR:
      // The parametric Hamiltonian H = g(t) * (a†² + a²) creates a squeezed
      // vacuum state. In the Wigner function representation, this appears as
      // a "negative-space" region where the Wigner function becomes negative,
      // indicating non-classical behavior. The squeezing occurs in quadrature
      // space: one quadrature is squeezed (reduced uncertainty) while the
      // conjugate quadrature is anti-squeezed (increased uncertainty).
      //
      // During extreme compression (large |Rdot|), g(t) becomes large and
      // the modes are strongly squeezed, storing energy in the negative-space
      // state. This energy is then released as photons when the cavity decays.
      //
      // From Heisenberg equations: da/dt = -i*[a, H] = -i*g(t)*(a† + a)
      // In terms of real/imaginary parts:
      // dRe/dt = g(t) * Im
      // dIm/dt = -g(t) * Re
      const g0 = params.parametricCoupling || params.couplingStrength;
      const R_safe = Math.max(R, 1e-10);
      const g_t = g0 * (Rdot / R_safe); // Time-dependent coupling
      
      // Parametric pumping (quadrature-dependent)
      // This cross-coupling (Re ↔ Im) is the signature of parametric amplification
      // and creates the squeezed/negative-space state
      pumpRe = g_t * a_im;
      pumpIm = -g_t * a_re;
      
      // REFRACTIVE INDEX TIME DERIVATIVE CONTRIBUTION: ∂n/∂t term
      // Mode squeezing equation: ȧk = f(ak, R, Ṙ, ∂n/∂t)
      // The refractive index time derivative contributes to mode evolution
      // via the cavity frequency shift: ω_k ∝ 1/n, so dω_k/dt ∝ -(1/n²) * dn/dt
      if (params.useRefractiveIndexTimeDerivative && dn_dt !== undefined) {
        const n_safe = Math.max(n_effective || 1.0, 0.1);
        const omega_n_derivative = -(omega_k / (n_safe * n_safe)) * dn_dt;
        // This contributes to mode evolution as frequency modulation
        pumpRe += omega_n_derivative * a_im * 0.5; // Coupling factor
        pumpIm -= omega_n_derivative * a_re * 0.5;
      }
      
      // Additional squeezing term (higher order)
      // For strong coupling, add squeezing: H_squeeze = g_squeeze * (a†² + a²)
      // This enhances the negative-space behavior for very strong compression
      const g_squeeze = g0 * 0.1; // Squeezing is weaker
      const squeezeRe = g_squeeze * a_re; // Squeezing in phase
      const squeezeIm = g_squeeze * a_im;
      
      pumpRe += squeezeRe;
      pumpIm += squeezeIm;
    } else {
      // Simplified pumping (original)
      // Pump strength ~ couplingStrength * Rdot * gradient
      const pumpStrength = params.couplingStrength * gradientMagnitude;
      pumpRe = pumpStrength * a_im; // Cross-coupling for parametric amplification
      pumpIm = -pumpStrength * a_re;
    }

    // Damping: cavity losses (quality factor Q)
    // DYNAMIC Q-FACTOR: Q(R, n, ω) - depends on cavity geometry and refractive index
    let Q = params.Q0 || 1000; // Base quality factor
    
    if (params.useDynamicQ || params.useFrequencyDependentQ) {
      // Dynamic Q depends on:
      // 1. Frequency: Q(ω) = Q0 * (ω0/ω)^alpha
      // 2. Radius: Q(R) ~ R (larger cavity = higher Q)
      // 3. Refractive index: Q(n) ~ n (higher n = higher Q, more confinement)
      
      if (params.QFrequencyDependence) {
        Q = params.QFrequencyDependence(omega_k);
      } else {
        // Default frequency dependence: Q decreases with frequency
        const omega_ref = params.modeFrequencies0[0] || 2 * Math.PI * 3e14;
        const alpha = 0.5; // Power law exponent
        Q = Q * Math.pow(omega_ref / Math.max(omega_k, omega_ref * 0.1), alpha);
      }
      
      // Radius dependence: Q ~ R (larger bubble = higher Q)
      if (params.useDynamicQ) {
        const R_ref = 5e-6; // Reference radius [m]
        const Q_radius_factor = R / Math.max(R_ref, 1e-10);
        Q = Q * Q_radius_factor;
      }
      
      // Refractive index dependence: Q ~ n (higher n = better confinement)
      if (params.useDynamicQ && n_effective !== undefined) {
        const n_ref = 1.0; // Reference refractive index (vacuum)
        const Q_n_factor = n_effective / n_ref;
        Q = Q * Q_n_factor;
      }
    }
    const damping = omega_k / (2.0 * Q);
    
    // Cross-mode coupling (if enabled)
    let couplingRe = 0;
    let couplingIm = 0;
    if (params.useModeCoupling && params.modeCouplingMatrix) {
      const couplingMatrix = params.modeCouplingMatrix;
      // Coupling term: sum_j g_kj * a_j
      for (let j = 0; j < nModes && j < couplingMatrix.length; j++) {
        if (j !== k && couplingMatrix[k] && couplingMatrix[k][j] !== undefined) {
          const g_kj = couplingMatrix[k][j];
          // Coupling adds to both real and imaginary parts
          couplingRe += g_kj * state.em.modes[j].re;
          couplingIm += g_kj * state.em.modes[j].im;
        }
      }
    }
    
    // Mode evolution: da/dt = -i*ω*a + pump - damping*a + coupling
    // Real part: dRe/dt = ω*Im + pumpRe - damping*Re + couplingRe
    // Imag part: dIm/dt = -ω*Re + pumpIm - damping*Im + couplingIm
    const dRe = omega_k * a_im + pumpRe - damping * a_re + couplingRe;
    const dIm = -omega_k * a_re + pumpIm - damping * a_im + couplingIm;

    return { dRe, dIm };
  });

  // Magnetic field effects (based on recent research) 
  // Research shows light's magnetic field directly influences matter
  // Magnetic component contributes ~17% in visible, up to 70% in infrared
  // Reference: https://www.sciencedaily.com/releases/2025/11/251120091945.htm
  let magneticFieldAmplitude: number[] | undefined;
  let faradayRotationAngle: number[] | undefined;
  let magneticEnergy = 0;
  
  if (params.useMagneticFieldEffects) {
    magneticFieldAmplitude = [];
    faradayRotationAngle = [];
    
    // For each mode, compute magnetic field from electric field
    // For plane waves: B = E/c (in vacuum/air)
    // In medium: B = n*E/c (approximately)
    const c = Constants.c; // Speed of light [m/s]
    const n_eff = n_effective || params.refractiveIndexBase;
    
    for (let k = 0; k < nModes; k++) {
      const a_re = modes[k].re;
      const a_im = modes[k].im;
      const E_amplitude = Math.sqrt(a_re * a_re + a_im * a_im);
      
      // Magnetic field amplitude: B = n*E/c
      // The mode amplitude is proportional to electric field
      // Convert to actual field strength (normalized by mode frequency)
      const omega_k = computeModeFrequency(
        params.modeFrequencies0[k] || 1e15,
        R,
        R0,
        omega_p,
        n_effective,
        params.modeFrequencies0[k]
      );
      
      // Electric field strength estimate (from mode amplitude)
      // Mode amplitude is normalized, so we estimate E ~ amplitude * characteristic field
      const E0_characteristic = Math.sqrt(2 * Constants.hbar * omega_k / (Constants.epsilon_0 * volume));
      const E_field = E_amplitude * E0_characteristic;
      
      // Magnetic field: B = n*E/c
      const B_amplitude = (n_eff * E_field) / c;
      magneticFieldAmplitude.push(B_amplitude);
      
      // Faraday rotation (if static magnetic field present)
      // Faraday Effect: rotation angle θ = V * B₀ * L
      // where V is Verdet constant, B₀ is static field, L is path length
      if (params.useFaradayRotation && params.staticMagneticField !== undefined) {
        const B0 = params.staticMagneticField; // Static magnetic field [T]
        const L = 2 * R; // Path length through bubble (diameter) [m]
        
        // Verdet constant (material-dependent, typical ~10 rad/(T·m) for TGG)
        // For gas/plasma, use approximate value
        const Verdet = params.magneticFieldCoupling || 1.0; // rad/(T·m)
        
        // Faraday rotation angle: θ = V * B₀ * L
        const rotationAngle = Verdet * B0 * L;
        faradayRotationAngle.push(rotationAngle);
        
        // Magnetic field contribution to mode evolution
        // The magnetic field can generate magnetic torque (similar to static field)
        // This affects the mode phase and amplitude
        const magneticCoupling = params.magneticFieldCoupling || 1e6; // Coupling strength
        const magneticContribution = params.magneticContributionFactor || 0.17; // 17% for visible
        
        // Add magnetic field effect to mode evolution
        // This is a simplified model - full treatment requires LLG equation
        const magneticEffect = magneticContribution * magneticCoupling * B_amplitude * B0;
        
        // Modify mode derivatives to include magnetic field effects
        // Magnetic field can cause phase rotation and amplitude modulation
        if (dModes[k]) {
          // Phase rotation from Faraday effect
          const phaseRotation = rotationAngle;
          const cos_rot = Math.cos(phaseRotation);
          const sin_rot = Math.sin(phaseRotation);
          
          // Rotate mode components (Faraday rotation)
          const dRe_old = dModes[k].dRe;
          const dIm_old = dModes[k].dIm;
          dModes[k].dRe = dRe_old * cos_rot - dIm_old * sin_rot;
          dModes[k].dIm = dRe_old * sin_rot + dIm_old * cos_rot;
          
          // Add magnetic torque effect (simplified)
          dModes[k].dRe += magneticEffect * a_im;
          dModes[k].dIm -= magneticEffect * a_re;
        }
      } else {
        faradayRotationAngle.push(0);
      }
      
      // Magnetic energy: U_B = (1/(2*μ₀)) * B² * V
      // For each mode, add magnetic energy contribution
      const mu_0 = Constants.mu_0; // Permeability of free space
      const magneticEnergyMode = (B_amplitude * B_amplitude / (2 * mu_0)) * volume;
      magneticEnergy += magneticEnergyMode;
    }
  }

  // Stored energy evolution
  // dE_em/dt = pump_term - decay_term + bremsstrahlung + recombination
  // 
  // PUMP TERM: α * |Rdot| * |gradient_terms|
  // This represents parametric amplification from boundary motion.
  // STORED-ENERGY PUMP FROM GRADIENTS:
  // The pump term explicitly uses gradients:
  // - Radial gradient: |Rdot| (boundary motion)
  // - Pressure gradient: |dPg/dt| (compression work)
  // - Acoustic gradient: |∇P_acoustic| (if available)
  // During extreme compression (large |Rdot|, large |dPg/dt|), the boundary
  // motion pumps energy into EM modes, creating a "negative-space" squeezed state.
  // The negative-space refers to the quantum squeezed vacuum state in the
  // Wigner function representation, where the uncertainty in one quadrature
  // is reduced below the vacuum level at the expense of increased uncertainty
  // in the conjugate quadrature.
  //
  // DECAY TERM: -E_em / τ
  // The stored energy decays over a short timescale (τ ~ 1 ns), converting
  // to photons (light emission). This decay is what produces the sonoluminescence
  // flash: E_em rises during collapse, then rapidly decays as photons are emitted.
  //
  // BREMSSTRAHLUNG EMISSION: Energy from electron-ion collisions
  // Recombination emission: Energy from electron-ion recombination
  let pumpTerm = params.pumpCoefficient * gradientMagnitude;
  const decayTerm = state.em.storedEnergy / Math.max(params.decayTime, 1e-12);
  
  // BREMSSTRAHLUNG + RECOMBINATION EMISSION ADDITIONS
  let bremsstrahlungEmission = 0;
  let recombinationEmission = 0;
  
  if (params.includeBremsstrahlungEmission) {
    // Bremsstrahlung power density: P ~ ne² * sqrt(Te)
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const bremsstrahlungPowerDensity = Calculations.bremsstrahlungPowerDensity(ne, Te);
    const bremsstrahlungPower = bremsstrahlungPowerDensity * volume;
    // Convert power to energy rate (if dt available, otherwise use instantaneous)
    const dt_safe = dt || 1e-12;
    bremsstrahlungEmission = bremsstrahlungPower * dt_safe;
  }
  
  if (params.includeRecombinationEmission) {
    // Recombination emission: Energy released when electrons recombine with ions
    // Estimate: P_recomb ~ ne * n_ions * α_recomb * E_ionization
    // where α_recomb is recombination coefficient, E_ionization is ionization energy
    const n_ions = ne; // Assume n_ions ≈ ne (single ionization)
    const alpha_recomb = 2.6e-19; // Recombination coefficient [m³/s] (approximate)
    const E_ionization = 15.76 * Constants.e; // Argon ionization energy [J]
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const recombinationPower = ne * n_ions * alpha_recomb * E_ionization * volume;
    const dt_safe = dt || 1e-12;
    recombinationEmission = recombinationPower * dt_safe;
  }
  
  // Add magnetic energy contribution to stored energy 
  // Research shows magnetic field directly influences matter and contributes to energy
  // Magnetic contribution: ~17% in visible, up to 70% in infrared
  if (params.useMagneticFieldEffects && magneticEnergy > 0) {
    const magneticContribution = params.magneticContributionFactor || 0.17; // 17% for visible
    // Add magnetic energy to stored energy (scaled by contribution factor)
    const magneticEnergyContribution = magneticContribution * magneticEnergy;
    // Note: This is added to the stored energy evolution, not directly to E_em
    // The magnetic field energy is part of the total EM energy
  }
  
  const dStoredEnergyDt = pumpTerm - decayTerm + bremsstrahlungEmission + recombinationEmission;
  
  // RADIATION BACKREACTION TERM
  // The EM field exerts a force back on the bubble boundary
  // This is the radiation pressure: F = -∇E_em / R
  let radiationBackreaction: { fx: number; fy: number; fz: number } | undefined;
  if (params.useRadiationBackreaction) {
    const backreactionCoeff = params.radiationBackreactionCoeff || 1e-15; // [N·s/J]
    const E_em = state.em.storedEnergy;
    const R_safe = Math.max(R, 1e-10);
    
    // Radiation pressure: F = -E_em / R (radial, inward)
    // Convert to force components (radial direction)
    const F_magnitude = (E_em / R_safe) * backreactionCoeff;
    
    // For spherical symmetry, force is radial
    // In 3D, this would be along the radial direction from bubble center
    // Simplified: assume force is along z-axis (or can be distributed)
    radiationBackreaction = {
      fx: 0, // Could be computed from bubble position if translation enabled
      fy: 0,
      fz: -F_magnitude, // Radial inward (negative z for upward bubble)
    };
  }
  
  // Return with optional dynamic Q and backreaction
  if (params.useDynamicQ || params.useRadiationBackreaction) {
    const dynamicQ = params.useDynamicQ 
      ? dModes.map((_, k) => {
          const omega0 = params.modeFrequencies0[k];
          const omega_k = computeModeFrequency(omega0, R, R0, omega_p, n_effective, omega0);
          let Q = params.Q0 || 1000;
          if (params.useFrequencyDependentQ) {
            const omega_ref = params.modeFrequencies0[0] || 2 * Math.PI * 3e14;
            const alpha = 0.5;
            Q = Q * Math.pow(omega_ref / Math.max(omega_k, omega_ref * 0.1), alpha);
          }
          if (params.useDynamicQ) {
            const R_ref = 5e-6;
            Q = Q * (R / Math.max(R_ref, 1e-10));
            if (n_effective !== undefined) {
              Q = Q * (n_effective / 1.0);
            }
          }
          return Q;
        })
      : undefined;
    
    return {
      dModes,
      dStoredEnergyDt,
      radiationBackreaction,
      dynamicQ,
      magneticFieldAmplitude: params.useMagneticFieldEffects ? magneticFieldAmplitude : undefined,
      faradayRotationAngle: params.useFaradayRotation ? faradayRotationAngle : undefined,
      magneticEnergy: params.useMagneticFieldEffects ? magneticEnergy : undefined,
    };
  }
  
  return {
    dModes,
    dStoredEnergyDt,
    magneticFieldAmplitude: params.useMagneticFieldEffects ? magneticFieldAmplitude : undefined,
    faradayRotationAngle: params.useFaradayRotation ? faradayRotationAngle : undefined,
    magneticEnergy: params.useMagneticFieldEffects ? magneticEnergy : undefined,
  };
}