// physics/bubbleTranslation.ts
// Bubble translation dynamics and Bjerknes forces

import { BubbleFullState, BubbleTranslationState } from "../model/types";
import { Constants } from "../core/units";

export interface BubbleTranslationParams {
  // Liquid density [kg/m³]
  rho: number;
  
  // Viscosity [Pa·s]
  mu: number;
  
  // Enable translation
  enableTranslation?: boolean;
  
  // Drag coefficient (Stokes drag: 6π for sphere)
  dragCoeff?: number;
  
  // Acoustic field parameters
  acousticGradient?: { x: number; y: number; z: number };  // ∇P_acoustic [Pa/m]
  acousticLaplacian?: number;  // ∇²P_acoustic [Pa/m²]
  
  // Buoyancy force
  enableBuoyancy?: boolean; // Enable buoyancy force calculation
  gravity?: number; // Gravity acceleration [m/s²] (default: 9.81)
  
  // Secondary Bjerknes force
  enableSecondaryBjerknes?: boolean; // Enable secondary Bjerknes force
  otherBubbles?: Array<{
    position: { x: number; y: number; z: number };
    volume: number;
  }>; // Other bubbles for interaction
}

export interface BubbleTranslationDerivatives {
  dx_dt: number;
  dy_dt: number;
  dz_dt: number;
  dvx_dt: number;
  dvy_dt: number;
  dvz_dt: number;
}

/**
 * Compute primary Bjerknes force
 * 
 * Primary Bjerknes force: F = -V * ∇P_acoustic
 * 
 * This force drives bubbles toward pressure nodes (low pressure regions)
 * or antinodes (high pressure regions) depending on bubble size.
 * 
 * For small bubbles: move to pressure antinodes
 * For large bubbles: move to pressure nodes
 */
function computePrimaryBjerknesForce(
  V: number,  // Bubble volume [m³]
  gradP: { x: number; y: number; z: number }  // Pressure gradient [Pa/m]
): { fx: number; fy: number; fz: number } {
  return {
    fx: -V * gradP.x,
    fy: -V * gradP.y,
    fz: -V * gradP.z,
  };
}

/**
 * Compute drag force (Stokes drag for small Reynolds number)
 * 
 * F_drag = -6π*μ*R*v
 * 
 * For larger Reynolds numbers, use:
 * F_drag = -0.5*C_d*ρ*A*v²
 */
function computeDragForce(
  R: number,      // Bubble radius [m]
  v: { x: number; y: number; z: number },  // Velocity [m/s]
  mu: number,     // Viscosity [Pa·s]
  rho: number,    // Density [kg/m³]
  dragCoeff?: number
): { fx: number; fy: number; fz: number } {
  const v_mag = Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  
  if (dragCoeff !== undefined) {
    // Custom drag coefficient
    const A = Math.PI * R * R;  // Cross-sectional area
    const F_mag = 0.5 * dragCoeff * rho * A * v_mag * v_mag;
    if (v_mag > 0) {
      return {
        fx: -F_mag * (v.x / v_mag),
        fy: -F_mag * (v.y / v_mag),
        fz: -F_mag * (v.z / v_mag),
      };
    }
  } else {
    // Stokes drag (low Reynolds number)
    const F_mag = 6.0 * Math.PI * mu * R * v_mag;
    if (v_mag > 0) {
      return {
        fx: -F_mag * (v.x / v_mag),
        fy: -F_mag * (v.y / v_mag),
        fz: -F_mag * (v.z / v_mag),
      };
    }
  }
  
  return { fx: 0, fy: 0, fz: 0 };
}

/**
 * Compute bubble translation derivatives
 * 
 * Position: dx/dt = vx, dy/dt = vy, dz/dt = vz
 * 
 * Velocity: dv/dt = F_total / m
 * 
 * Forces:
 * 1. Primary Bjerknes force: F = -V * ∇P_acoustic
 * 2. Drag force: F_drag = -6π*μ*R*v (Stokes)
 * 3. Buoyancy (if gravity included)
 */
export function computeBubbleTranslationDerivatives(
  state: BubbleFullState,
  params: BubbleTranslationParams
): BubbleTranslationDerivatives {
  if (!params.enableTranslation || !state.translation) {
    return {
      dx_dt: 0, dy_dt: 0, dz_dt: 0,
      dvx_dt: 0, dvy_dt: 0, dvz_dt: 0,
    };
  }
  
  const { R } = state.hydro;
  const { x, y, z, vx, vy, vz } = state.translation;
  const { rho, mu, dragCoeff } = params;
  
  // Position derivatives (kinematics)
  const dx_dt = vx;
  const dy_dt = vy;
  const dz_dt = vz;
  
  // Bubble volume
  const V = (4.0 / 3.0) * Math.PI * R * R * R;
  
  // Bubble mass (gas + liquid equivalent)
  // Simplified: use liquid density for added mass
  const m = rho * V;  // Effective mass (includes added mass)
  
  // Forces
  let fx_total = 0;
  let fy_total = 0;
  let fz_total = 0;
  
  // 1. Primary Bjerknes force (if acoustic gradient provided)
  if (params.acousticGradient) {
    const F_bjerknes = computePrimaryBjerknesForce(V, params.acousticGradient);
    fx_total += F_bjerknes.fx;
    fy_total += F_bjerknes.fy;
    fz_total += F_bjerknes.fz;
  }
  
  // 2. Drag force
  const F_drag = computeDragForce(R, { x: vx, y: vy, z: vz }, mu, rho, dragCoeff);
  fx_total += F_drag.fx;
  fy_total += F_drag.fy;
  fz_total += F_drag.fz;
  
  // 3. Buoyancy force (if enabled)
  if (params.enableBuoyancy) {
    const g = params.gravity || 9.81; // Gravity [m/s²]
    
    // Compute gas density from ideal gas law: ρ_gas = P * M_avg / (R * T)
    const { Pg, T } = state.gas;
    const { numberDensity } = state.species;
    
    // Compute average molecular weight
    // M_avg = Σ(n_i * M_i) / Σ(n_i)
    const speciesMolarMasses: Record<string, number> = {
      H2O: 0.018015,  // kg/mol
      O2: 0.031999,
      N2: 0.028014,
      Ar: 0.039948,
      Xe: 0.131293,
      H: 0.001008,
      O: 0.015999,
      OH: 0.017007,
      N: 0.014007,
    };
    
    let totalMoles = 0;
    let totalMass = 0;
    
    for (const [species, density] of Object.entries(numberDensity)) {
      const n_moles = density * V; // Number of moles
      const M = speciesMolarMasses[species] || 0.03; // Default ~air
      totalMoles += n_moles;
      totalMass += n_moles * M;
    }
    
    // Gas density: ρ_gas = total_mass / V = (Σ n_i * M_i) / V
    const rho_gas = totalMoles > 0 ? totalMass / V : 0;
    
    // Alternative: use ideal gas law if we have pressure and temperature
    // ρ_gas = P * M_avg / (R * T)
    const M_avg = totalMoles > 0 ? totalMass / totalMoles : 0.03;
    const rho_gas_ideal = (Pg * M_avg) / (Constants.R_gas * T);
    
    // Use ideal gas law result (more accurate for high pressure)
    const rho_gas_final = rho_gas_ideal;
    
    // Buoyancy force: F_buoyancy = (ρ_liquid - ρ_gas) * V * g
    const F_buoyancy = (rho - rho_gas_final) * V * g;
    
    // Buoyancy acts upward (positive z direction)
    fz_total += F_buoyancy;
  }
  
  // 4. Secondary Bjerknes force (bubble-bubble interactions)
  if (params.enableSecondaryBjerknes && params.otherBubbles && params.acousticLaplacian !== undefined) {
    const laplacianP = params.acousticLaplacian;
    
    for (const otherBubble of params.otherBubbles) {
      // Distance vector from this bubble to other bubble
      const dx = otherBubble.position.x - x;
      const dy = otherBubble.position.y - y;
      const dz = otherBubble.position.z - z;
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      
      if (r > 1e-10) { // Avoid division by zero
        // Secondary Bjerknes force: F = -ρ * V₁ * V₂ * ∇²P / (4π * r²)
        // Direction: along line connecting bubbles
        const F_secondary = computeSecondaryBjerknesForce(
          V,
          otherBubble.volume,
          r,
          laplacianP,
          rho
        );
        
        // Force components
        const F_mag = Math.abs(F_secondary);
        fx_total += F_mag * (dx / r);
        fy_total += F_mag * (dy / r);
        fz_total += F_mag * (dz / r);
      }
    }
  }
  
  // Velocity derivatives (dynamics)
  const dvx_dt = fx_total / m;
  const dvy_dt = fy_total / m;
  const dvz_dt = fz_total / m;
  
  return { dx_dt, dy_dt, dz_dt, dvx_dt, dvy_dt, dvz_dt };
}

/**
 * Compute secondary Bjerknes force between two bubbles
 * 
 * Secondary Bjerknes force: F = -ρ * V₁ * V₂ * ∇²P / (4π * r²)
 * 
 * This is the interaction force between bubbles.
 * For now, we implement the force on a single bubble from the field.
 */
export function computeSecondaryBjerknesForce(
  V: number,  // Bubble volume [m³]
  V_other: number,  // Other bubble volume [m³]
  r: number,  // Distance between bubbles [m]
  laplacianP: number,  // ∇²P_acoustic [Pa/m²]
  rho: number  // Liquid density [kg/m³]
): number {
  const r_safe = Math.max(r, 1e-10);
  return -(rho * V * V_other * laplacianP) / (4.0 * Math.PI * r_safe * r_safe);
}

