// physics/acoustic.ts
import { AcousticState, BubbleFullState } from "../model/types";

export interface AcousticParams {
  // Single frequency mode (backward compatible)
  Pa?: number;   // amplitude [Pa]
  omega?: number; // frequency [rad/s]
  
  // Multi-frequency mode
  frequencies?: Array<{
    amplitude: number; // [Pa]
    frequency: number; // [rad/s]
    phase: number;     // initial phase [rad]
  }>;
  
  // Non-sinusoidal driving
  waveform?: 'sinusoidal' | 'square' | 'sawtooth' | 'custom';
  customWaveform?: (phase: number) => number; // Custom waveform function
  
  // Acoustic field gradients (for Bjerknes forces)
  enableGradients?: boolean; // Compute spatial gradients
  waveVector?: { x: number; y: number; z: number }; // Wave vector k [1/m]
  standingWave?: boolean; // Standing wave pattern
  nodePosition?: { x: number; y: number; z: number }; // Node position for standing wave
  
  // Jet-driven cavitation (pistol shrimp mechanism) 
  jetDriven?: boolean; // Enable jet-driven cavitation
  jetVelocity?: number; // Water jet velocity [m/s] (typical: 30-60 m/s for pistol shrimp)
  jetDuration?: number; // Jet duration [s] (typical: ~0.1-1 ms)
  jetNozzleArea?: number; // Effective nozzle area [m²] (typical: ~1 mm²)
  jetDirection?: { x: number; y: number; z: number }; // Jet direction vector (normalized)
}

export interface AcousticDerivatives {
  dPhaseDt: number;
  Pacoustic: number; // value of acoustic pressure at this instant
  gradient?: { x: number; y: number; z: number }; // ∇P_acoustic [Pa/m]
  laplacian?: number; // ∇²P_acoustic [Pa/m²]
  // Jet-driven cavitation 
  jetPressure?: number; // Pressure from water jet [Pa]
  jetActive?: boolean; // Whether jet is currently active
  shockwavePressure?: number; // Shockwave pressure from bubble collapse [Pa]
  shockwaveIntensity?: number; // Sound intensity in dB (reference: 1e-12 W/m²)
}

/**
 * Compute acoustic state with support for:
 * - Single frequency (backward compatible)
 * - Multi-frequency driving
 * - Non-sinusoidal waveforms
 * - Spatial gradients using actual bubble position
 */
export function computeAcousticState(
  t: number,
  acoustic: AcousticState,
  params: AcousticParams,
  bubblePosition?: { x: number; y: number; z: number } // Optional: actual bubble position
): AcousticDerivatives {
  // Determine primary frequency for phase evolution
  let primaryOmega = params.omega || 0;
  if (params.frequencies && params.frequencies.length > 0) {
    primaryOmega = params.frequencies[0].frequency;
  }
  
  const dPhaseDt = primaryOmega;
  
  // Compute acoustic pressure
  let Pacoustic = 0;
  let jetPressure = 0;
  let jetActive = false;
  let shockwavePressure = 0;
  let shockwaveIntensity = 0;
  
  // Jet-driven cavitation (pistol shrimp mechanism) 
  if (params.jetDriven && params.jetVelocity !== undefined) {
    const jetDuration = params.jetDuration || 1e-4; // Default 0.1 ms
    const rho_water = 998.2; // kg/m³ (water density)
    
    // Jet is active during the jet duration
    jetActive = t < jetDuration;
    
    if (jetActive) {
      // Bernoulli principle: P_jet = (1/2) * ρ * v²
      // This creates a low-pressure zone that initiates cavitation
      // The pressure difference drives bubble formation
      const v_jet = params.jetVelocity;
      jetPressure = -0.5 * rho_water * v_jet * v_jet; // Negative = low pressure (cavitation)
      
      // Add to acoustic pressure (negative pressure creates cavitation)
      Pacoustic += jetPressure;
    }
  }
  
  // Standard acoustic driving (if not jet-driven or in addition to jet)
  if (params.frequencies && params.frequencies.length > 0) {
    // Multi-frequency mode
    for (const freq of params.frequencies) {
      const phase = freq.frequency * t + freq.phase;
      Pacoustic -= freq.amplitude * Math.sin(phase);
    }
  } else if (params.Pa !== undefined && params.omega !== undefined) {
    // Single frequency mode (backward compatible)
    Pacoustic = -params.Pa * Math.sin(acoustic.phase);
  }
  
  // Apply waveform transformation if specified
  if (params.waveform && params.waveform !== 'sinusoidal') {
    const normalizedPhase = (acoustic.phase % (2 * Math.PI)) / (2 * Math.PI);
    
    switch (params.waveform) {
      case 'square':
        Pacoustic = Pacoustic > 0 ? Math.abs(Pacoustic) : -Math.abs(Pacoustic);
        break;
      case 'sawtooth':
        Pacoustic *= 2 * (normalizedPhase - 0.5);
        break;
      case 'custom':
        if (params.customWaveform) {
          Pacoustic = params.customWaveform(acoustic.phase);
        }
        break;
    }
  }
  
  // Compute gradients if enabled
  let gradient: { x: number; y: number; z: number } | undefined;
  let laplacian: number | undefined;
  
  if (params.enableGradients) {
    // For a plane wave: P = P₀ * sin(k·r - ωt)
    // Gradient: ∇P = P₀ * k * cos(k·r - ωt)
    // Laplacian: ∇²P = -P₀ * k² * sin(k·r - ωt)
    
    const k = params.waveVector || { x: 0, y: 0, z: 0 };
    const k_mag = Math.sqrt(k.x * k.x + k.y * k.y + k.z * k.z);
    
    if (k_mag > 0) {
      // Compute phase at bubble position (use actual position if available)
      const bubblePos = bubblePosition || { x: 0, y: 0, z: 0 };
      const k_dot_r = k.x * bubblePos.x + k.y * bubblePos.y + k.z * bubblePos.z;
      const phase_spatial = k_dot_r - primaryOmega * t;
      
      // Gradient
      const cos_phase = Math.cos(phase_spatial);
      if (params.frequencies && params.frequencies.length > 0) {
        const P0 = params.frequencies[0].amplitude;
        gradient = {
          x: P0 * k.x * cos_phase,
          y: P0 * k.y * cos_phase,
          z: P0 * k.z * cos_phase,
        };
        
        // Laplacian
        const sin_phase = Math.sin(phase_spatial);
        laplacian = -P0 * k_mag * k_mag * sin_phase;
      } else if (params.Pa !== undefined) {
        const cos_phase_local = Math.cos(acoustic.phase);
        gradient = {
          x: params.Pa * k.x * cos_phase_local,
          y: params.Pa * k.y * cos_phase_local,
          z: params.Pa * k.z * cos_phase_local,
        };
        laplacian = -params.Pa * k_mag * k_mag * Math.sin(acoustic.phase);
      }
    }
    
    // Standing wave pattern (if specified)
    if (params.standingWave && params.nodePosition) {
      const node = params.nodePosition;
      const bubblePos = bubblePosition || { x: 0, y: 0, z: 0 };
      
      // Standing wave: P = P₀ * sin(k·(r - r_node)) * sin(ωt)
      // where r_node is the node position
      // 
      // Gradient: ∇P = P₀ * k * cos(k·(r - r_node)) * sin(ωt)
      // Laplacian: ∇²P = -P₀ * k² * sin(k·(r - r_node)) * sin(ωt)
      
      // Position relative to node
      const r_rel = {
        x: bubblePos.x - node.x,
        y: bubblePos.y - node.y,
        z: bubblePos.z - node.z,
      };
      
      // Wave vector dot product with relative position
      const k_dot_r_rel = k.x * r_rel.x + k.y * r_rel.y + k.z * r_rel.z;
      
      // Standing wave factors
      const sin_spatial = Math.sin(k_dot_r_rel);
      const cos_spatial = Math.cos(k_dot_r_rel);
      const sin_temporal = Math.sin(primaryOmega * t);
      
      if (gradient) {
        // Proper standing wave gradient
        if (params.frequencies && params.frequencies.length > 0) {
          const P0 = params.frequencies[0].amplitude;
          gradient = {
            x: P0 * k.x * cos_spatial * sin_temporal,
            y: P0 * k.y * cos_spatial * sin_temporal,
            z: P0 * k.z * cos_spatial * sin_temporal,
          };
        } else if (params.Pa !== undefined) {
          gradient = {
            x: params.Pa * k.x * cos_spatial * sin_temporal,
            y: params.Pa * k.y * cos_spatial * sin_temporal,
            z: params.Pa * k.z * cos_spatial * sin_temporal,
          };
        }
        
        // Update laplacian for standing wave
        if (params.frequencies && params.frequencies.length > 0) {
          const P0 = params.frequencies[0].amplitude;
          laplacian = -P0 * k_mag * k_mag * sin_spatial * sin_temporal;
        } else if (params.Pa !== undefined) {
          laplacian = -params.Pa * k_mag * k_mag * sin_spatial * sin_temporal;
        }
      }
    }
  }
  
  // Shockwave from bubble collapse (pistol shrimp) 
  // The shockwave is generated when the bubble collapses rapidly
  // Intensity can reach 218 dB for pistol shrimp
  // We compute this based on bubble collapse rate (if available from state)
  // For now, we'll compute it in the diagnostics module where we have access to bubble state
  
  return { 
    dPhaseDt, 
    Pacoustic, 
    gradient, 
    laplacian,
    jetPressure: params.jetDriven ? jetPressure : undefined,
    jetActive: params.jetDriven ? jetActive : undefined,
    shockwavePressure: params.jetDriven ? shockwavePressure : undefined,
    shockwaveIntensity: params.jetDriven ? shockwaveIntensity : undefined,
  };
}
