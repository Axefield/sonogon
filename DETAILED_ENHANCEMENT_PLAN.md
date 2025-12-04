# Detailed Physics Enhancement Plan

## Overview
This plan focuses on adding detailed physics to replace simplified approximations throughout the model.

## Areas Requiring Detailed Implementation

### 1. Hydrodynamics Details
**Current Issues**:
- Keller-Miksis uses hardcoded gamma=1.4 instead of from params
- No shape oscillations (only spherical assumption)
- Simplified dPg/dt in compressibility term

**Enhancements Needed**:
- Use actual gamma from thermo params
- Add shape oscillation modes (P2, P4 spherical harmonics)
- Proper time derivative calculation in Keller-Miksis
- Add Bjerknes force for bubble translation

### 2. Thermodynamics Details
**Current Issues**:
- Heat loss calculation is very simplified (hardcoded coefficient)
- Van der Waals implementation is basic (no iterative solution)
- No detailed heat capacity (Cp, Cv) calculations
- No thermal conductivity effects

**Enhancements Needed**:
- Detailed heat transfer model (conduction, convection, radiation)
- Proper van der Waals iterative solution
- Temperature-dependent heat capacity
- Thermal boundary layer effects

### 3. Plasma Physics Details
**Current Issues**:
- Collision cross sections are hardcoded approximations
- Recombination rate is simplified (just n_e²)
- Ionization rate uses simple relaxation (not detailed rate equations)
- No multi-step ionization
- No detailed collision integrals

**Enhancements Needed**:
- Temperature-dependent collision cross sections
- Detailed recombination (radiative, three-body, dielectronic)
- Proper ionization rate equations (not just relaxation)
- Collision integrals for energy exchange
- Multi-step ionization processes

### 4. EM Cavity Details
**Current Issues**:
- Mode coupling is simplified (no cross-mode interactions)
- Pumping term is simplified (just Rdot * gradient)
- Quality factor Q is hardcoded
- No detailed cavity geometry effects

**Enhancements Needed**:
- Cross-mode coupling (mode mixing)
- Detailed parametric pumping (proper Hamiltonian)
- Frequency-dependent Q factor
- Cavity geometry effects (spherical harmonics)
- Quantum corrections (squeezing, entanglement)

### 5. Chemical Reactions Details
**Current Issues**:
- Rate constants are simplified (Arrhenius only)
- No pressure-dependent rates
- Three-body reactions are simplified
- No detailed reaction mechanisms

**Enhancements Needed**:
- Pressure-dependent rate constants (Lindemann mechanism)
- Detailed three-body collision rates
- Reaction mechanism networks
- Temperature-dependent pre-exponential factors
- Quantum tunneling corrections

### 6. Species Transport Details
**Current Issues**:
- No diffusion effects
- No mass transport
- Species just scale with volume

**Enhancements Needed**:
- Diffusion through bubble wall
- Mass transport equations
- Species-dependent diffusion coefficients
- Concentration gradients

### 7. Internal Energy Details
**Current Issues**:
- Relaxation times are simplified
- Energy exchange is linear relaxation only
- No detailed mode-specific physics

**Enhancements Needed**:
- Detailed Landau-Teller relaxation
- Mode-specific relaxation times
- Quantum effects for vibrational modes
- Detailed electronic state populations

### 8. Acoustic Details
**Current Issues**:
- No acoustic field coupling
- No standing wave effects
- No bubble translation

**Enhancements Needed**:
- Acoustic field gradients
- Standing wave patterns
- Bubble translation dynamics
- Secondary Bjerknes forces

## Implementation Plan

### Phase 1: Critical Physics Details
1. Fix Keller-Miksis to use actual gamma
2. Detailed heat transfer model
3. Proper ionization rate equations
4. Detailed recombination rates

### Phase 2: Enhanced Models
5. Shape oscillations
6. Cross-mode EM coupling
7. Pressure-dependent reaction rates
8. Detailed collision integrals

### Phase 3: Advanced Effects
9. Quantum corrections
10. Diffusion effects
11. Acoustic field coupling
12. Multi-step ionization

