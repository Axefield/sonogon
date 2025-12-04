# Kinematic Sciences Enhancement Plan

## Overview
Focus on enhancing the kinematic and dynamic aspects of bubble motion, including shape oscillations, translation, and acoustic field coupling.

## Current State
- Only spherical bubble (R, Rdot)
- No shape deformations
- No bubble translation
- Uniform acoustic field (no gradients)

## Planned Enhancements

### 1. Shape Oscillations
**Goal**: Model non-spherical bubble shapes using spherical harmonics

**Implementation**:
- Add P2 (quadrupole) and P4 (hexadecapole) modes
- Shape: R(θ,φ) = R₀ + a₂*P₂(cos θ) + a₄*P₄(cos θ)
- Mode amplitudes: a₂, a₄ (and their velocities)
- Restoring forces from surface tension
- Damping from viscosity

**State Variables**:
- a₂, a₂_dot (quadrupole mode)
- a₄, a₄_dot (hexadecapole mode)

### 2. Bubble Translation
**Goal**: Model bubble position and motion in acoustic field

**Implementation**:
- Position: x, y, z (or r, θ, φ in spherical)
- Velocity: vx, vy, vz
- Acceleration from Bjerknes forces
- Drag forces

**State Variables**:
- x, y, z (position)
- vx, vy, vz (velocity)

### 3. Bjerknes Forces
**Goal**: Forces on bubble from acoustic field gradients

**Primary Bjerknes Force**:
- F = -V * ∇P_acoustic
- Drives bubble to pressure nodes/antinodes

**Secondary Bjerknes Force**:
- Interaction between bubbles
- F = -ρ * V₁ * V₂ * ∇²P / (4π * r²)

### 4. Acoustic Field Gradients
**Goal**: Spatial variation of acoustic pressure

**Implementation**:
- Standing wave patterns
- Gradient computation: ∇P_acoustic
- Laplacian: ∇²P_acoustic

### 5. Shape-Radial Coupling
**Goal**: Couple shape oscillations to radial dynamics

**Implementation**:
- Shape affects effective radius
- Radial motion affects shape stability
- Energy exchange between modes

