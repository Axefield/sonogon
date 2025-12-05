"use strict";
// analysis/diagnostics.ts
// Diagnostic tools for monitoring bubble dynamics and detecting extreme conditions
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeEnergyBudget = computeEnergyBudget;
exports.detectExtremeGradients = detectExtremeGradients;
exports.estimateRefractiveIndex = estimateRefractiveIndex;
exports.computePlasmaDiagnostics = computePlasmaDiagnostics;
exports.computeStabilityMetrics = computeStabilityMetrics;
exports.computeAtomicSubatomicDisturbances = computeAtomicSubatomicDisturbances;
exports.computeModeSqueezingMetrics = computeModeSqueezingMetrics;
const observables_1 = require("./observables");
const units_1 = require("../core/units");
/**
 * Reaction enthalpies [J/mol] at standard conditions
 * Positive = endothermic (requires energy)
 * Negative = exothermic (releases energy)
 */
const REACTION_ENTHALPIES = {
    // Reaction 0: H2O ↔ H + OH
    H2O_dissociation: 500e3 * units_1.Constants.R_gas, // ~500 kJ/mol
    // Reaction 1: O2 ↔ 2O
    O2_dissociation: 498e3 * units_1.Constants.R_gas, // ~498 kJ/mol
    // Reaction 2: H + O ↔ OH (recombination, exothermic)
    H_O_recombination: -428e3 * units_1.Constants.R_gas, // ~-428 kJ/mol (releases energy)
};
/**
 * Compute chemical energy change from reaction progress
 *
 * Chemical energy change = Σ(ΔH_i · dξ_i/dt · dt)
 * where ΔH_i is the reaction enthalpy and dξ_i/dt is the reaction rate
 */
function computeChemicalEnergy(state, statePrev, dt) {
    const { reactions } = state;
    const { reactions: reactionsPrev } = statePrev;
    const { T } = state.gas;
    // Compute reaction progress changes
    const dxi0 = reactions.xi[0] - reactionsPrev.xi[0];
    const dxi1 = reactions.xi.length > 1 ? reactions.xi[1] - reactionsPrev.xi[1] : 0;
    const dxi2 = reactions.xi.length > 2 ? reactions.xi[2] - reactionsPrev.xi[2] : 0;
    // Convert reaction progress to number of moles (simplified: use volume)
    const { R } = state.hydro;
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const volumePrev = (4.0 / 3.0) * Math.PI * statePrev.hydro.R * statePrev.hydro.R * statePrev.hydro.R;
    const avgVolume = (volume + volumePrev) / 2.0;
    // Approximate number of moles (using ideal gas law)
    const n_moles = (state.gas.Pg * avgVolume) / (units_1.Constants.R_gas * T);
    // Energy change from each reaction
    // Reaction 0: H2O ↔ H + OH (endothermic)
    const E0 = REACTION_ENTHALPIES.H2O_dissociation * dxi0 * n_moles;
    // Reaction 1: O2 ↔ 2O (endothermic)
    const E1 = REACTION_ENTHALPIES.O2_dissociation * dxi1 * n_moles;
    // Reaction 2: H + O ↔ OH (exothermic, releases energy)
    const E2 = REACTION_ENTHALPIES.H_O_recombination * dxi2 * n_moles;
    // Total chemical energy change (positive = energy absorbed, negative = energy released)
    return E0 + E1 + E2;
}
/**
 * Compute energy budget for the bubble system
 *
 * Tracks energy flows to detect conservation and identify
 * dominant energy pathways.
 */
function computeEnergyBudget(state, statePrev, dt, acousticPower) {
    const { R, Rdot } = state.hydro;
    const { Pg, T } = state.gas;
    const { ne, Te } = state.plasma;
    const { em } = state;
    // Acoustic work input (simplified: acoustic power * time)
    const acousticWorkIn = acousticPower * dt;
    // Compression work: P * dV
    // dV = 4*π*R²*dR = 4*π*R²*Rdot*dt
    const dV = 4.0 * Math.PI * R * R * Rdot * dt;
    const compressionWork = Pg * dV;
    // Heat loss (simplified: proportional to surface area and temperature difference)
    const T_ambient = 293.15; // 20°C
    const heatLossCoeff = 100; // W/(m²·K) - approximate
    const surfaceArea = 4.0 * Math.PI * R * R;
    const heatLoss = heatLossCoeff * (T - T_ambient) * surfaceArea * dt;
    // Light emission
    const emission = (0, observables_1.estimateEmission)(state);
    const lightEmission = emission.totalPower * dt;
    // Viscous dissipation (simplified: ~mu * (Rdot/R)² * volume)
    const mu = 1e-3; // Water viscosity [Pa·s]
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const strainRate = Math.abs(Rdot / Math.max(R, 1e-10));
    const viscousDissipation = mu * strainRate * strainRate * volume * dt;
    // Chemical energy: compute from reaction progress and enthalpies
    const chemicalEnergy = computeChemicalEnergy(state, statePrev, dt);
    // EM stored energy
    const emStoredEnergy = em.storedEnergy;
    // Total energy (should be approximately conserved over cycles)
    const totalEnergy = acousticWorkIn -
        compressionWork -
        heatLoss -
        lightEmission -
        viscousDissipation -
        chemicalEnergy -
        emStoredEnergy;
    return {
        acousticWorkIn,
        compressionWork,
        heatLoss,
        lightEmission,
        viscousDissipation,
        chemicalEnergy,
        emStoredEnergy,
        totalEnergy,
    };
}
/**
 * Detect extreme gradient events
 *
 * Identifies when the system enters non-classical regimes
 * with extreme spatial or temporal gradients.
 */
function detectExtremeGradients(statePrev, stateNext, dt) {
    const metrics = (0, observables_1.computeGradients)(statePrev, stateNext, dt);
    const spatialGradients = (0, observables_1.computeSpatialGradients)(stateNext);
    // Additional spatial gradient thresholds
    const extremePressureGradient = spatialGradients.pressureGradient > 1e18; // Pa/m
    const extremeDensityGradient = spatialGradients.densityGradient > 1e40; // 1/m⁴
    const detected = metrics.extremeGradient || extremePressureGradient || extremeDensityGradient;
    // Compute atomic/subatomic level disturbances when extreme gradients are detected
    // This connects extreme gradients to atomic-level disturbances that probe
    // deeper into subatomic structure, similar to how extreme conditions at CERN
    // reveal exotic hadron states
    const atomicDisturbances = detected
        ? computeAtomicSubatomicDisturbances(stateNext)
        : undefined;
    return {
        detected,
        metrics,
        spatialGradients,
        atomicDisturbances,
    };
}
/**
 * Estimate refractive index from bubble state
 *
 * Combines contributions from:
 * - Neutral gas (density-dependent)
 * - Plasma (frequency-dependent, can be negative)
 * - Ionization effects
 */
function estimateRefractiveIndex(state, opticalFrequency) {
    const { R } = state.hydro;
    const { T } = state.gas;
    const { ne } = state.plasma;
    const { numberDensity } = state.species;
    // Neutral gas refractive index (simplified: n ≈ 1 + α * ρ)
    // α is polarizability, ρ is density
    const totalNeutralDensity = Object.values(numberDensity).reduce((sum, n) => sum + n, 0);
    const alpha = 1e-29; // Approximate polarizability [m³]
    const n_neutral = 1.0 + alpha * totalNeutralDensity;
    // Plasma refractive index: n² = 1 - (ωp² / ω²)
    // For ωp > ω, n becomes imaginary (cutoff)
    const omega_p = units_1.Calculations.plasmaFrequency(ne);
    const n_plasma_squared = 1.0 - (omega_p * omega_p) / (opticalFrequency * opticalFrequency);
    const n_plasma = n_plasma_squared > 0 ? Math.sqrt(n_plasma_squared) : 0;
    // Effective refractive index (simplified combination)
    const n_effective = n_neutral + (n_plasma - 1.0) * 0.5; // Weighted average
    // Radial gradient (simplified: ~n / R)
    const R_safe = Math.max(R, 1e-10);
    const n_gradient = Math.abs(n_effective) / R_safe;
    return {
        n_effective,
        n_plasma,
        n_neutral,
        n_gradient,
    };
}
/**
 * Compute plasma diagnostics
 *
 * Analyzes plasma properties and their relationship to
 * optical frequencies and bubble size.
 */
function computePlasmaDiagnostics(state, opticalFrequency = 2 * Math.PI * 3e14 // ~500 nm
) {
    const { R } = state.hydro;
    const { T } = state.gas;
    const { ne, Te } = state.plasma;
    const plasmaFrequency = units_1.Calculations.plasmaFrequency(ne);
    const debyeLength = units_1.Calculations.debyeLength(ne, T);
    // Electron mean free path (simplified)
    // λ_mfp ≈ v_th / ν_coll, where v_th = sqrt(k_B*Te/m_e)
    const v_th = Math.sqrt((units_1.Constants.k_B * Te) / units_1.Constants.m_e);
    const nu_coll = 1e13; // Approximate collision frequency [1/s]
    const meanFreePath = v_th / nu_coll;
    const plasmaCutoff = plasmaFrequency >= opticalFrequency;
    const collisional = meanFreePath < R * 0.1; // Mean free path << bubble size
    return {
        plasmaFrequency,
        debyeLength,
        opticalFrequency,
        plasmaCutoff,
        meanFreePath,
        collisional,
    };
}
/**
 * Compute stability metrics
 *
 * Detects limit cycles, runaway heating, and energy conservation
 * to assess system stability.
 */
function computeStabilityMetrics(states, times) {
    if (states.length < 10) {
        // Need enough data points
        return {
            isStable: false,
            limitCycleDetected: false,
            runawayHeating: false,
            energyConserved: false,
        };
    }
    // Check for runaway heating
    const temperatures = states.map((s) => s.gas.T);
    const maxTemp = Math.max(...temperatures);
    const minTemp = Math.min(...temperatures);
    const tempTrend = temperatures[temperatures.length - 1] - temperatures[0];
    const runawayHeating = tempTrend > 1000 && maxTemp > 1e5; // Rapid temperature increase
    // Check for limit cycle (periodic behavior)
    // Simplified: look for repeating patterns in radius
    const radii = states.map((s) => s.hydro.R);
    const maxR = Math.max(...radii);
    const minR = Math.min(...radii);
    const amplitude = maxR - minR;
    // Count zero crossings of (R - R_mean) to estimate period
    const R_mean = radii.reduce((sum, r) => sum + r, 0) / radii.length;
    let zeroCrossings = 0;
    for (let i = 1; i < radii.length; i++) {
        if ((radii[i - 1] - R_mean) * (radii[i] - R_mean) < 0) {
            zeroCrossings++;
        }
    }
    const periodEstimate = zeroCrossings > 0
        ? (times[times.length - 1] - times[0]) / (zeroCrossings / 2)
        : undefined;
    const limitCycleDetected = zeroCrossings >= 4 && amplitude > 0.1 * R_mean; // At least 2 full cycles
    // Energy conservation check (simplified)
    // Energy should be approximately constant over a cycle
    const energies = states.map((s) => {
        const volume = (4.0 / 3.0) * Math.PI * s.hydro.R * s.hydro.R * s.hydro.R;
        return s.gas.Pg * volume + s.em.storedEnergy; // Simplified energy
    });
    const energyVariance = energies.reduce((sum, e) => {
        const mean = energies.reduce((s, v) => s + v, 0) / energies.length;
        return sum + (e - mean) * (e - mean);
    }, 0) / energies.length;
    const energyMean = energies.reduce((sum, e) => sum + e, 0) / energies.length;
    const energyConserved = energyVariance / (energyMean * energyMean) < 0.1; // <10% variance
    const isStable = !runawayHeating && limitCycleDetected && energyConserved;
    return {
        isStable,
        limitCycleDetected,
        runawayHeating,
        energyConserved,
        periodEstimate,
    };
}
/**
 * Nuclear and particle physics constants for subatomic disturbance calculations
 */
const NUCLEAR_PHYSICS = {
    // Nuclear matter density (saturation density)
    NUCLEAR_MATTER_DENSITY: 2.3e17, // kg/m³
    // Temperature conversions
    K_TO_MEV: 1.0 / (1.16e10), // 1 MeV = 1.16e10 K
    MEV_TO_K: 1.16e10, // 1 MeV = 1.16e10 K
    // Strong force range
    STRONG_FORCE_RANGE: 1e-15, // 1 fm (typical strong force range)
    // Quark deconfinement temperature
    DECONFINEMENT_TEMPERATURE_MEV: 150, // MeV (~150 MeV for QCD phase transition)
    DECONFINEMENT_TEMPERATURE_K: 1.74e12, // K
    // Nuclear binding energies (typical values)
    BINDING_ENERGY_PER_NUCLEON_MEV: 8.0, // MeV (typical for medium-mass nuclei)
    IONIZATION_POTENTIAL_Ar_MEV: 15.76e-6, // MeV (Argon ionization, converted)
    // Strong coupling constant (varies with energy scale)
    // At low energies (~1 GeV): α_s ≈ 0.3-1.0
    // At high energies (~100 GeV): α_s ≈ 0.1
    ALPHA_STRONG_LOW: 0.5, // At ~1 GeV scale
    ALPHA_STRONG_HIGH: 0.1, // At ~100 GeV scale
    // LHC/CERN comparison values
    LHC_CENTER_OF_MASS_ENERGY_MEV: 1.3e7, // 13 TeV = 1.3e7 MeV
    LHC_COLLISION_TEMPERATURE_MEV: 300, // Estimated from heavy-ion collisions
    // Exotic hadron masses (reference)
    TETRAQUARK_MASS_RANGE_MEV: [3800, 7000], // MeV (typical all-charm tetraquarks)
    PENTAQUARK_MASS_RANGE_MEV: [4200, 4500], // MeV (typical pentaquarks)
    HEXAQUARK_MASS_RANGE_MEV: [1800, 2400], // MeV (dibaryons, e.g., d*)
    HYBRID_MESON_MASS_RANGE_MEV: [1400, 2000], // MeV (quark + gluon hybrid)
    GLUEBALL_MASS_RANGE_MEV: [1500, 3000], // MeV (pure gluon states)
    // Nuclear radius constant
    NUCLEAR_RADIUS_CONSTANT: 1.2e-15, // m (r0 ≈ 1.2 fm)
    // Fermi momentum (nuclear matter)
    FERMI_MOMENTUM_MEV_C: 250, // MeV/c (typical for nuclear matter)
    // Quark masses (MeV/c²)
    QUARK_MASS_UP_MEV: 2.3, // MeV
    QUARK_MASS_DOWN_MEV: 4.8, // MeV
    QUARK_MASS_STRANGE_MEV: 95, // MeV
    QUARK_MASS_CHARM_MEV: 1270, // MeV
    QUARK_MASS_BOTTOM_MEV: 4180, // MeV
    QUARK_MASS_TOP_MEV: 173000, // MeV
    // QCD scale (confinement scale)
    LAMBDA_QCD_MEV: 200, // MeV (typical QCD scale)
    // Nuclear fusion
    GAMOW_PEAK_ENERGY_MEV: 0.1, // MeV (typical for light nuclei)
    FUSION_CROSS_SECTION_REF: 1e-47, // m² (reference at Gamow peak)
    // Planck scale (quantum gravity)
    PLANCK_DENSITY: 5.16e96, // kg/m³
    PLANCK_LENGTH: 1.616e-35, // m
    PLANCK_TIME: 5.391e-44, // s
    PLANCK_ENERGY_MEV: 1.22e19, // MeV (Planck energy)
    // QCD vacuum energy density (bag constant)
    BAG_CONSTANT_MEV4: 150, // MeV⁴ (typical value)
    BAG_CONSTANT_J_M3: 150 * Math.pow(1.602e-13, 4) / Math.pow(1.97e-16, 3), // J/m³ (converted)
    // Gluon field strength (typical)
    GLUON_FIELD_STRENGTH_REF: 1e15, // T (tesla equivalent, very rough)
    // Color charge (QCD)
    COLOR_CHARGE_STRENGTH: 1.0, // Normalized (strong coupling handles scaling)
};
/**
 * Compute atomic/subatomic level disturbance metrics
 *
 * At extreme gradients during collapse, conditions approach regimes where:
 * - Nuclear densities are reached (~10¹⁸ kg/m³)
 * - Temperatures reach MeV scales (nuclear binding energies)
 * - Strong force effects may become relevant
 * - Atomic structure is probed at fundamental levels
 *
 * This function tracks conditions where exotic hadron-like states might form
 * or where quark-level physics might manifest, similar to the tetraquark
 * discoveries at CERN where extreme conditions reveal multi-quark bound states.
 *
 * The connection: extreme gradients in sonoluminescence create transient
 * conditions where atomic-level disturbances probe deeper into subatomic structure.
 *
 * EXPANDED VERSION: Now includes detailed nuclear physics, quark deconfinement
 * conditions, exotic hadron formation probabilities, comparison to accelerator conditions,
 * quark flavor tracking, quantum field effects, nuclear fusion conditions, momentum space
 * distributions, and Planck-scale physics.
 */
function computeAtomicSubatomicDisturbances(state) {
    const { R, Rdot } = state.hydro;
    const { Pg, T } = state.gas;
    const { ne } = state.plasma;
    // ========== MASS DENSITY CALCULATIONS ==========
    const volume = (4.0 / 3.0) * Math.PI * R * R * R;
    const n_total = Object.values(state.species.numberDensity).reduce((sum, n) => sum + n, 0);
    // Average molecular weight (weighted by species density)
    // For Ar-dominated: ~40 amu; for mixed: compute weighted average
    const m_avg = 40 * 1.66e-27; // kg (argon atom mass, simplified)
    const massDensity = n_total * m_avg; // kg/m³
    const nuclearDensityRatio = massDensity / NUCLEAR_PHYSICS.NUCLEAR_MATTER_DENSITY;
    const nuclearDensityReached = nuclearDensityRatio > 0.01; // 1% of nuclear density
    // ========== TEMPERATURE METRICS ==========
    const temperatureMeV = T * NUCLEAR_PHYSICS.K_TO_MEV;
    const mevTemperatureReached = temperatureMeV > 0.1; // 100 keV threshold
    const temperatureRatioToDeconfinement = temperatureMeV / NUCLEAR_PHYSICS.DECONFINEMENT_TEMPERATURE_MEV;
    const approachingQuarkDeconfinement = temperatureRatioToDeconfinement > 0.1; // 10% of deconfinement temp
    const deconfinementProximity = Math.min(temperatureRatioToDeconfinement, 1.0);
    // ========== STRONG FORCE CALCULATIONS ==========
    // Strong coupling constant depends on energy scale (temperature)
    // Use running coupling: α_s decreases with energy
    const energyScaleMeV = Math.max(temperatureMeV, 1.0); // Minimum 1 MeV
    // Simplified running: α_s(Q) ≈ α_s(Q0) * log(Q0/Λ) / log(Q/Λ)
    // For rough estimate, interpolate between low and high energy values
    const logEnergy = Math.log(energyScaleMeV);
    const logLow = Math.log(1000); // 1 GeV reference
    const logHigh = Math.log(100000); // 100 GeV reference
    const alpha_s = logEnergy < logLow
        ? NUCLEAR_PHYSICS.ALPHA_STRONG_LOW
        : logEnergy > logHigh
            ? NUCLEAR_PHYSICS.ALPHA_STRONG_HIGH
            : NUCLEAR_PHYSICS.ALPHA_STRONG_LOW +
                (NUCLEAR_PHYSICS.ALPHA_STRONG_HIGH - NUCLEAR_PHYSICS.ALPHA_STRONG_LOW) *
                    (logEnergy - logLow) / (logHigh - logLow);
    const strongForceRelevant = nuclearDensityRatio > 0.1 &&
        temperatureMeV > 0.5 &&
        R < 10 * NUCLEAR_PHYSICS.STRONG_FORCE_RANGE;
    // Strong force field estimate: F ~ α_s * (energy density) * (range)²
    const energyDensity = (3.0 / 2.0) * n_total * units_1.Constants.k_B * T; // J/m³
    const strongForceFieldEstimate = alpha_s * energyDensity *
        NUCLEAR_PHYSICS.STRONG_FORCE_RANGE * NUCLEAR_PHYSICS.STRONG_FORCE_RANGE;
    // ========== QUARK-LEVEL CONDITIONS ==========
    const quarkLevelConditions = nuclearDensityRatio > 0.1 &&
        temperatureMeV > 1.0 &&
        Pg > 1e15;
    // ========== EXOTIC HADRON FORMATION (EXPANDED) ==========
    // Tetraquark formation probability (heuristic)
    // Based on: density, temperature, and compression rate
    // Similar to CMS conditions: high density, MeV temperatures, transient compression
    const densityFactor = Math.min(nuclearDensityRatio, 1.0);
    const temperatureFactor = Math.min(temperatureMeV / 10.0, 1.0); // Normalize to ~10 MeV
    const compressionFactor = Math.min(Math.abs(Rdot) / 10000, 1.0); // Normalize to ~10 km/s
    const tetraquarkFormationProbability = densityFactor * temperatureFactor * compressionFactor * 0.1; // Max 10% (very rough)
    // Pentaquark formation (typically rarer)
    const pentaquarkFormationProbability = tetraquarkFormationProbability * 0.1; // ~10x rarer
    // Hexaquark formation (dibaryons, e.g., d*)
    // Typically require even higher densities and specific conditions
    const hexaquarkFormationProbability = tetraquarkFormationProbability * 0.01; // ~100x rarer than tetraquarks
    // Hybrid meson formation (quark + gluon hybrid states)
    // Requires strong gluon field coupling
    const gluonCouplingFactor = Math.min(alpha_s / 0.5, 1.0); // Normalize to typical α_s
    const hybridMesonFormationProbability = densityFactor * temperatureFactor * gluonCouplingFactor * 0.05;
    // Glueball formation (pure gluon states)
    // Very rare, requires extreme gluon field conditions
    const glueballFormationProbability = gluonCouplingFactor * temperatureFactor * 0.001; // Very rare
    // Multi-quark state formation rate (estimated)
    // Rate ~ (density)² * (temperature) * (collision cross section)
    const collisionCrossSection = Math.PI * NUCLEAR_PHYSICS.STRONG_FORCE_RANGE * NUCLEAR_PHYSICS.STRONG_FORCE_RANGE;
    const thermalVelocity = Math.sqrt(3.0 * units_1.Constants.k_B * T / m_avg); // m/s
    const multiQuarkStateFormationRate = n_total * n_total * thermalVelocity * collisionCrossSection *
        Math.exp(-NUCLEAR_PHYSICS.BINDING_ENERGY_PER_NUCLEON_MEV * NUCLEAR_PHYSICS.MEV_TO_K * units_1.Constants.k_B / (units_1.Constants.k_B * T)); // Boltzmann factor
    const exoticStateFormationPossible = nuclearDensityRatio > 0.05 &&
        temperatureMeV > 0.5 &&
        Math.abs(Rdot) > 1000;
    // ========== QUARK FLAVOR TRACKING ==========
    // Determine which quark flavors can be produced at this energy scale
    const lightQuarkEnergy = temperatureMeV; // Energy available for quark production
    const strangeQuarkThreshold = temperatureMeV > NUCLEAR_PHYSICS.QUARK_MASS_STRANGE_MEV * 2; // Need ~2x mass for pair production
    const charmQuarkThreshold = temperatureMeV > NUCLEAR_PHYSICS.QUARK_MASS_CHARM_MEV * 2;
    const heavyQuarkProductionPossible = charmQuarkThreshold || temperatureMeV > NUCLEAR_PHYSICS.QUARK_MASS_BOTTOM_MEV * 2;
    const quarkFlavorRelevant = temperatureMeV > NUCLEAR_PHYSICS.QUARK_MASS_UP_MEV;
    // ========== QUANTUM FIELD EFFECTS ==========
    // QCD vacuum energy (bag constant)
    const qcdVacuumEnergy = NUCLEAR_PHYSICS.BAG_CONSTANT_J_M3; // J/m³
    // Vacuum fluctuations become relevant when energy density approaches QCD scale
    // (energyDensity already computed above for strong force calculations)
    const vacuumFluctuationsRelevant = energyDensity > 0.1 * qcdVacuumEnergy;
    // Casimir effect estimate (simplified: E ~ hbar * c / R⁴)
    const casimirEffectEstimate = units_1.Constants.hbar * units_1.Constants.c / (R * R * R * R); // J (very rough)
    // Gluon field strength estimate
    // F ~ α_s * (energy density) / (strong force range)
    const gluonFieldStrength = alpha_s * energyDensity / NUCLEAR_PHYSICS.STRONG_FORCE_RANGE; // N/m² ≈ T (rough)
    // Color charge density (normalized)
    // Estimate based on quark density and strong coupling
    const colorChargeDensity = n_total * alpha_s; // Normalized (very rough)
    // ========== NUCLEAR FUSION CONDITIONS ==========
    // Fusion requires: high density, high temperature, sufficient energy to overcome Coulomb barrier
    // Gamow factor: exp(-2π * Z1*Z2*e² / (hbar * v))
    const Z1 = 18; // Argon atomic number (simplified)
    const Z2 = 18;
    const coulombBarrier = Z1 * Z2 * units_1.Constants.e * units_1.Constants.e / (4.0 * Math.PI * units_1.Constants.epsilon_0 * NUCLEAR_PHYSICS.STRONG_FORCE_RANGE); // J
    const gamowFactor = Math.exp(-2.0 * Math.PI * Z1 * Z2 * units_1.Constants.e * units_1.Constants.e /
        (units_1.Constants.hbar * thermalVelocity * 4.0 * Math.PI * units_1.Constants.epsilon_0)); // Penetration probability
    // Fusion cross-section (simplified: σ ~ σ_ref * Gamow * exp(-E_Gamow / kT))
    const E_Gamow = NUCLEAR_PHYSICS.GAMOW_PEAK_ENERGY_MEV * 1.602e-13; // J
    const fusionCrossSection = NUCLEAR_PHYSICS.FUSION_CROSS_SECTION_REF * gamowFactor *
        Math.exp(-E_Gamow / (units_1.Constants.k_B * T));
    // Fusion rate: rate ~ n² * v * σ
    const fusionRate = n_total * n_total * thermalVelocity * fusionCrossSection;
    const fusionConditionsApproached = nuclearDensityRatio > 0.1 &&
        temperatureMeV > 0.1 &&
        gamowFactor > 1e-10;
    // ========== MOMENTUM SPACE DISTRIBUTIONS ==========
    // Fermi energy: E_F = (hbar² / (2m)) * (3π² * n)^(2/3)
    const fermiMomentum = units_1.Constants.hbar * Math.pow(3.0 * Math.PI * Math.PI * n_total, 1.0 / 3.0); // kg·m/s
    const fermiEnergy = (fermiMomentum * fermiMomentum) / (2.0 * m_avg); // J
    // Degeneracy pressure: P_deg = (3π²)^(2/3) * (hbar² / (5m)) * n^(5/3)
    const degeneracyPressure = Math.pow(3.0 * Math.PI * Math.PI, 2.0 / 3.0) *
        (units_1.Constants.hbar * units_1.Constants.hbar / (5.0 * m_avg)) *
        Math.pow(n_total, 5.0 / 3.0); // Pa
    // Quantum degeneracy parameter: n * λ³, where λ = h / (2π * m * v_th) is thermal de Broglie wavelength
    const thermalDeBroglieWavelength = units_1.Constants.h / (2.0 * Math.PI * m_avg * thermalVelocity); // m
    const degeneracyParameter = n_total * thermalDeBroglieWavelength * thermalDeBroglieWavelength * thermalDeBroglieWavelength;
    // Quantum statistics relevant when degeneracy parameter > 1
    const quantumStatisticsRelevant = degeneracyParameter > 1.0;
    // ========== PLANCK-SCALE PHYSICS ==========
    const planckDensityRatio = massDensity / NUCLEAR_PHYSICS.PLANCK_DENSITY;
    const planckLengthRatio = R / NUCLEAR_PHYSICS.PLANCK_LENGTH;
    // Characteristic time scale: τ ~ R / v_sound or R / v_thermal
    const v_sound = Math.sqrt(units_1.Constants.k_B * T / m_avg); // Sound speed estimate
    const characteristicTime = R / Math.max(v_sound, Math.abs(Rdot));
    const planckTimeRatio = characteristicTime / NUCLEAR_PHYSICS.PLANCK_TIME;
    const approachingPlanckScale = planckDensityRatio > 1e-30 || // Even tiny fraction is interesting
        planckLengthRatio < 1e10; // Within 10 orders of magnitude
    // Quantum gravity effects possible when approaching Planck scale
    const quantumGravityEffectsPossible = planckDensityRatio > 1e-20 ||
        planckLengthRatio < 1e5;
    // ========== NUCLEAR STRUCTURE EFFECTS ==========
    const nuclearBindingEnergy = NUCLEAR_PHYSICS.BINDING_ENERGY_PER_NUCLEON_MEV * NUCLEAR_PHYSICS.MEV_TO_K * units_1.Constants.k_B; // J
    const thermalEnergyPerParticle = (3.0 / 2.0) * units_1.Constants.k_B * T; // J
    const nuclearBindingEnergyRatio = thermalEnergyPerParticle / nuclearBindingEnergy;
    const nuclearStructureAffected = nuclearBindingEnergyRatio > 0.1; // Thermal energy > 10% of binding energy
    // Estimated nuclear radius if matter were compressed to nuclear density
    // Using liquid drop model: r = r0 * A^(1/3), where A is mass number
    // For our case, estimate A from density
    const A_estimate = Math.pow(massDensity / (NUCLEAR_PHYSICS.NUCLEAR_MATTER_DENSITY * 0.6), 3); // Rough estimate
    const nuclearRadius = NUCLEAR_PHYSICS.NUCLEAR_RADIUS_CONSTANT * Math.pow(A_estimate, 1.0 / 3.0);
    // ========== ACCELERATOR COMPARISONS ==========
    // LHC: 13 TeV center-of-mass energy per collision
    // Per particle: ~6.5 TeV = 6.5e6 MeV
    const lhcEnergyPerParticleMeV = NUCLEAR_PHYSICS.LHC_CENTER_OF_MASS_ENERGY_MEV / 2.0;
    const lhcEnergyRatio = temperatureMeV / lhcEnergyPerParticleMeV; // Very small, but tracks approach
    // RHIC: ~200 GeV per nucleon = 2e5 MeV
    const rhicEnergyPerParticleMeV = 2e5;
    const rhicEnergyRatio = temperatureMeV / rhicEnergyPerParticleMeV;
    // SPS: ~17 GeV per nucleon = 1.7e4 MeV
    const spsEnergyPerParticleMeV = 1.7e4;
    const spsEnergyRatio = temperatureMeV / spsEnergyPerParticleMeV;
    const comparableToLHC = lhcEnergyRatio > 1e-6; // Even 1 ppm is interesting
    const acceleratorConditions = {
        comparableToLHC: lhcEnergyRatio > 1e-6,
        comparableToRHIC: rhicEnergyRatio > 1e-4,
        comparableToSPS: spsEnergyRatio > 1e-3,
    };
    // ========== QCD PHASE DIAGRAM ==========
    // Determine QCD phase based on temperature and density
    let qcdPhase = 'unknown';
    if (temperatureRatioToDeconfinement > 1.0 && nuclearDensityRatio > 0.1) {
        qcdPhase = 'quark-gluon-plasma';
    }
    else if (temperatureRatioToDeconfinement > 0.5 && nuclearDensityRatio > 0.05) {
        qcdPhase = 'crossover'; // In crossover region
    }
    else if (nuclearDensityReached || mevTemperatureReached) {
        qcdPhase = 'hadronic'; // Still in hadronic phase but extreme
    }
    const qcdPhaseDiagramPosition = {
        temperatureRatio: temperatureRatioToDeconfinement,
        densityRatio: nuclearDensityRatio,
    };
    // ========== PROBING DEPTH (EXPANDED) ==========
    let probingDepth = 'none';
    if (approachingPlanckScale || quantumGravityEffectsPossible) {
        probingDepth = 'planck';
    }
    else if (quarkLevelConditions || approachingQuarkDeconfinement || qcdPhase === 'quark-gluon-plasma') {
        probingDepth = 'quark';
    }
    else if (strongForceRelevant || nuclearStructureAffected) {
        probingDepth = 'nuclear';
    }
    else if (nuclearDensityReached || mevTemperatureReached) {
        probingDepth = 'atomic';
    }
    const atomicStructureProbing = probingDepth !== 'none';
    // ========== TIME SCALES (EXPANDED) ==========
    // Nuclear time scale: ~10⁻²² s (typical nuclear process)
    const nuclearTimeScale = 1e-22; // s
    // Strong interaction time scale: ~10⁻²³ s (strong force interaction)
    const strongInteractionTimeScale = 1e-23; // s
    // Fusion time scale: ~1 / fusion_rate (if fusion is occurring)
    const fusionTimeScale = fusionRate > 0 ? 1.0 / fusionRate : Infinity; // s
    // Planck time scale
    const planckTimeScale = NUCLEAR_PHYSICS.PLANCK_TIME; // s
    return {
        // Basic nuclear conditions
        nuclearDensityReached,
        nuclearDensityRatio,
        massDensity,
        // Temperature metrics
        mevTemperatureReached,
        temperatureMeV,
        temperatureRatioToDeconfinement,
        // Strong force metrics
        strongForceRelevant,
        strongForceFieldEstimate,
        strongCouplingConstant: alpha_s,
        strongForceRange: NUCLEAR_PHYSICS.STRONG_FORCE_RANGE,
        // Quark-level conditions
        quarkLevelConditions,
        approachingQuarkDeconfinement,
        deconfinementProximity,
        // Exotic state formation (expanded)
        exoticStateFormationPossible,
        tetraquarkFormationProbability,
        pentaquarkFormationProbability,
        hexaquarkFormationProbability,
        hybridMesonFormationProbability,
        glueballFormationProbability,
        multiQuarkStateFormationRate,
        // Quark flavor tracking
        quarkFlavorRelevant,
        lightQuarkEnergy,
        charmQuarkThreshold,
        strangeQuarkThreshold,
        heavyQuarkProductionPossible,
        // Quantum field effects
        qcdVacuumEnergy,
        vacuumFluctuationsRelevant,
        casimirEffectEstimate,
        gluonFieldStrength,
        colorChargeDensity,
        // Nuclear fusion conditions
        fusionConditionsApproached,
        fusionCrossSection,
        gamowFactor,
        fusionRate,
        // Momentum space distributions
        fermiEnergy,
        degeneracyPressure,
        quantumStatisticsRelevant,
        fermiMomentum,
        degeneracyParameter,
        // Planck-scale physics
        approachingPlanckScale,
        planckDensityRatio,
        planckLengthRatio,
        planckTimeRatio,
        quantumGravityEffectsPossible,
        // Nuclear structure effects
        nuclearStructureAffected,
        nuclearBindingEnergyRatio,
        nuclearRadius,
        // Comparison to accelerators
        comparableToLHC,
        lhcEnergyRatio,
        acceleratorConditions,
        // QCD phase diagram
        qcdPhase,
        qcdPhaseDiagramPosition,
        // Atomic structure probing
        atomicStructureProbing,
        probingDepth,
        // Time scales (expanded)
        nuclearTimeScale,
        strongInteractionTimeScale,
        fusionTimeScale,
        planckTimeScale,
    };
}
/**
 * Compute mode squeezing metrics
 *
 * Analyzes EM mode amplitudes to detect parametric amplification
 * and mode squeezing effects. Reports pumped energy into modes vs released energy.
 *
 * This quantifies the negative-space behavior: energy pumped into squeezed states
 * during compression vs energy released as photons during decay.
 */
function computeModeSqueezingMetrics(state, statePrev, dt) {
    const { modes, storedEnergy } = state.em;
    const { R, Rdot } = state.hydro;
    const modeAmplitudes = modes.map((m) => Math.sqrt(m.re * m.re + m.im * m.im));
    const totalModeEnergy = modes.reduce((sum, m) => sum + (m.re * m.re + m.im * m.im), 0);
    // Squeezing detected if mode amplitudes are large relative to stored energy
    // (indicating parametric amplification)
    const squeezingDetected = totalModeEnergy > storedEnergy * 0.1 && storedEnergy > 1e-20;
    // Pump efficiency: ratio of stored energy to input (simplified)
    const pumpEfficiency = storedEnergy > 0 ? totalModeEnergy / storedEnergy : 0;
    // PUMPED ENERGY: Energy pumped into modes via parametric amplification
    // During compression (Rdot < 0), parametric Hamiltonian pumps energy into modes
    // Estimate: pumped energy ≈ coupling_strength * |Rdot| * |gradient| * dt
    let pumpedEnergy = 0;
    if (statePrev && dt !== undefined) {
        const dE_em = storedEnergy - statePrev.em.storedEnergy;
        // If E_em increased, that's pumped energy (positive contribution)
        // If E_em decreased, that's released energy (negative contribution)
        if (dE_em > 0) {
            pumpedEnergy = dE_em;
        }
    }
    else {
        // Estimate from current state: pumped energy ~ coupling * |Rdot| * storedEnergy
        // This is a rough estimate when previous state not available
        const couplingEstimate = Math.abs(Rdot) / Math.max(R, 1e-10);
        pumpedEnergy = couplingEstimate * storedEnergy * 1e-9; // Rough estimate
    }
    // RELEASED ENERGY: Energy released from modes as photons
    // This is the decay term: released = E_em / τ * dt
    const decayTime = 1e-9; // 1 ns
    const releasedEnergy = dt !== undefined
        ? (storedEnergy / decayTime) * dt
        : storedEnergy / decayTime * 1e-9; // Default dt estimate
    // NET ENERGY CHANGE: pumped - released
    const netEnergyChange = pumpedEnergy - releasedEnergy;
    return {
        modeAmplitudes,
        totalModeEnergy,
        squeezingDetected,
        pumpEfficiency,
        pumpedEnergy,
        releasedEnergy,
        netEnergyChange,
    };
}
//# sourceMappingURL=diagnostics.js.map