/**
 * Fundamental physical constants (SI units)
 */
export declare const Constants: {
    readonly R_gas: 8.314462618;
    readonly k_B: 1.380649e-23;
    readonly N_A: 6.02214076e+23;
    readonly e: 1.602176634e-19;
    readonly m_e: 9.1093837015e-31;
    readonly m_p: 1.67262192369e-27;
    readonly epsilon_0: 8.8541878128e-12;
    readonly mu_0: 0.00000125663706212;
    readonly c: 299792458;
    readonly h: 6.62607015e-34;
    readonly hbar: 1.054571817e-34;
    readonly sigma_SB: 5.670374419e-8;
    readonly alpha_s_QCD: 0.1181;
    readonly nuclearBindingEnergyPerNucleon: number;
    readonly deconfinementTemperature: 150000000;
    readonly nuclearMatterDensity: 230000000000000000;
    readonly nuclearRadiusConstant: 1.2e-15;
    readonly fermiMomentumNuclear: 1.36e-19;
    readonly verdet_TGG_633nm: 134;
    readonly verdet_TGG_532nm: 50;
    readonly verdet_TGG_1064nm: 40;
    readonly verdet_TGG_1550nm: 8;
};
/**
 * Common physical properties (SI units)
 */
export declare const Properties: {
    readonly water: {
        readonly density: 998.2;
        readonly viscosity: 0.001002;
        readonly surfaceTension: 0.0728;
        readonly vaporPressure: 2339;
        readonly speedOfSound: 1482;
    };
    readonly air: {
        readonly density: 1.204;
        readonly viscosity: 0.00001825;
    };
    readonly argon: {
        readonly molecularWeight: 0.039948;
        readonly gamma: 1.67;
    };
    readonly xenon: {
        readonly molecularWeight: 0.131293;
        readonly gamma: 1.67;
    };
};
/**
 * Reference values for normalization
 */
export declare const ReferenceValues: {
    readonly R0_ref: 0.000001;
    readonly t0_ref: 0.000001;
    readonly P0_ref: 101325;
    readonly T0_ref: 293.15;
    readonly rho0_ref: 998.2;
    readonly E0_ref: 1e-18;
};
/**
 * Unit conversion utilities
 */
export declare const Units: {
    /**
     * Convert pressure from Pa to normalized units
     */
    pressureToNormalized(P: number): number;
    /**
     * Convert pressure from normalized units to Pa
     */
    pressureFromNormalized(P_norm: number): number;
    /**
     * Convert temperature from K to normalized units
     */
    temperatureToNormalized(T: number): number;
    /**
     * Convert temperature from normalized units to K
     */
    temperatureFromNormalized(T_norm: number): number;
    /**
     * Convert length from m to normalized units
     */
    lengthToNormalized(L: number): number;
    /**
     * Convert length from normalized units to m
     */
    lengthFromNormalized(L_norm: number): number;
    /**
     * Convert time from s to normalized units
     */
    timeToNormalized(t: number): number;
    /**
     * Convert time from normalized units to s
     */
    timeFromNormalized(t_norm: number): number;
    /**
     * Convert number density from 1/m³ to normalized units
     */
    numberDensityToNormalized(n: number): number;
    /**
     * Convert number density from normalized units to 1/m³
     */
    numberDensityFromNormalized(n_norm: number): number;
};
/**
 * Helper functions for common calculations
 */
export declare const Calculations: {
    /**
     * Ideal gas law: P = n * k_B * T
     */
    idealGasPressure(n: number, T: number): number;
    /**
     * Number density from pressure and temperature: n = P / (k_B * T)
     */
    numberDensityFromPressure(P: number, T: number): number;
    /**
     * Plasma frequency: ωp = sqrt(ne * e² / (ε0 * me))
     */
    plasmaFrequency(ne: number): number;
    /**
     * Debye length: λD = sqrt(ε0 * k_B * T / (ne * e²))
     */
    debyeLength(ne: number, T: number): number;
    /**
     * Blackbody power density: P = σ * T⁴
     */
    blackbodyPowerDensity(T: number): number;
    /**
     * Bremsstrahlung power density (simplified): P ~ ne² * sqrt(Te)
     */
    bremsstrahlungPowerDensity(ne: number, Te: number): number;
    /**
     * Verdet constant for Terbium-Gallium-Garnet (TGG) as function of wavelength
     *
     * TGG was used in the Hebrew University research on magnetic field effects in light.
     * Verdet constant determines Faraday rotation: θ = V(λ) * B₀ * L
     *
     * @param wavelength Wavelength in meters [m]
     * @returns Verdet constant in rad/(T·m)
     *
     * Reference values (from Northrop Grumman and literature):
     * - 532 nm (green): ~50 rad/(T·m)
     * - 633 nm (red): ~134 rad/(T·m) (negative sign indicates direction)
     * - 1064 nm (near-IR): ~40 rad/(T·m)
     * - 1550 nm (IR): ~8 rad/(T·m)
     *
     * For other wavelengths, uses approximate scaling: V(λ) ~ 1/λ² (dispersive)
     */
    verdetConstantTGG(wavelength: number): number;
    /**
     * Verdet constant for gas/plasma (much smaller than TGG)
     *
     * For gases, Verdet constant is typically ~0.01-0.1 rad/(T·m)
     * For plasma, depends on electron density and temperature
     *
     * @param wavelength Wavelength in meters [m]
     * @param electronDensity Electron density [1/m³] (for plasma)
     * @returns Verdet constant in rad/(T·m)
     */
    verdetConstantGas(wavelength: number, electronDensity?: number): number;
};
//# sourceMappingURL=units.d.ts.map