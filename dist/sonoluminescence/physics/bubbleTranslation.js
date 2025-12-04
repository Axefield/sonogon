"use strict";
// physics/bubbleTranslation.ts
// Bubble translation dynamics and Bjerknes forces
Object.defineProperty(exports, "__esModule", { value: true });
exports.computeBubbleTranslationDerivatives = computeBubbleTranslationDerivatives;
exports.computeSecondaryBjerknesForce = computeSecondaryBjerknesForce;
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
function computePrimaryBjerknesForce(V, // Bubble volume [m³]
gradP // Pressure gradient [Pa/m]
) {
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
function computeDragForce(R, // Bubble radius [m]
v, // Velocity [m/s]
mu, // Viscosity [Pa·s]
rho, // Density [kg/m³]
dragCoeff) {
    const v_mag = Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    if (dragCoeff !== undefined) {
        // Custom drag coefficient
        const A = Math.PI * R * R; // Cross-sectional area
        const F_mag = 0.5 * dragCoeff * rho * A * v_mag * v_mag;
        if (v_mag > 0) {
            return {
                fx: -F_mag * (v.x / v_mag),
                fy: -F_mag * (v.y / v_mag),
                fz: -F_mag * (v.z / v_mag),
            };
        }
    }
    else {
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
function computeBubbleTranslationDerivatives(state, params) {
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
    const m = rho * V; // Effective mass (includes added mass)
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
    // 3. Buoyancy (optional, if gravity included)
    // F_buoyancy = (ρ_liquid - ρ_gas) * V * g
    // For now, assume neutral buoyancy (gas density << liquid)
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
function computeSecondaryBjerknesForce(V, // Bubble volume [m³]
V_other, // Other bubble volume [m³]
r, // Distance between bubbles [m]
laplacianP, // ∇²P_acoustic [Pa/m²]
rho // Liquid density [kg/m³]
) {
    const r_safe = Math.max(r, 1e-10);
    return -(rho * V * V_other * laplacianP) / (4.0 * Math.PI * r_safe * r_safe);
}
//# sourceMappingURL=bubbleTranslation.js.map