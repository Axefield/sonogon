# Advanced ODE Integrators for Sonoluminescence

## Overview

The sonoluminescence model requires robust numerical integration due to extreme gradients near bubble collapse. This document describes the integrator implementations available in `core/integrator.ts`.

## The Challenge: Extreme Gradients

During bubble collapse, the system exhibits:

- **Extreme temporal gradients**: |∂R/∂t| can exceed sound speed, |∂T/∂t| ~ 10¹² K/s
- **Multiple timescales**: Acoustic period (~50 μs) vs. collapse (~100 ns) vs. light emission (~1 ns)
- **Stiff behavior**: Fast processes (plasma formation) coupled with slow processes (acoustic cycle)

A fixed-step integrator would either:
- Use small steps everywhere (inefficient)
- Use large steps (misses collapse dynamics, becomes unstable)

**Solution**: Adaptive step size control with error estimation.

## Available Integrators

### 1. Fixed-Step RK4 (`integrateRK4`)

**Use when**: Quick simulations, known smooth behavior, or when adaptive control is not needed.

```typescript
const result = integrateRK4(rhs, x0, {
  dt: 1e-9,    // Fixed time step
  tMax: 1e-5   // Maximum time
});
```

**Pros**:
- Simple and fast
- Predictable step count
- No overhead from error estimation

**Cons**:
- Cannot adapt to extreme gradients
- May require very small steps globally
- Can miss important dynamics

### 2. Dormand-Prince 5(4) (`integrateDormandPrince`)

**Use when**: Production simulations with extreme gradients, unknown behavior, or when accuracy is critical.

```typescript
const result = integrateDormandPrince(rhs, x0, {
  dt: 1e-9,           // Initial time step
  tMax: 1e-5,
  dtMin: 1e-15,       // Minimum step (prevents infinite loops)
  dtMax: 1e-6,        // Maximum step (prevents missing dynamics)
  tolerance: 1e-6,     // Error tolerance
  maxSteps: 1e6       // Safety limit
});
```

**How it works**:
1. Computes both 5th-order and 4th-order solutions
2. Estimates local truncation error from the difference
3. Adjusts step size to keep error below tolerance
4. Automatically uses smaller steps near collapse, larger steps during expansion

**Pros**:
- **5th order accuracy** (high precision)
- **Automatic step size control** (efficient)
- **Error estimation** (quality control)
- **Robust** (handles extreme gradients)

**Cons**:
- More expensive per step (7 function evaluations vs. 4)
- Step size varies (less predictable)
- Can reject steps (slight overhead)

### 3. Adaptive with Event Detection (`integrateAdaptive`)

**Use when**: You need to detect specific events (e.g., minimum radius, extreme gradients).

```typescript
// Event function: returns 0 when event occurs
const minRadiusEvent = (t: number, x: Float64Array) => {
  const state = mapper.fromVector(x, t);
  return state.hydro.R - 1e-6; // Event when R < 1 micron
};

const result = integrateAdaptive(rhs, x0, {
  dt: 1e-9,
  tMax: 1e-5,
  tolerance: 1e-6,
  adaptive: true  // Use Dormand-Prince
}, [minRadiusEvent]);

// Access detected events
if (result.events) {
  result.events.forEach(event => {
    console.log(`Event ${event.eventIndex} at t=${event.time}s`);
  });
}
```

**Features**:
- Detects when event functions cross zero
- Refines event location using bisection
- Records event times and states

## Recommended Usage

### For Sonoluminescence Simulations

```typescript
import { integrateAdaptive } from './core/integrator';
import { SonoluminescenceModel } from './model/sonoluminescenceModel';
import { DefaultStateVectorMapper } from './core/statevector';

const mapper = new DefaultStateVectorMapper();
const model = new SonoluminescenceModel(mapper, params);
const x0 = mapper.toVector(initialState);

// Use adaptive integrator with event detection
const result = integrateAdaptive(
  (t, x) => model.rhs(t, x),
  x0,
  {
    dt: 1e-9,        // Start with 1 ns steps
    tMax: 5e-5,      // 5 acoustic periods
    dtMin: 1e-15,    // Allow down to femtoseconds
    dtMax: 1e-7,     // But not larger than 100 ns
    tolerance: 1e-6, // 0.0001% relative error
    adaptive: true,
    maxSteps: 1e7    // Safety limit
  },
  [
    // Detect minimum radius
    (t, x) => {
      const state = mapper.fromVector(x, t);
      return state.hydro.R - 0.5e-6; // 0.5 micron
    },
    // Detect extreme temperature
    (t, x) => {
      const state = mapper.fromVector(x, t);
      return 20000 - state.gas.T; // 20,000 K threshold
    }
  ]
);

// Check statistics
console.log(`Steps: ${result.stats.steps}`);
console.log(`Rejected: ${result.stats.rejectedSteps}`);
console.log(`Min dt: ${result.stats.minDt}s`);
console.log(`Max dt: ${result.stats.maxDt}s`);
```

## Step Size Adaptation

The adaptive integrator adjusts step size based on local truncation error:

```
error_norm = ||x5 - x4|| / (tolerance * scale)
```

If `error_norm ≤ 1`: Accept step, increase dt
If `error_norm > 1`: Reject step, decrease dt

The step size adjustment:
```
dt_new = dt_old * safety * (1/error_norm)^(1/order)
```

Where:
- `safety = 0.9` (conservative factor)
- `order = 5` for Dormand-Prince

## Performance Characteristics

### Fixed-Step RK4
- **Function evaluations per step**: 4
- **Order**: 4
- **Step size**: Fixed
- **Best for**: Smooth regions, known behavior

### Dormand-Prince 5(4)
- **Function evaluations per step**: 7 (but steps are larger)
- **Order**: 5
- **Step size**: Adaptive (typically 10-1000x larger than fixed-step)
- **Best for**: Extreme gradients, unknown behavior, production runs

### Typical Step Sizes

During bubble expansion (smooth):
- Fixed-step: `dt = 1e-9` (1 ns)
- Adaptive: `dt ≈ 1e-7` to `1e-6` (100 ns to 1 μs)

During bubble collapse (extreme gradients):
- Fixed-step: `dt = 1e-9` (may be too large!)
- Adaptive: `dt ≈ 1e-12` to `1e-15` (1 ps to 1 fs)

**Result**: Adaptive integrator is often **10-100x faster** overall while maintaining accuracy.

## Error Control

The tolerance parameter controls accuracy:

- `tolerance = 1e-3`: Fast, less accurate (~0.1% error)
- `tolerance = 1e-6`: Balanced (default, ~0.0001% error)
- `tolerance = 1e-9`: High accuracy (~0.0000001% error, slower)

For sonoluminescence, `1e-6` is typically sufficient.

## Event Detection

Event functions should be smooth and cross zero:

```typescript
// Good: smooth crossing
const event = (t, x) => {
  const state = mapper.fromVector(x, t);
  return state.hydro.R - threshold;
};

// Bad: discontinuous
const badEvent = (t, x) => {
  const state = mapper.fromVector(x, t);
  return state.hydro.R < threshold ? -1 : 1; // Never crosses zero!
};
```

The integrator uses bisection to refine event locations to high precision.

## Limitations and Future Improvements

### Current Limitations

1. **No implicit methods**: For extremely stiff systems, implicit methods (e.g., BDF) might be needed
2. **Simple event detection**: Uses bisection; could use more sophisticated root finding
3. **No dense output**: Interpolation between steps is linear

### Potential Future Additions

1. **Implicit methods**: BDF or Rosenbrock methods for stiff systems
2. **Dense output**: Hermite interpolation for smooth output
3. **Parallelization**: GPU acceleration for large state vectors
4. **Symplectic methods**: For energy-conserving systems
5. **Multi-step methods**: Adams-Bashforth/Moulton for smooth regions

## Comparison with Other Integrators

| Method | Order | Adaptive | Stiff? | Best For |
|--------|-------|----------|--------|----------|
| RK4 (fixed) | 4 | No | No | Smooth, known behavior |
| Dormand-Prince | 5 | Yes | No | Extreme gradients, general use |
| BDF (future) | Variable | Yes | Yes | Very stiff systems |
| Symplectic (future) | 4 | No | No | Energy-conserving systems |

For sonoluminescence, **Dormand-Prince 5(4)** is the recommended choice.

## References

1. Dormand, J. R., & Prince, P. J. (1980). A family of embedded Runge-Kutta formulae. *Journal of Computational and Applied Mathematics*, 6(1), 19-26.

2. Hairer, E., Norsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems*. Springer.

3. Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing*. Cambridge University Press.

