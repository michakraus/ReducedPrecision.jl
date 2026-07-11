# Methodology

## Precisions

Every experiment is run at three precisions:

```julia
const PRECISIONS = (Float16, Float32, Float64)
```

## Methods

The methods are split into a *geometric* (symplectic) group and a *non-geometric* group. The line
style in every plot follows this classification — **solid** for geometric, **dashed** for
non-geometric.

| Method | Group | Order | Notes |
|:--|:--|:--:|:--|
| Symplectic Euler A | geometric | 1 | symplectic partitioned Euler |
| Symplectic Euler B | geometric | 1 | symplectic partitioned Euler |
| Implicit Midpoint  | geometric | 2 | symplectic, symmetric |
| Explicit Euler     | non-geometric | 1 | energy-increasing |
| Implicit Euler     | non-geometric | 1 | energy-dissipating |
| Explicit Midpoint  | non-geometric | 2 | |
| Crank–Nicolson     | non-geometric | 2 | trapezoidal rule |
| RK4                | non-geometric | 4 | explicit Runge–Kutta |

### One problem form for all methods

The literal method list cannot be run on a single problem form: the special `ExplicitEuler` /
`ImplicitEuler` of `GeometricIntegratorsBase` are ODE-only, while `SymplecticEulerA` / `B` require a
partitioned (PODE/HODE) problem. The resolution is to use the **partitioned form** for every
problem and to represent explicit/implicit Euler by their numerically identical Runge–Kutta tableau
twins `ExplicitEulerRK` / `ImplicitEulerRK`, which auto-promote to partitioned Runge–Kutta on a
PODE/HODE. A single partitioned problem then runs the entire method set, giving a fair
geometric-vs-non-geometric comparison.

## Problems

| Problem | Form | Precision handling | Reference solution |
|:--|:--|:--|:--|
| Harmonic oscillator | `podeproblem` | `::Type{T}` constructor | analytic `exact_solution` |
| Pendulum | `podeproblem` | `::Type{T}` constructor | Float64 `Gauss(8)` |
| Double pendulum | `hodeproblem` | hand-built T-typed inputs | Float64 `Gauss(8)` |
| Toda lattice (N = 16) | `hodeproblem` | hand-built T-typed inputs | Float64 `Gauss(8)` |

The double pendulum and Toda lattice are generated symbolically by `EulerLagrange` and have no
`::Type{T}` constructor, so `make_problem(T)` builds the initial conditions, timespan, timestep and
parameters at precision `T` explicitly. The Toda lattice additionally carries a lattice size `N`
(here `N = 16`, kept small so the sweep — in particular the high-order reference and the implicit
solves — stays tractable) and a Hamiltonian `hamiltonian(t, q, p, params, N)` that takes `N`, so a
closure is passed to the energy-error routine.

## Scenarios

Each problem is run in two scenarios:

* **Short / fine step** — harmonic oscillator, pendulum and Toda lattice use `Δt = 0.1`, `t ≤ 100`;
  the double pendulum uses `Δt = 0.01`, `t ≤ 10` (its natural timescale is much shorter).
* **Long / coarse step** — `Δt = 1`, `t ≤ 10 000`, to stress long-time stability.

## Type-purity verification

Precision is set by the problem in two independent type parameters: the **data type**
(`datatype`, from `eltype(q₀)`) drives the state and cache allocations, and the **time type**
(`timetype`, from `promote(t₀, t₁, Δt)`) drives the tableau coefficients. To stay pure in `T`,
*both* the initial conditions and the timespan/timestep must be `T`.

For every successful run, `verify_precision` asserts that

* `datatype(prob) === timetype(prob) === T`,
* `datatype(sol) === timetype(sol) === T`, and
* the element type of the stored `q` and `p` arrays is `T`,

which proves that no library in the stack silently promotes to `Float64`. Integration failures
(e.g. half-precision blow-ups) are caught per run so a single failure never aborts the sweep.

## Error metrics

* **Relative energy error** `|(H(qₙ, pₙ) − H₀) / H₀|`, computed in the run's own precision.
* **Solution error** — the Euclidean norm of the state difference against the reference solution on
  the same time grid.

## Plotting conventions

* Every figure has **one panel per precision** (Float16 / Float32 / Float64) and is produced twice,
  once for the **Euler** method group (`_euler`) and once for the **other** methods (`_other`), so
  each figure stays readable.
* Error plots use a **shared logarithmic y-axis** across all three panels, with the upper limit
  **capped at `1e5`** so runaway (non-geometric) errors are clipped rather than dominating the
  scale; the x-axis is fitted **exactly** to the integration interval.
* Solution plots show the **2D trajectory** of each method (phase space `(q, p)`, or configuration
  space for multi-degree-of-freedom systems) with the **reference** drawn as a black backdrop.
* Colours are consistent per method; the legend is a horizontal row beneath the panels.
