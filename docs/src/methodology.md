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

### Partitioned Gauss(2) midpoint variants

A third comparison group holds four flavours of the 2-stage Gauss (partitioned midpoint) rule that
are algebraically the same method but differ in implementation detail:

| Method | Construction |
|:--|:--|
| `PRK Gauss(2)`  | `PartitionedTableau(Gauss(2))` — the Gauss tableau duplicated for `q` and `p` |
| `SPRK Gauss(2)` | `SymplecticPartitionedTableau(Gauss(2))` — the `p`-tableau is the symplectic conjugate, so symplecticity holds to floating-point accuracy by construction |
| `PRK Gauss(2), â=b̂=ĉ=0`  | as `PRK`, with the rounding-error compensation coefficients zeroed |
| `SPRK Gauss(2), â=b̂=ĉ=0` | as `SPRK`, with the rounding-error compensation coefficients zeroed |

These are all symplectic (drawn solid) and let the study isolate how the symplectic-vs-duplicated
tableau construction and the compensated-summation coefficients `â, b̂, ĉ` affect energy conservation
in reduced precision.

### Variational integrators for the degenerate-Lagrangian problems

The Lotka–Volterra problems are degenerate Lagrangian systems (IODE/LODE), on which the explicit,
symplectic-Euler and DIRK methods above are undefined. They are compared instead with several
flavours of the implicit midpoint rule that do apply: `Implicit Midpoint`, `VPRK(Gauss(1))`,
`PMVImidpoint` and (2D only) `CMDVI`. `VPRK(Gauss(1))` is rebuilt at the run precision so it stays
type-pure; `CMDVI` is omitted for the 4D system, where it fails to converge.

### Nonlinear solver

The implicit methods' stage equations are solved with the trust-region **`DogLeg`** solver of
`SimpleSolvers` (rather than a line-search Newton iteration), which is more robust in reduced
precision. Explicit methods carry no solver. The solver is a keyword of `run_study` /
`integrate_bounded` (default `DogLeg()`), so `Newton()` can be selected for comparison.

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
| Lotka–Volterra 2D | `lodeproblem` (`LotkaVolterra2dSingular`) | hand-built T-typed inputs | Float64 `Gauss(8)` |
| Lotka–Volterra 4D | `lodeproblem` (`LotkaVolterra4dLagrangian`, `A_quasicanonical_reduced`) | hand-built T-typed inputs | Float64 `Gauss(8)` |

The double pendulum, Toda lattice and both Lotka–Volterra problems are generated symbolically by
`EulerLagrange` and have no `::Type{T}` constructor, so `make_problem(T)` builds the initial
conditions, timespan, timestep and parameters at precision `T` explicitly. The Toda lattice
additionally carries a lattice size `N` (here `N = 16`, kept small so the sweep — in particular the
high-order reference and the implicit solves — stays tractable) and a Hamiltonian
`hamiltonian(t, q, p, params, N)` that takes `N`, so a closure is passed to the energy-error routine.

The Lotka–Volterra problems are **degenerate Lagrangian** systems posed as LODEs; the 4D case uses
the quasi-canonical reduced gauge matrix `A_quasicanonical_reduced` (with the exact one-form `B`) so
the discrete system is non-singular.

## Scenarios

Each problem is run in two scenarios, a fine short-horizon run and a coarser one:

* **Short / fine step** — harmonic oscillator, pendulum and Toda lattice use `Δt = 0.1`, `t ≤ 100`;
  the double pendulum uses `Δt = 0.01`, `t ≤ 10` (its natural timescale is much shorter); the
  Lotka–Volterra problems use `Δt = 0.01`, `t ≤ 10`.
* **Coarse step** — harmonic oscillator, pendulum and Toda lattice use `Δt = 1`, `t ≤ 10 000`; the
  double pendulum uses `Δt = 0.1`, `t ≤ 1000`; the Lotka–Volterra problems use `Δt = 0.1`,
  `t ≤ 100`.

The output figure filenames encode the timestep (e.g. `…_dt=0.1_…`), so the two scenarios of a
problem are distinguished by `Δt` rather than by a "longtime" label.

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

* Every figure has **one panel per precision** (Float16 / Float32 / Float64) and is produced once
  per method group. The four Hamiltonian problems use three groups — **Euler** (`_euler`), **other**
  (`_other`), and the **partitioned Gauss(2)** midpoint variants (`_midpoint`); the Lotka–Volterra
  problems use a single **variational** group (`_variational`). Scripts pass their group set to the
  plotting routines, which colour methods consistently within it.
* Error plots use a **shared logarithmic y-axis** across all three panels, with the upper limit
  **capped at `1e5`** so runaway (non-geometric) errors are clipped rather than dominating the
  scale; the x-axis is fitted **exactly** to the integration interval.
* Solution plots show the **2D trajectory** of each method (phase space `(q, p)`, or configuration
  space for multi-degree-of-freedom systems) with the **reference** drawn as a black backdrop.
* Colours are consistent per method; the legend is a horizontal row beneath the panels.
