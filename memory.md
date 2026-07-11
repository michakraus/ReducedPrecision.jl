# ReducedPrecision — project summary

Analysis of geometric (symplectic) vs. non-geometric integrators run in varying floating-point
precision (Float16, Float32, Float64) on example problems from GeometricProblems.

## What was built

`ReducedPrecision` is a proper Julia package (`/Users/mkraus/Datashare/Julia/ReducedPrecision`):

- **`Project.toml`** — package name/uuid + `CairoMakie`, `GeometricSolutions`, `GeometricBase`
  added; all Geometric* deps dev-linked to sibling directories.
- **`src/ReducedPrecision.jl`** — reusable pipeline:
  - `PRECISIONS` and a method registry (`GEOMETRIC_METHODS` / `NONGEOMETRIC_METHODS`).
  - `run_study(make_problem)` — runs every method × precision, catching per-run failures.
  - `verify_precision` / `assert_precision` — the precision-purity gate.
  - `energy_error` (reuses `compute_invariant_error`), `solution_error`.
  - `plot_energy_error` / `plot_solution_error` (CairoMakie, one panel per precision).
- **`scripts/{harmonic_oscillator,pendulum,double_pendulum}.jl`** — short-step drivers
  (Δt = 0.1, t ≤ 100).
- **`scripts/{harmonic_oscillator,pendulum,double_pendulum}_longtime.jl`** — coarse-step /
  long-horizon drivers (Δt = 1, t ≤ 10 000) for a long-time stability study.
- **`plots/`** — generated figures. Each study writes, per method group, an energy-error, a
  solution-error, and a 2D solution (trajectory) figure.

Run any problem with: `julia --project=. scripts/<problem>.jl`.

## Plotting

Three plot types, all sharing the grid layout (one panel per precision):
- `plot_energy_error` — relative energy error vs. time (log y).
- `plot_solution_error` — state error vs. reference over time (log y).
- `plot_solution` — 2D trajectory of each method: phase space `(q, p)` for a one-DOF system,
  configuration space `(q₁, q₂)` otherwise (override via the `coords`/`xlabel`/`ylabel` kwargs).

Conventions:
- **Method groups.** Every plot function writes *two* files, appending `_euler` / `_other` to
  the given `path`: the Euler group (`EULER_METHODS` = explicit/implicit Euler + symplectic
  Euler A/B) and the rest (`OTHER_METHODS` = implicit/explicit midpoint + Crank-Nicolson + RK4).
  This grouping is orthogonal to the geometric/non-geometric flag, which still sets the line
  style (geometric = solid, non-geometric = dashed). Colors are globally consistent per method
  (Wong 8-colour palette `_PALETTE`, one entry per method in `ALL_METHODS`).
- **Legend below.** The legend is a horizontal row beneath the panels (not a side column);
  figure width is `430·np` and height `500`, keeping the per-panel size while making room. On
  the solution plots the legend also carries a `reference` entry.
- **Exact x-limits.** Time-series panels fit the x-axis exactly to the problem's `timespan`
  `[t₀, t₁]` (via `timespan(prob)`, not the accumulated grid endpoint, which rounds slightly
  short/long at low precision and would drop the final tick).
- **Overflow-safe y-limits, capped at 1e5.** Non-geometric energy errors can reach ~1e308
  (finite), which would both overflow Makie's log autolimit padding to `Inf` and dwarf the
  scale; `_plot_grid` sets explicit limits with the upper limit capped at **1e5**, so runaway
  lines are clipped at the top.
- **Reference trajectory.** `plot_solution` draws the bounded `reference` as a black backdrop in
  every panel and fits the axes to it, so each method's deviation is visible and runaway methods
  are clipped to the region of the true solution.
- `timevalues` uses the nominal grid `t₀ + n·Δt` (not the stored clock `sol.t`), because the
  low-precision clock saturates at long horizons (Float16: `t + Δt == t` past the representable
  integer range).

## Key design resolution

The literal method list can't run on one problem form (special `ExplicitEuler`/`ImplicitEuler` are
ODE-only; `SymplecticEulerA/B` are PODE/HODE-only). Resolution: use the **partitioned form** (PODE
for oscillator/pendulum, HODE for double pendulum) plus the numerically-identical RK Euler twins
`ExplicitEulerRK` / `ImplicitEulerRK`, which auto-promote to partitioned RK — so a single problem
form runs the full geometric-vs-non-geometric comparison.

Methods compared:
- Geometric: Symplectic Euler A, Symplectic Euler B, Implicit Midpoint.
- Non-geometric: Explicit Euler, Implicit Euler, Explicit Midpoint, Crank-Nicolson, RK4
  (explicit Runge–Kutta of order 4).

Solution-error reference: analytic `exact_solution` for the harmonic oscillator; Float64 `Gauss(8)`
(same grid) for the pendulum and double pendulum.

## Verification results

**Precision purity confirmed** — `verify_precision` passes for every successful run across all three
problems: `datatype`, `timetype`, and the stored `q`/`p` element types all equal the requested
precision. No library (GeometricIntegrators/Base, Solutions, Equations, Base, SimpleSolvers)
implicitly promotes to Float64, including for the double pendulum's hand-built Float16/Float32
construction.

**Plots show the expected physics:** symplectic methods keep bounded/oscillating energy error while
explicit Euler blows up and implicit Euler dissipates; the energy-conserving methods' error floor
drops cleanly per precision (e.g. oscillator: Float16 ~1e-2 → Float32 ~1e-6 → Float64 ~1e-15).

## Genuine reduced-precision findings (not bugs)

- **Float16 + long horizon:** with `t` up to 1000, Float16 cannot resolve successive time stamps
  (ulp(1000) ≈ 0.5 ≫ Δt = 0.1), breaking the implicit methods' Hermite initial guess
  (`t₀ == t₁`). Horizon capped at `t = 100` (nt = 1000, still many periods) so the full
  method × precision matrix is populated.
- **Float16 double pendulum:** the three implicit methods fail with "NaN in direction vector" — a
  real Float16 instability for this stiff, dimensional (g = 9.8), chaotic system. `run_study`
  catches these per-run so the sweep completes and reports them as skips.

All runs are wrapped so a single failure never aborts the study.

## Long-time variant findings (Δt = 1, t ≤ 10 000)

- **Harmonic oscillator & pendulum:** the contrast is dramatic at the coarse step — explicit
  Euler and explicit midpoint diverge exponentially (energy → ~1e300 in Float64, clipped at the
  top of the plot), while symplectic Euler A/B stay bounded over the full horizon and implicit
  midpoint / Crank-Nicolson stay near machine level. The energy floor still drops per precision.
  In Float16 some implicit methods fail (time-grid saturation → identical successive times).
- **Double pendulum:** Δt = 1 is roughly one natural period per step, far too coarse. All implicit
  methods fail at every precision ("NaN in direction vector") and the Gauss(8) reference also
  fails, so only the energy-error plot is produced (the solution-error plot is skipped by the
  guard). The explicit/symplectic methods blow up almost immediately. The plot documents this
  breakdown rather than a meaningful comparison — the short-step script is the informative one
  for this problem.
