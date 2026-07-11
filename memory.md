# ReducedPrecision — project summary

Analysis of geometric (symplectic) vs. non-geometric integrators run in varying floating-point
precision (Float16, Float32, Float64) on example problems from GeometricProblems.

## What was built

`ReducedPrecision` is a proper Julia package (`/Users/mkraus/Datashare/Julia/ReducedPrecision`):

- **`Project.toml`** — package name/uuid + deps (`CairoMakie` + the Geometric* ecosystem),
  resolved from the **registry** (no `[sources]`; `[compat]` pins the working versions). For local
  work against the sibling checkouts, dev-link them into the git-ignored Manifest with
  `pkg> dev ../GeometricBase ../GeometricEquations …`. Test-only `Test` via `[extras]`/`[targets]`.
- **`src/`** — reusable pipeline, split into logical units and stitched together by
  `ReducedPrecision.jl` (usings/exports/`include`s only):
  - `methods.jl` — `PRECISIONS`, `MethodSpec`, the method registries, and the plotting groups
    (`EULER_METHODS` / `OTHER_METHODS` / `METHOD_GROUPS`).
  - `study.jl` — the `Run` type and `run_study(make_problem)` (runs every method × precision,
    catching per-run failures; the problem is built once per precision and reused across methods —
    important for the EulerLagrange-generated double pendulum and Toda lattice). Integration goes
    through `integrate_bounded` (step-by-step, replicating GeometricIntegrators' own loop, so
    results match `integrate` for well-behaved runs) with a **divergence guard**: if the state
    goes non-finite or exceeds `bound` (default `1e3`) it stops early and NaN-fills the tail, so
    runaway (typically non-convergent implicit / coarse-step explicit) runs don't waste the rest of
    the sweep or pollute the plots. `run.diverged` records the stop step; `nothing` if bounded.
  - `diagnostics.jl` — `verify_precision` / `assert_precision` (purity gate), `timevalues`,
    `energy_error` (reuses `compute_invariant_error`), `solution_error`.
  - `plotting.jl` — `plot_energy_error` / `plot_solution_error` / `plot_solution` (CairoMakie).
- **`test/runtests.jl`** — unit tests (registry, run_study + purity across all precisions,
  diagnostics, per-run failure capture, and that the plotting routines write both group figures).
  Run with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- **`scripts/{harmonic_oscillator,pendulum,double_pendulum,toda_lattice}.jl`** — short-step
  drivers (HO & pendulum: Δt = 0.1, t ≤ 100; double pendulum: Δt = 0.01, t ≤ 10; Toda lattice:
  Δt = 0.1, t ≤ 100).
- **`scripts/{harmonic_oscillator,pendulum,double_pendulum,toda_lattice}_longtime.jl`** —
  coarse-step / long-horizon drivers (Δt = 1, t ≤ 10 000) for a long-time stability study.
- **`plots/`** — generated figures. Each study writes, per method group, an energy-error, a
  solution-error, and a 2D solution (trajectory) figure.
- Every plot title carries the run parameters, e.g. `… (Δt = 0.1, t ≤ 100)`, so short- and
  long-time figures are distinguishable at a glance.

Run any problem with: `julia --project=. scripts/<problem>.jl`.

- **`docs/`** — Documenter.jl site summarising all experiments and findings (Home, Methodology, a
  page per problem, Findings). Build with `julia --project=docs docs/make.jl` *after* running the
  scripts — `make.jl` copies `plots/*.png` into `docs/src/figures/` (git-ignored) and embeds them.
  Docs depend only on `Documenter` (figures are pre-generated, not built via `@example`).
- **CI (`.github/workflows/`)** — split by concern:
  - `CI.yml` runs **only the tests** (matrix: Julia LTS `1.10` + latest stable `1` × ubuntu/macOS/
    windows; no `arch` pin since macOS runners are arm64). Deps resolve from the registry.
  - `Documenter.yml` **builds and deploys the docs**: `julia-buildpkg` instantiates the main
    project, a step runs all eight experiment scripts to (re)generate `plots/`, then
    `docs/make.jl` embeds them and deploys. Figures are never committed — regenerated each build.

## Plotting

Three plot types, all sharing the grid layout (one panel per precision):
- `plot_energy_error` — relative energy error vs. time (log y).
- `plot_solution_error` — state error vs. reference over time (log y).
- `plot_solution` — 2D trajectory of each method: phase space `(q, p)` for a one-DOF system,
  configuration space `(q₁, q₂)` otherwise (override via the `coords`/`xlabel`/`ylabel` kwargs).

Conventions:
- **Method groups.** Every plot function writes *two* files, appending `_euler` / `_other` to
  the given `path`: the Euler group (`EULER_METHODS` = symplectic Euler A/B + explicit/implicit
  Euler) and the rest (`OTHER_METHODS`, in draw/legend order: RK4, Explicit Midpoint, Implicit
  Midpoint, Crank-Nicolson). The order each group is *listed* is the order methods are drawn and
  appear in the legend; colours stay tied to the method name, so reordering does not recolour.
  This grouping is orthogonal to the geometric/non-geometric flag, which still sets the line
  style (geometric = solid, non-geometric = dashed). Colours are globally consistent per method
  (Wong 8-colour palette `_PALETTE`, one entry per method in `ALL_METHODS`).
- **Legend below.** The legend is a horizontal row beneath the panels (not a side column);
  figure width is `430·np` and height `500`, keeping the per-panel size while making room. On
  the solution plots the legend also carries a `reference` entry.
- **Exact x-limits.** Time-series panels fit the x-axis exactly to the problem's `timespan`
  `[t₀, t₁]` (via `timespan(prob)`, not the accumulated grid endpoint, which rounds slightly
  short/long at low precision and would drop the final tick).
- **Shared, overflow-safe y-limits, capped at 1e5.** All panels of an error figure share one
  y-range (`_shared_ylims`, computed from every panel's finite data) so precisions are directly
  comparable. Non-geometric energy errors can reach ~1e308 (finite), which would both overflow
  Makie's log autolimit padding to `Inf` and dwarf the scale, so the upper limit is capped at
  **1e5** and runaway lines are clipped at the top.
- **Reference trajectory.** `plot_solution` draws the bounded `reference` *first*, as a black
  backdrop (legend entry "Reference"), so every method lies on top of it; the axes are fitted to
  the reference so runaway methods are clipped to the region of the true solution.
- `timevalues` uses the nominal grid `t₀ + n·Δt` (not the stored clock `sol.t`), because the
  low-precision clock saturates at long horizons (Float16: `t + Δt == t` past the representable
  integer range).

## Key design resolution

The literal method list can't run on one problem form (special `ExplicitEuler`/`ImplicitEuler` are
ODE-only; `SymplecticEulerA/B` are PODE/HODE-only). Resolution: use the **partitioned form** (PODE
for oscillator/pendulum, HODE for double pendulum and Toda lattice) plus the numerically-identical
RK Euler twins `ExplicitEulerRK` / `ImplicitEulerRK`, which auto-promote to partitioned RK — so a
single problem form runs the full geometric-vs-non-geometric comparison.

Methods compared:
- Geometric: Symplectic Euler A, Symplectic Euler B, Implicit Midpoint.
- Non-geometric: Explicit Euler, Implicit Euler, Explicit Midpoint, Crank-Nicolson, RK4
  (explicit Runge–Kutta of order 4).

Problems (4): harmonic oscillator, pendulum (both have a `::Type{T}` `podeproblem`), double
pendulum and Toda lattice (both EulerLagrange-generated, no `::Type{T}` — `make_problem(T)`
hand-builds T-typed `q₀`, `p₀`, timespan, timestep and parameters). Two wrinkles specific to the
Toda lattice: a lattice size `N` (using `N = 16`, not the default 200, to keep the sweep
tractable) and a `hamiltonian(t, q, p, params, N)` that takes `N`, so the script passes a
`(t,q,p,params) -> hamiltonian(…, N)` closure to `plot_energy_error`; its trajectory plot uses the
phase space of the first lattice site.

Solution-error reference: analytic `exact_solution` for the harmonic oscillator; Float64 `Gauss(8)`
(same grid) for the pendulum, double pendulum and Toda lattice.

**Registry-vs-local caveat (resolved).** The registered `GeometricProblems` can lag the local
checkouts, causing signature drift. This bit the pendulum: registered ≤ 0.6.24 had only a
parameter-free `hamiltonian(t, q, p)`, while `energy_error` (via `compute_invariant_error`) passes
`params`, so it errored on CI (registry) but not locally (dev-linked). Resolved by releasing
`GeometricProblems` **0.6.25**, which restores the `(t, q, p, params)` method: `[compat]` now
requires `0.6.25` and the pendulum scripts pass `hamiltonian` directly again (like the others).
Heuristic: if a script errors only on CI but not locally, suspect this kind of signature drift first.

## Verification results

**Precision purity confirmed** — `verify_precision` passes for every successful run across all four
problems: `datatype`, `timetype`, and the stored `q`/`p` element types all equal the requested
precision. No library (GeometricIntegrators/Base, Solutions, Equations, Base, SimpleSolvers)
implicitly promotes to Float64, including for the double pendulum's and Toda lattice's hand-built
Float16/Float32 construction.

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
- **Toda lattice:** behaves like the oscillator/pendulum — implicit midpoint / Crank-Nicolson keep
  energy bounded (~1e-5) while explicit midpoint and RK4 drift; the Gauss(8) reference converges
  even at Δt = 1, so the full plot set is produced. In Float16 the two implicit methods fail on the
  long-horizon time-grid saturation; the short scenario runs every method at every precision (the
  bounded bump initial data keeps the exponentials well-behaved).
