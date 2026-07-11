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
  - `methods.jl` — `PRECISIONS`, `MethodSpec`, the method registries, and the plotting groups.
    `ALL_METHODS` (12) = `GEOMETRIC_METHODS` (3) + `NONGEOMETRIC_METHODS` (5) + `MIDPOINT_METHODS`
    (4 partitioned Gauss(2) variants), grouped for plotting by `METHOD_GROUPS`
    (`euler` / `other` / `midpoint`). The degenerate-Lagrangian Lotka–Volterra comparison uses
    separate sets `LV2D_METHODS` (4, incl. `CMDVI`) / `LV4D_METHODS` (3, no `CMDVI`) and
    `LV2D_GROUPS` / `LV4D_GROUPS` (a single `variational` group each). The `GaussVPRK` wrapper
    rebuilds `VPRK(Gauss(1))` at the run precision (its own `initmethod` otherwise bakes in Float64).
  - `study.jl` — the `Run` type and `run_study(make_problem; …, solver = DogLeg())` (runs every
    method × precision, catching per-run failures; the problem is built once per precision and
    reused across methods — important for the EulerLagrange-generated problems). Implicit methods
    use the trust-region **`DogLeg`** solver by default (more robust in reduced precision than
    line-search Newton; pass `solver = Newton()` to compare); explicit methods carry no solver (gated
    on `isimplicit`). Integration goes through `integrate_bounded` (step-by-step, replicating
    GeometricIntegrators' own loop, so results match `integrate` for well-behaved runs) with a
    **divergence guard**: if the state goes non-finite or exceeds `bound` (default `1e3`) it stops
    early and NaN-fills the tail. `run.diverged` records the stop step; `nothing` if bounded.
  - `diagnostics.jl` — `verify_precision` / `assert_precision` (purity gate), `timevalues`,
    `energy_error` (reuses `compute_invariant_error`), `solution_error`.
  - `plotting.jl` — `plot_energy_error` / `plot_solution_error` / `plot_solution` (CairoMakie).
- **`test/runtests.jl`** — unit tests (registry incl. the three plotting groups, run_study + purity
  across all precisions, diagnostics, per-run failure capture, and that the plotting routines write
  per-group figures incl. a custom group set). Run with `julia --project=. -e 'using Pkg; Pkg.test()'`.
- **`scripts/{harmonic_oscillator,pendulum,double_pendulum,toda_lattice,lotka_volterra_2d,
  lotka_volterra_4d}.jl`** — short-step drivers (HO & pendulum & Toda: Δt = 0.1, t ≤ 100; double
  pendulum & both Lotka–Volterra: Δt = 0.01, t ≤ 10).
- **`scripts/{…}_longtime.jl`** — coarse-step drivers: HO/pendulum/Toda Δt = 1, t ≤ 10 000; double
  pendulum Δt = 0.1, t ≤ 1000; Lotka–Volterra Δt = 0.1, t ≤ 100.
- **`plots/`** — generated figures. Each study writes, per method group, an energy-error, a
  solution-error, and a 2D solution (trajectory) figure. **Filenames encode the timestep**
  (`…_dt=<Δt>_<group>.png`), so a problem's two scenarios are distinguished by `Δt` (not a
  "longtime" label). Every plot title also carries the run parameters, e.g. `… (Δt = 0.1, t ≤ 100)`.

Run any problem with: `julia --project=. scripts/<problem>.jl`.

- **`docs/`** — Documenter.jl site summarising all experiments and findings (Home, Methodology, a
  page per problem, Findings). Build with `julia --project=docs docs/make.jl` *after* running the
  scripts — `make.jl` copies `plots/*.png` into `docs/src/figures/` (git-ignored) and embeds them.
  Docs depend only on `Documenter` (figures are pre-generated, not built via `@example`).
- **CI (`.github/workflows/`)** — split by concern:
  - `CI.yml` runs **only the tests** (matrix: Julia LTS `1.10` + latest stable `1` × ubuntu/macOS/
    windows; no `arch` pin since macOS runners are arm64). Deps resolve from the registry.
  - `Documenter.yml` **builds and deploys the docs**: `julia-buildpkg` instantiates the main
    project, a step runs all twelve experiment scripts to (re)generate `plots/`, then
    `docs/make.jl` embeds them and deploys. Figures are never committed — regenerated each build.

## Plotting

Three plot types, all sharing the grid layout (one panel per precision):
- `plot_energy_error` — relative energy error vs. time (log y).
- `plot_solution_error` — state error vs. reference over time (log y).
- `plot_solution` — 2D trajectory of each method: phase space `(q, p)` for a one-DOF system,
  configuration space `(q₁, q₂)` otherwise (override via the `coords`/`xlabel`/`ylabel` kwargs).

Conventions:
- **Method groups.** Each plot function takes a `groups` kwarg (default `METHOD_GROUPS`) and writes
  one file per group, appending the group label to `path`. The Hamiltonian problems use three
  groups — `_euler`, `_other`, `_midpoint` (partitioned Gauss(2) variants) — and the Lotka–Volterra
  scripts pass their own single `_variational` group. The order each group is *listed* is the draw/
  legend order. Grouping is orthogonal to the geometric/non-geometric flag, which sets the line
  style (geometric = solid, non-geometric = dashed). Colours are assigned by a method's position in
  the flattened `groups` list (12-colour palette `_PALETTE`), so a given group set colours its
  methods consistently across its figures.
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

Methods compared (Hamiltonian problems):
- Geometric: Symplectic Euler A, Symplectic Euler B, Implicit Midpoint.
- Non-geometric: Explicit Euler, Implicit Euler, Explicit Midpoint, Crank-Nicolson, RK4.
- Partitioned Gauss(2) midpoint variants (`MIDPOINT_METHODS`): `PartitionedTableau(Gauss(2))`,
  `SymplecticPartitionedTableau(Gauss(2))`, and both with the compensation coefficients `â,b̂,ĉ`
  zeroed — all IPRK, built at the run precision via a `tableau(method, T)` accessor.

Degenerate-Lagrangian (Lotka–Volterra) comparison: `Implicit Midpoint`, `VPRK(Gauss(1))`,
`PMVImidpoint`, and (2D only) `CMDVI`. `CMDVI` diverges on the 4D Lagrangian, so LV4d uses only the
first three.

Problems (6): harmonic oscillator, pendulum (both `::Type{T}` `podeproblem`); double pendulum, Toda
lattice, Lotka–Volterra 2D (`LotkaVolterra2dSingular`, `lodeproblem`) and Lotka–Volterra 4D
(`LotkaVolterra4dLagrangian`, `lodeproblem` with `A_quasicanonical_reduced` + `B`) — all
EulerLagrange-generated, no `::Type{T}`, so `make_problem(T)` hand-builds T-typed inputs. Toda
wrinkles: lattice size `N = 16` (not the default 200) and a `hamiltonian(…, N)` closure; its
trajectory plot uses the first lattice site's phase space. The Lotka–Volterra problems are
degenerate Lagrangian (LODE) systems, so only variational / implicit-midpoint methods apply.

Solution-error reference: analytic `exact_solution` for the harmonic oscillator; Float64 `Gauss(8)`
(same grid) for all others.

**Registry-vs-local caveat (resolved).** The registered `GeometricProblems` can lag the local
checkouts, causing signature drift. This bit the pendulum: registered ≤ 0.6.24 had only a
parameter-free `hamiltonian(t, q, p)`, while `energy_error` (via `compute_invariant_error`) passes
`params`, so it errored on CI (registry) but not locally (dev-linked). Resolved by releasing
`GeometricProblems` **0.6.25**, which restores the `(t, q, p, params)` method: `[compat]` now
requires `0.6.25` and the pendulum scripts pass `hamiltonian` directly again (like the others).
Heuristic: if a script errors only on CI but not locally, suspect this kind of signature drift first.

## Verification results

**Precision purity confirmed** — `verify_precision` passes for every successful run across all six
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
- **Float16 double pendulum:** the implicit methods fail with "NaN in direction vector" — a real
  Float16 instability for this stiff, dimensional (g = 9.8), chaotic system. Using the trust-region
  `DogLeg` solver (the default) instead of `Newton` improves robustness generally but does **not**
  rescue these Float16 solves. `run_study` catches the failures per-run so the sweep completes.

All runs are wrapped so a single failure never aborts the study.

## Coarse-step variant findings

- **Harmonic oscillator & pendulum:** the contrast is dramatic at the coarse step — explicit
  Euler and explicit midpoint diverge exponentially (energy → ~1e300 in Float64, clipped at the
  top of the plot), while symplectic Euler A/B stay bounded over the full horizon and implicit
  midpoint / Crank-Nicolson stay near machine level. The energy floor still drops per precision.
  In Float16 some implicit methods fail (time-grid saturation → identical successive times).
- **Double pendulum:** the coarse step is now **Δt = 0.1** (t ≤ 1000), replacing the earlier
  far-too-coarse Δt = 1. The explicit and symplectic-Euler methods blow up quickly (guard-truncated)
  while the implicit midpoint / Crank-Nicolson rules stay bounded much longer; reduced-precision and
  Float16 non-convergence effects are prominent. The short-step (Δt = 0.01) run remains the more
  informative one for this stiff, chaotic problem.
- **Toda lattice:** behaves like the oscillator/pendulum — implicit midpoint / Crank-Nicolson keep
  energy bounded (~1e-5) while explicit midpoint and RK4 drift; the Gauss(8) reference converges
  even at Δt = 1, so the full plot set is produced. In Float16 the two implicit methods fail on the
  long-horizon time-grid saturation; the short scenario runs every method at every precision (the
  bounded bump initial data keeps the exponentials well-behaved).

## Comparison-group findings

- **Partitioned Gauss(2) midpoint variants:** the four algebraically-equivalent forms coincide on
  the linear harmonic oscillator but separate on the nonlinear problems, where the
  symplectic-vs-duplicated tableau construction and the compensation coefficients `â,b̂,ĉ` leave a
  visible imprint on the energy-error fine structure (most so at Float64).
- **Lotka–Volterra variational integrators:** the implicit-midpoint flavours agree to their common
  order and differ in reduced-precision energy behaviour. `CMDVI` integrates the 2D (singular-gauge)
  system fine but **fails on the 4D Lagrangian**: with `A_quasicanonical_reduced` its Jacobian is
  singular, and with other gauge matrices it diverges (~step 356) — so LV4d drops `CMDVI`. LV4d uses
  the quasi-canonical reduced gauge `A_quasicanonical_reduced` + the exact one-form `B` (the default
  `B = 0` makes that `A` singular). At Float16, `VPRK`/`PMVI` hit the usual half-precision
  `NaN`-direction breakdown while `Implicit Midpoint` (and `CMDVI` on 2D) still run.
