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
    `ALL_METHODS` (12) = `GEOMETRIC_METHODS` (3) + `NONGEOMETRIC_METHODS` (5) + `GAUSS2_METHODS`
    (4 partitioned Gauss(2) variants), grouped for plotting by `METHOD_GROUPS`
    (`euler` / `other` / `gauss2`). The degenerate-Lagrangian Lotka–Volterra comparison uses
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
- **`scripts/{…}_longtime.jl`** — coarse-step drivers: HO/pendulum Δt = 1, t ≤ 10 000; Toda Δt = 1,
  t ≤ 100 (same horizon as its short run); double pendulum Δt = 0.1, t ≤ 10 (same horizon as its
  short run); Lotka–Volterra Δt = 0.1, t ≤ 100. (The double-pendulum coarse run uses the line-search
  `Newton` solver with a `Backtracking` linesearch and `max_iterations = 100`, not the default
  `DogLeg`; the short run matches.)
- **`plots/`** — generated figures. Each study writes, per method group, an energy-error, a
  solution-error, and a 2D solution (trajectory) figure. **Filenames encode the timestep**
  (`…_dt_<Δt>_<group>.png`), so a problem's two scenarios are distinguished by `Δt` (not a
  "longtime" label). Every plot title also carries the run parameters, e.g. `… (Δt = 0.1, t ≤ 100)`.

Run any problem with: `julia --project=. scripts/<problem>.jl`, or regenerate every figure at once
with `julia --project=. scripts/run_all.jl` (it auto-discovers `scripts/*.jl`, activates and
instantiates the project, and `include`s each script into its own throwaway module in a *single*
Julia session — so the shared packages compile once; per-script failures are collected and reported
at the end with a non-zero exit).

- **`docs/`** — Documenter.jl site summarising all experiments and findings (Home, Methodology, a
  page per problem, Findings). Build with `julia --project=docs docs/make.jl` *after* running the
  scripts — `make.jl` copies `plots/*.png` into `docs/src/figures/` (git-ignored) and embeds them.
  Docs depend only on `Documenter` (figures are pre-generated, not built via `@example`).
- **CI (`.github/workflows/`)** — split by concern:
  - `CI.yml` runs **only the tests** (matrix: Julia LTS `1.10` + latest stable `1` × ubuntu/macOS/
    windows; no `arch` pin since macOS runners are arm64). Deps resolve from the registry.
  - `Documenter.yml` **builds and deploys the docs**: `julia-buildpkg` instantiates the main
    project, a step runs all twelve experiment scripts (via `julia --project=. scripts/run_all.jl`)
    to (re)generate `plots/`, then
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
  groups — `_euler`, `_other`, `_gauss2` (partitioned Gauss(2) variants) — and the Lotka–Volterra
  scripts pass their own single `_variational` group. The order each group is *listed* is the draw/
  legend order. Grouping is orthogonal to the geometric/non-geometric flag, which sets the line
  style (geometric = solid, non-geometric = dashed). Colours are assigned by each method's position
  *within its group* from a high-contrast palette (`_PALETTE`), so every per-group figure uses the
  vivid leading colours — keeping the four near-coincident Gauss(2) variants distinguishable.
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
- Partitioned Gauss(2) variants (`GAUSS2_METHODS`): `PartitionedTableau(Gauss(2))`,
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
for all others. For the non-analytic problems the **coarse** scenarios compute the reference at the
**fine** step (`Δt_ref` = the short scenario's Δt, via a `make_reference(T)` closure) and
`solution_error` subsamples it onto the coarse grid (`round.(Int, range(1, nref; length=nsol))`),
so the reference is trustworthy independent of the coarse step.

**Registry-vs-local caveat (resolved).** The registered `GeometricProblems` can lag the local
checkouts, causing signature drift. This bit the pendulum: registered ≤ 0.6.24 had only a
parameter-free `hamiltonian(t, q, p)`, while `energy_error` (via `compute_invariant_error`) passes
`params`, so it errored on CI (registry) but not locally (dev-linked). Resolved by releasing
`GeometricProblems` **0.6.25**, which restores the `(t, q, p, params)` method: `[compat]` now
requires `0.6.25` and the pendulum scripts pass `hamiltonian` directly again (like the others).
Heuristic: if a script errors only on CI but not locally, suspect this kind of signature drift first.

**GeometricProblems v0.7.0 migration.** `[compat]` bumped to `0.7` (no dependency-bound changes;
0.7.0 only altered the problem modules' own API). Two breaking changes touched this repo:
1. The `podeproblem(::Type{T})` / `lodeproblem(::Type{T})` precision constructor was removed — the
   first positional argument is now the initial condition, not the datatype. The harmonic-oscillator
   and pendulum scripts (and `test/runtests.jl`) now hand-build `T`-typed initial conditions from the
   module defaults, e.g. `podeproblem(T.(HO.q₀), T.(HO.p₀); …)` — matching the pattern the double
   pendulum / Toda / Lotka–Volterra scripts already used.
2. `default_parameters` changed from a `const` NamedTuple to a type-parameterized function
   `default_parameters(::Type{T}=Float64)`. Every `map(T, MOD.default_parameters)` became
   `MOD.default_parameters(T)` (double pendulum, Toda, both Lotka–Volterra, and their longtime
   variants). All 167 tests pass and all 12 scripts regenerate their figures cleanly under 0.7.0.

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

- **Float16 + long horizon (time-grid saturation).** Float16 cannot resolve successive time stamps
  once `ulp(t) ≥ Δt` (`t + Δt == t`), which makes the implicit methods' Hermite initial guess throw
  `ArgumentError: t₀ and t₁ … identical`. The onset depends on Δt: for Δt = 1 it is `t ≈ 2048` (integers
  stop being exactly representable); for Δt = 0.1 it is already `t ≈ 100`.
  - **`capped_final_time(T, t₁, Δt)`** (in `study.jl`, exported) is **Δt-aware**: it walks the exact
    grid the solution stores (`GeometricSolutions.TimeSeries` backs it with `tbegin:Δt:tend`) and returns
    the largest final time whose grid is strictly increasing in `T` — i.e. the last point before the
    first saturation collision (or `t₁` if the grid stays resolvable, always so for Float32/Float64 at
    these horizons). It also backs off the endpoint when the rebuilt `StepRangeLen` would pin a colliding
    last point. Resulting Float16 caps: **128.0** for the short scripts (Δt = 0.1) and **2048.0** for the
    longtime scripts (Δt = 1); both grids are verified collision-free.
  - Applied in **all four** HO/pendulum scripts (`harmonic_oscillator.jl`, `pendulum.jl`, and their
    `*_longtime.jl`): **every implicit method now runs at Float16** (Implicit Midpoint / Euler /
    Crank-Nicolson / all four Gauss(2) variants — previously all `[skip]` with the `t₀ == t₁` error),
    while Float32/Float64 keep their full horizon (t ≤ 1000 short, t ≤ 10_000 longtime). Explicit
    Euler/midpoint still diverge — genuine instability, not saturation (guard-truncated).
  - Support machinery: `solution_error` derives the reference-refinement factor from the two **timesteps**
    (not the length ratio), so a capped Float16 run (shorter horizon) is compared against the matching
    leading portion of the full-horizon reference; `_plot_grid` uses **per-panel** x-limits so the shorter
    Float16 panel is not stretched to the longer precisions' horizon.
  - **Upstream `reset!` fix (GeometricIntegratorsBase 0.4.0, commit `b60654d`) — verified.** The step
    reset changed from *accumulating* (`solstep.t += Δt`, done inside `integrate!(solstep,int)`) to
    *setting* the canonical grid time (`reset!(solstep, timesteps(sol)[n])`, done in the outer loop).
    `study.jl`'s `integrate_bounded` loop mirrors this (`reset!(solstep, timesteps(sol)[n])` before
    `integrate!`) — it REQUIRES 0.4.0 (on 0.3.x `reset!` took a Δt and would mis-advance). Verified via
    `scripts/experiments/verify_reset_fix.jl`: the fix is correct and all 167 tests + scripts pass under
    0.4.0, **but it does NOT remove the need for the Float16 cap**. The long-horizon failure is a Float16
    *representability* limit, not an accumulation artifact: the canonical grid `0:1:10000` collected in
    Float16 itself has 5678 collisions (first at t = 2048, where consecutive integers stop being
    representable), so at the full t ≤ 10000 horizon Implicit Midpoint / Crank-Nicolson still throw
    `t₀ and t₁ … identical`. Capped at t ≤ 2000 the grid has 0 collisions and everything runs. So the
    fix (correct clock tracking) and `capped_final_time` (Float16 ≤ 2000) are complementary; both kept.
  - **Upstream deps released; committed against the registry.** GeometricIntegratorsBase 0.4.0 and
    GeometricIntegrators 0.16.5+ (registry resolves 0.16.6) are released; `[compat]` is `GIB = "0.4"`,
    `GI = "0.16.5"`. All `[sources]` dev-links removed and the transient `NonlinearIntegrators` dep (added
    by `Pkg.develop`, never actually used) dropped. `study.jl`'s `integrate_bounded` loop uses the 0.4.0
    API `reset!(solstep, timesteps(sol)[n])` (set canonical grid time; the old `integrate!(solstep,int)`
    internal `t += Δt` reset is gone in 0.4.0). Bump GIB compat past 0.4 only after re-checking that loop.
- **Float16 double pendulum:** the implicit methods fail with "NaN in direction vector" — a real
  Float16 instability for this stiff, dimensional (g = 9.8), chaotic system. Using the trust-region
  `DogLeg` solver (the default) instead of `Newton` improves robustness generally but does **not**
  rescue these Float16 solves. `run_study` catches the failures per-run so the sweep completes.
- **Initial guess: MidpointExtrapolation vs the default HermiteExtrapolation.** `integrate_bounded` /
  `run_study` now take an `initialguess` kwarg (default `nothing` → the method default, which is
  `HermiteExtrapolation()` for every implicit RK/variational method). Investigated whether
  `MidpointExtrapolation()` improves convergence, esp. in Float16 (see
  `scripts/experiments/iguess_extrapolation.jl`, which instruments per-step nonlinear-iteration
  counts via `SimpleSolvers.iteration_number(solverstate(int))`). Findings:
  - **Float16 double pendulum (t ≤ 10, no time-grid saturation): genuine win.** Midpoint *rescues*
    `Implicit Midpoint` and `Implicit Euler`, which complete all 1000 steps where Hermite throws the
    "NaN detected" solver exception. (`Crank-Nicolson` and `PRK Gauss(2)` still fail under both.) So
    switching the **initial guess** does what switching the **solver** (DogLeg) could not.
  - **Float16 long-horizon HO/pendulum (t ≤ 1000): no help.** There the failure is the *saturated
    Float16 clock* (`t₀ == t₁`), not the guess: Hermite errors with "t₀ and t₁ … identical", Midpoint
    just swaps that for a NaN in the solve. Neither completes — the horizon still has to be capped.
  - **Float32/Float64: neutral-to-worse; do NOT adopt globally.** For single-stage methods it is
    comparable; for the multi-stage `PRK Gauss(2)` it is a much worse guess (iteration count explodes,
    nearly every step hits the cap), and it can break `Crank-Nicolson` outright (NaN even at Float64).
  - **Conclusion:** keep `HermiteExtrapolation` as the default; reach for `MidpointExtrapolation` as a
    targeted, per-precision option for stiff Float16 solves (double pendulum `Implicit Midpoint`/`Euler`)
    where it converts a hard failure into a usable run.
  - **Applied:** `run_study`'s `initialguess` kwarg now also accepts a callable `T -> iguess|nothing`
    resolved per precision. Both double-pendulum scripts pass
    `initialguess = T -> T === Float16 ? MidpointExtrapolation() : nothing`, so `Implicit Midpoint`
    and `Implicit Euler` now complete at Float16 (were `[skip] … NaN detected`) while Float32/Float64
    keep the Hermite default (no regression to the multi-stage methods). `Crank-Nicolson` and the four
    Gauss(2) variants still fail at Float16 under both guesses (genuine half-precision breakdown).
    `MidpointExtrapolation`/`HermiteExtrapolation` are re-exported from `ReducedPrecision`.

All runs are wrapped so a single failure never aborts the study.

## Coarse-step variant findings

- **Harmonic oscillator & pendulum:** the contrast is dramatic at the coarse step — explicit
  Euler and explicit midpoint diverge exponentially (energy → ~1e300 in Float64, clipped at the
  top of the plot), while symplectic Euler A/B stay bounded over the full horizon and implicit
  midpoint / Crank-Nicolson stay near machine level. The energy floor still drops per precision.
  In Float16 some implicit methods fail (time-grid saturation → identical successive times).
- **Double pendulum:** both scenarios share the horizon **t ≤ 10** and differ only in Δt (short
  Δt = 0.01, coarse Δt = 0.1). The explicit and symplectic-Euler methods blow up quickly
  (guard-truncated) while the implicit midpoint / Crank-Nicolson rules stay bounded longer;
  reduced-precision and Float16 non-convergence effects are prominent. The short-step (Δt = 0.01)
  run remains the more informative one for this stiff, chaotic problem. Both runs use the
  line-search `Newton` solver with a `Backtracking` linesearch and `max_iterations = 100`.
- **Toda lattice:** behaves like the oscillator/pendulum — implicit midpoint / Crank-Nicolson keep
  energy bounded (~1e-5) while explicit methods diverge early (guard-truncated) and RK4 drifts
  mildly; the Gauss(8) reference converges even at Δt = 1, so the full plot set is produced. Both
  scenarios now share the t ≤ 100 horizon (short Δt = 0.1, coarse Δt = 1), so the Float16 implicit
  methods no longer hit the long-horizon time-grid saturation — every method runs at every precision
  (the bounded bump initial data keeps the exponentials well-behaved).

## Comparison-group findings

- **Partitioned Gauss(2) variants:** the four algebraically-equivalent forms coincide on
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
