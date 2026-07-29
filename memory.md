# ReducedPrecision — project summary

Analysis of geometric (symplectic) vs. non-geometric integrators run in varying floating-point
precision (BFloat16, Float16, Float32, Float64) on example problems from GeometricProblems.

## What was built

`ReducedPrecision` is a proper Julia package (`/Users/mkraus/Datashare/Julia/ReducedPrecision`):

- **`Project.toml`** — package name/uuid + deps (`CairoMakie`, `BFloat16s`, `NaNMath` + the Geometric*
  ecosystem), all resolved from the **registry** (no `[sources]`; `[compat]` pins the working versions).
  `GeometricIntegratorsBase` has a `0.4.2` lower bound because that is the release that added
  `NormalizedHermiteExtrapolation`, which `src/initial_guess.jl` needs. For local work against the
  sibling checkouts, dev-link them into the git-ignored Manifest with
  `pkg> dev ../GeometricBase ../GeometricEquations …`. Test-only `Test` via `[extras]`/`[targets]`.
- **`src/`** — reusable pipeline, split into logical units and stitched together by
  `ReducedPrecision.jl` (usings/exports/`include`s only):
  - `bfloat16_compat.jl` — `Base`/`NaNMath` methods BFloat16s.jl v0.6.1 lacks but the stack needs
    (`rem` — hence `fld`/`div`/`mod`, on which float-range and therefore `Solution` construction
    depend — plus `Integer(::BFloat16)`, `BFloat16(::BigInt)`, `sincos`, `NaNMath.log`).
  - `initial_guess.jl` — the tableau-driven, clock-free `initial_guess!` for `IPRK` on PODE/HODE.
  - `methods.jl` — `PRECISIONS` (4, BFloat16 first), `MethodSpec`, the method registries, and the
    plotting groups. `ALL_METHODS` (12) = `GEOMETRIC_METHODS` (4) + `NONGEOMETRIC_METHODS` (4) +
    `GAUSS2_METHODS` (4 partitioned Gauss(2) variants), grouped for plotting by `METHOD_GROUPS`
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
  lotka_volterra_4d}.jl`** — short-step drivers (HO & pendulum: Δt = 0.1, t ≤ 1000; Toda: Δt = 0.1, t ≤ 100; double
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
- Geometric: Symplectic Euler A, Symplectic Euler B, Implicit Midpoint (`Gauss(1)`),
  Implicit Runge-Kutta 4 (`Gauss(2)`).
- Non-geometric: Explicit Euler, Implicit Euler, Explicit Midpoint, Explicit Runge-Kutta 4 (`RK4`).
- The `other` plotting group is these last four as a **2 × 2**: explicit vs implicit at order 2, then
  at order 4 (`Explicit Midpoint`, `Explicit Runge-Kutta 4`, `Implicit Midpoint`,
  `Implicit Runge-Kutta 4`). Crank-Nicolson was dropped when this group was reorganised — it does not
  fit the 2 × 2 — and `ImplicitMidpoint()` was replaced by `Gauss(1)` so both implicit rules sit on
  the Runge-Kutta code path and share the tableau-driven initial guess in `src/initial_guess.jl`.
  Note `Gauss(2)` reduces to `PRK Gauss(2)` on a partitioned problem: `Implicit Runge-Kutta 4` and
  `PRK Gauss(2)` are the **same integrator**, verified bit-identical at every precision.
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
   variants). All 261 tests pass and all 12 scripts regenerate their figures cleanly under 0.7.0.

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

**The headline result of the BFloat16 work (2026-07-29): most of what this section used to call a
"genuine Float16 limit" was not one.** Three separate quantities were being carried or compared at the
working precision when they did not need to be, and each failure they caused looked like a hardware
limit until the quantity was fixed:

1. **The clock.** A `T`-typed time variable stops advancing once `ulp(t) ≥ Δt`, after which the Hermite
   initial guess throws `t₀ == t₁`. Onsets measured by `capped_final_time`: BFloat16 t ≈ 2 / 16 / 256 at
   Δt = 0.01 / 0.1 / 1.0; Float16 t ≈ 16 / 128 / 2048. `integrate_bounded` now advances the solution
   step in a **fixed local time frame** (`reset_local!`): after each history shift the times are
   relabelled so every step looks like the second one (history at `0, Δt`, target at `2Δt`). Legitimate
   because the integrators take Δt from the problem, never from a clock difference, and all six problems
   are autonomous. Verified: explicit methods bit-identical to the old global-clock runs, implicit ones
   differ only at solver tolerance. This retired the `capped_final_time` horizon caps in all four
   HO/pendulum scripts (the function stays as the diagnostic that measures the onset).
2. **The initial guess.** Sharper still: the stage node the Hermite guess extrapolates to is a *tableau
   constant* `c[i]`, but upstream reconstructed it as `history[1].t + Δt·c[i]` and then differenced it
   back down. `src/initial_guess.jl` overrides `initial_guess!` for `IPRK` on PODE/HODE to pass `c[i]`
   straight to `extrapolate!` with `NormalizedHermiteExtrapolation`, so no clock value enters the guess
   at all. Result: the partitioned RK methods are now **completely clock-independent** — at BFloat16
   past the saturation point, `localclock = true` and `false` give bit-identical results. Needs
   GeometricIntegratorsBase ≥ 0.4.2, the release that added `NormalizedHermiteExtrapolation`.
3. **The solver tolerance.** `GeometricIntegratorsBase.default_options` pins `f_abstol = 8eps()` —
   `8eps(Float64)` — at every precision. Since convergence is `rfₐ ≤ f_abstol + f_reltol·‖F(x₀)‖`, and a
   good guess makes `‖F(x₀)‖` small, a half-precision residual could never satisfy it: every implicit
   solve burnt all 1000 iterations *per step*, silently. `run_study` now defaults to
   `solver_tolerances(T) = (f_abstol = 8eps(T),)`. Float64 bit-identical, others unchanged to round-off,
   2.6–40× faster, and `run_all.jl`'s log dropped from **900 MB** to a few MB.
4. **The same tolerance in the Float64 references.** The `Gauss(8)` reference for the 4D Lotka–Volterra
   system could not reach `8eps(Float64)` either — not from precision but from *conditioning* (it is the
   degenerate one, with `A_quasicanonical_reduced`) — and so capped out on every step for a solution the
   study treats as ground truth. `reference_solution(problem, method)` scales the floor by `√n` for `n`
   stage unknowns (dimensionally right: the residual is an `l2` norm over `n` components). **563×
   faster** on LV4D (36.9 s → 0.066 s for 500 steps, 270 capped → 0), references agree to 1.4e-12. All
   ten `Gauss(8)` call sites in `scripts/` now use it. Probed the rest: pendulum (n=16), double pendulum
   (n=32), Toda (n=**256**) and LV2D (n=16) all reached `8eps` already, so this is not a size effect —
   `√n` is simply a principled floor that happens to supply the 5.7× LV4D needs.
   Upstream gap worth reporting: `assess_convergence` has **no stagnation exit** — every path requires
   the residual test, and `f_settled = rfₛ ≤ ‖F‖·f_suctol` cannot fire once `F` jitters at its own floor.

**Consequence for the method matrix:** with (1)–(3) in place, `HermiteExtrapolation` converges for every
Gauss-family method at every precision on every Hamiltonian problem/scenario — so the old
`MidpointExtrapolation` workaround in the double-pendulum scripts is **removed** (Midpoint is now
*worse*: it causes 5 failures at Float16 and 1 at BFloat16 where Hermite has none). The old "Float16
double pendulum implicit methods fail with NaN in direction" finding is superseded: its real cause was
the differenced Hermite interval, since Float16 resolves only ~0.004 near t = 5, comparable to Δt = 0.01.

**What remains genuinely precision-limited:**

- **The round-off floor on the achievable energy error** — the honest, unavoidable cost, and the thing
  the study is actually about. Measured for implicit midpoint at Δt = 0.1: HO `1.1e-2` / `8.9e-3` /
  `2.1e-6` / `5.4e-15` and pendulum `2.5e-1` / `3.7e-2` / `6.2e-4` / `6.2e-4` across
  BFloat16 / Float16 / Float32 / Float64. Note the pendulum's Float32 = Float64: that rule is
  truncation-limited there, so which floor a curve sits on is method-dependent.
- **BFloat16 is a factor ≈ 8 worse than Float16** — exactly the ratio of their `eps`. Its wider exponent
  range is useless for these bounded Hamiltonian systems, so it is simply the coarsest of the four.
- **The degenerate-Lagrangian (Lotka–Volterra) variational integrators break down at BFloat16.**
  `VPRK Gauss(1)`, `PMVI Midpoint` and `CMDVI` all throw "NaN detected in direction₁ vector"; the local
  frame lets `Implicit Midpoint` get further (diverging around step 624 at Δt = 0.1) but does not rescue
  them. These are singular systems with a `log(q)` one-form, so 8 significand bits is genuinely too
  coarse. This path is not covered by the tableau-driven guess, which targets PODE/HODE only.
- **Problem-dependent robustness.** The Toda lattice, whose bump initial data keeps the state bounded,
  runs every method at every precision.

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
