# Changelog

All notable changes to ReducedPrecision.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Categories: **Added**; **Changed**; **Removed**; **Bug fixes** = code defects; **Findings** =
results the experiments establish, where a change moved or superseded one; **Documentation**;
**Tests**; **Repository hygiene**.

## [Unreleased] — targeting 0.3.0

Update to GeometricIntegrators 0.18.4, GeometricIntegratorsBase 0.6.4, RungeKutta 0.6.1 and
SimpleSolvers 0.13.2 (CompactBasisFunctions 0.3.1 and QuadratureRules 0.2.1 indirect;
GenericLinearAlgebra, FastTransforms, FFTW, DSP and MKL leave the dependency graph with them).
**No source change was needed**, but the energy-error floors move, so the figures the documentation
quotes are re-measured under **Findings** below.

### Changed
- **Dependency bounds**: `GeometricIntegrators` `0.17` → `0.18`, `GeometricIntegratorsBase` `0.5.1`
  → `0.6.3`, `RungeKutta` `0.5` → `0.6`, `SimpleSolvers` `0.10` → `0.12.1, 0.13`.

  The `SimpleSolvers` entry is a *range* rather than `0.13`, matching what GeometricIntegrators and
  GeometricIntegratorsBase declare: SimpleSolvers 0.13.1 raises its own Julia floor to 1.11, so on
  the 1.10 LTS the resolver has to fall back to 0.12.2. Pinning `0.13` alone would make this package
  uninstallable on the LTS, which is why `julia = "1.10"` is unchanged.
- **`SymplecticEulerA`, `SymplecticEulerB` and the Lotka–Volterra `ImplicitMidpoint` now name
  GeometricIntegratorsBase's methods unambiguously.** GeometricIntegrators 0.18.0 renamed its own
  four Runge–Kutta types to `SymplecticEulerARK`, `SymplecticEulerBRK`, `ImplicitMidpointRK` and
  `CrankNicolsonRK` precisely because both packages exported the unsuffixed names, and
  GeometricIntegrators reexports GeometricIntegratorsBase. The three bare names in `src/methods.jl`
  are unchanged and now resolve to one binding each.
- **The whole stack again precompiles on Julia 1.13.** RungeKutta 0.5.23 pulled in
  GenericLinearAlgebra 0.4.0, whose `LinearAlgebra.eigencopy_oftype` method for `UpperHessenberg`
  overwrites the one Julia 1.13 added to LinearAlgebra itself — a hard error during precompilation.
  RungeKutta 0.6 drops that dependency, having moved its quadrature to QuadratureRules. This is also
  why the pre-bump numbers below could not be re-measured on this machine for an A/B: the old
  dependency set does not load on the Julia available here.

### Findings
- **The half-precision and `Float64` energy-error floors move; the truncation-limited ones do not.**
  Mean relative energy error over the second half of the `Δt = 0.1`, `t ≤ 1000` run, before → after:

  | method / problem | BFloat16 | Float16 | Float32 | Float64 |
  |:--|:--|:--|:--|:--|
  | implicit midpoint, harmonic oscillator | `1.1e-2` → **`5.9e-2`** | `8.9e-3` → **`7.9e-3`** | `2.1e-6` → `2.1e-6` | `5.4e-15` → **`1.3e-14`** |
  | implicit midpoint, pendulum | `2.5e-1` → **`3.6e-1`** | `3.7e-2` → `3.7e-2` | `6.2e-4` → `6.3e-4` | `6.2e-4` → `6.2e-4` |
  | `SPRK Gauss(2)`, pendulum | — | — | `2.2e-5` → **`1.9e-5`** | `2.7e-7` → `2.7e-7` |

  The largest movement is a factor 5.4, on the coarsest format. The two figures that do *not* move at
  all are the pendulum's implicit-midpoint `Float32`/`Float64` pair and that rule's `Float64` value —
  exactly the ones the study identifies as truncation- rather than round-off-limited. `Float32`
  values move in the second significant digit at most. So the bump does not disturb any conclusion;
  it shifts the numbers that were always a property of the arithmetic.

  **The cause is the nonlinear solve, not the tableaux.** `Gauss(1)`'s coefficients are ½, 1, ½,
  exact in binary, so RungeKutta 0.6's more accurate quadrature cannot move them. What moved is
  SimpleSolvers: 0.13 makes `LapackLU` the default linear solver, and its `getrf` pivot order differs
  from the old scalar `LU`; 0.11 and 0.12 changed when a solve stops, through the
  `f_stall_window = 50` criterion GeometricIntegratorsBase 0.6 sets and the `all(isfinite, …)`
  direction test that replaces the bare `isnan` one.

  Re-measured with `scripts/experiments/findings_energy_floor.jl` (added, see below). The "before"
  column is the figure the documentation carried, not a re-run: the pre-bump dependency set does not
  load on this machine's Julia, so this is not an A/B of one variable.
- **The failure set is unchanged.** All 12 experiment scripts complete, and precision purity passes
  for every successful run. The per-run failures are the documented ones and only those: the
  degenerate-Lagrangian half precisions on the two Lotka–Volterra problems (`Implicit Midpoint`,
  `VPRK Gauss(1)`, `PMVI Midpoint` and `CMDVI`), plus `Implicit Euler` at `Float16` on the coarse-step
  pendulum, still the only failure among the four Hamiltonian problems.

### Added
- **`scripts/experiments/findings_energy_floor.jl`**, which reproduces the energy-error floors
  quoted in `docs/src/findings.md` — implicit midpoint and the two fourth-order Gauss(2) variants, on
  the harmonic oscillator and the pendulum, at all four precisions, as the mean relative energy
  error over the second half of the `Δt = 0.1`, `t ≤ 1000` run. It exists so that the next
  dependency bump can tell a moved figure from a stale one in one command rather than by re-running
  the full sweep and reading figures.

### Documentation
- The re-measured floors are carried into `docs/src/findings.md`, `docs/src/index.md`,
  `docs/src/harmonic_oscillator.md` and `docs/src/pendulum.md`. The findings table now says which of
  its entries are expected to move with a dependency release and which are not, and names the script
  that reproduces it.
- **`Implicit Euler` at `Float16` on the coarse pendulum now fails with a different message.**
  `docs/src/pendulum.md` said it "throws a `NaN` in the Newton direction"; it throws
  `NonlinearSolverException`, *"non-finite direction₁ vector"*. SimpleSolvers 0.11 widened that guard
  from `isnan` to `all(isfinite, …)`, so an overflowed direction is now caught as well as a NaN one.
  Same failure, at the same place, reported more precisely.

### Tests
- The suite passes unedited: 296 of 296, with `--check-bounds=auto`.

## [0.2.0] — 2026-08-08

Update to GeometricIntegrators 0.17, GeometricIntegratorsBase 0.5.1, SimpleSolvers 0.10.1 and
GeometricProblems 0.8.2 (EulerLagrange 0.5.1 indirect). Two workarounds this package carried are
deleted, those releases having taken both upstream, and a third upstream change makes a wider
`BFloat16` compatibility shim necessary. Measured results are unchanged wherever the generated code
is unchanged.

**Breaking**, hence the minor bump: `solver_tolerances` and `reference_solution` are removed from the
exported API.

### Changed
- **Dependency bounds**: `GeometricIntegrators` `0.16.5` → `0.17`, `GeometricIntegratorsBase`
  `0.4.2` → `0.5.1`, `GeometricProblems` `0.7` → `0.8`, `SimpleSolvers` `0.9` → `0.10`.
- **`default_options(method)` → `default_options(method, problem)`.** GeometricIntegratorsBase 0.5.1
  standardises `default_options`, `solversize` and `nullvectorsize` to a `(method, problem)`
  argument order. `integrate_bounded` and `reference_solution` called the one-argument form.
- **Caller options are merged into `default_options` upstream rather than replacing it**, so
  `integrate_bounded` no longer splices `default_options(...)...` in by hand to preserve
  `min_iterations = 1`.
- **Toda lattice scripts use `hamiltonian(t, q, p, params)` directly.** GeometricProblems 0.8.1
  hand-writes the Toda vector fields and adds the four-argument method that recovers the lattice
  size from the state, so the size-fixing closure is dropped. Verified bit-identical to the
  five-argument form on this lattice.

### Removed
- **`solver_tolerances` and `reference_solution`.** GeometricIntegratorsBase 0.5.1 defaults
  `f_abstol` to `max(8, solversize(method, problem)) · eps(datatype(problem))` — precision- *and*
  size-scaled at once, which is exactly what the two overrides supplied separately. This sweep's
  partitioned problems report `solversize = 0` and so land on the `max(8, …)` floor at exactly
  `8eps(T)`, bit-for-bit what `solver_tolerances(T)` gave; the `Gauss(8)` reference on the
  degenerate 4D Lotka–Volterra system reports 32, giving `32eps ≈ 7.1e-15` and converging without
  capping. Retaining the local overrides would have fought the new default. The scripts call plain
  `integrate(problem, Gauss(8))` for their references, and `solveropts` now defaults to empty.

### Bug fixes
- **`NaNMath.cos(::BFloat16)` and nine siblings.** GeometricProblems 0.8.0 passes `nanmath = true`
  to all 25 of its symbolic-generation call sites, so an EulerLagrange-generated vector field reaches
  the NaNMath variant of *every* elementary function it contains, not only the `log` in the
  Lotka–Volterra one-form ϑ. NaNMath defines its domain-guarded variants for
  `Union{Float16,Float32,Float64}` alone, so the whole `BFloat16` column of both double-pendulum
  studies (24 runs) threw `MethodError`. `src/bfloat16_compat.jl` now shims all eleven —
  `sin cos tan asin acos acosh atanh log log2 log10 log1p` — mirroring NaNMath's own guards and
  delegating to `Base`. `NaNMath.sqrt`/`pow`/`max`/`min` need no shim, being generic over
  `AbstractFloat`/`Real`; `lgamma` is omitted, `Base` having none to delegate to.

### Findings
- **The Lotka–Volterra half-precision breakdown is reported differently, and diagnosed the same.**
  Where the `VPRK`/`PMVI` solves used to abort with `DomainError: log was called with a negative
  real argument`, `nanmath = true` makes the model return `NaN` instead, so the same breakdown
  surfaces one level up as `NaN detected in direction₁ vector!`. The iterate still leaves the
  positive orthant; the same five runs still fail.
- **The upstream gap this project had recorded is closed.** `assess_convergence` had no stagnation
  exit, so a solve that had converged as far as its arithmetic allowed could only stop by exhausting
  `max_iterations`. SimpleSolvers 0.10 stops a stalled solve after `max_stalls = 2` steps with one
  actionable message, bounds a line search by its own `linesearch_max_iterations`, and lets a line
  search share its solver's `Options` — so `verbosity = 0` silences one. Over the full suite that
  takes the solver warnings from **27 871 to 3 796** and the log from **5.4 MB to 0.77 MB**, with
  "Solver took 1000 iterations" going 665 → 0. There is no `muffle`/`QuietLogger` scaffolding in
  this repository to retire; that helper lived in GeometricIntegrators' own tests.
- **Benchmarks rerun over all twelve scripts and compared against the previous resolution.** Total
  wall clock 397.1 s → 377.3 s; the same 18 runs fail and no new ones do; 367 of 422 energy-error
  and 366 of 422 solution-error maxima agree to within 1 %. The differences are confined to the
  double pendulum and the 4D Lotka–Volterra system, the two problems whose vector fields
  EulerLagrange 0.5 generates differently (`simplify = false` by default plus common-subexpression
  elimination, so results can move in the last bit, which chaos and degeneracy amplify at reduced
  precision), plus two `Implicit Euler` `BFloat16` entries where `max_stalls` stops a stalling solve
  earlier. The Toda lattice, whose vector fields were rewritten upstream wholesale, agrees
  everywhere to better than 1 % — a useful cross-check that the codegen change is the cause.

### Documentation
- Comments and prose throughout `src/`, `scripts/`, `test/`, `README.md` and `docs/src/` describe
  the current state rather than the sequence of changes that produced it. The methodology page's
  "Solver tolerances" section is restructured around the two scaling factors and the evidence for
  each, instead of a before/after narrative.
- `README.md` gains a **Development** section: dev-linking the sibling `Geometric*` checkouts, the
  signature-drift heuristic, and what the two CI workflows do.
- `memory.md` became this file.

### Tests
- New `BFloat16 compatibility shims` testset covering all eleven NaNMath functions (value agreement
  with `Float32` and the `NaN` domain guards) and the four `Base` gaps.
- The solver-tolerance testset now pins the *upstream* scaling as a tripwire — `8eps(T)` at every
  precision on the floor, `16eps` for a 16-unknown stage system — and asserts that a `solveropts`
  override reaches the solver and is merged into `default_options` rather than replacing it.

## [0.1.0] — 2026-07-29

Initial development, from the first commit on 2026-07-11 to the BFloat16 work on 2026-07-29. The
package settles into a reusable `make_problem(T)`-driven sweep over methods × precisions, six
problems, a Documenter site, and CI.

### Added
- **The package and its sweep.** `run_study(make_problem; methods, precisions, …)` builds the
  problem once per precision and runs every method against it, catching per-run failures so one
  failure never aborts a sweep. `Run` records the method, precision, problem, solution, error
  message and divergence step.
- **Four precisions.** `Float16`/`Float32`/`Float64` from the start; `BFloat16` on 2026-07-29,
  ordered first in `PRECISIONS` as the coarsest (8 significand bits against `Float16`'s 11).
- **Six problems**: harmonic oscillator, pendulum, double pendulum, Toda lattice (added
  2026-07-11 as the fourth Hamiltonian case, at `N = 16` rather than the default 200 to keep the
  sweep tractable), and the degenerate-Lagrangian Lotka–Volterra 2D (`LotkaVolterra2dSingular`) and
  4D (`LotkaVolterra4dLagrangian`, quasi-canonical reduced gauge `A` plus the exact one-form `B`;
  the default `B = 0` makes that `A` singular). Each has a short-step and a coarse-step scenario.
- **Twelve methods on the Hamiltonian problems**, in three plotting groups. Geometric: symplectic
  Euler A/B, implicit midpoint (`Gauss(1)`), implicit RK4 (`Gauss(2)`). Non-geometric: explicit and
  implicit Euler, explicit midpoint, `RK4`. Plus four partitioned `Gauss(2)` variants
  (`GAUSS2_METHODS`, added 2026-07-11): duplicated versus symplectic-conjugate tableau, each with
  and without the rounding-compensation coefficients `â, b̂, ĉ`, all built at the run precision
  through a `tableau(method, T)` accessor.
- **A separate variational comparison for the Lotka–Volterra problems** (`LV2D_METHODS`,
  `LV4D_METHODS`), on which the explicit, symplectic-Euler and DIRK methods are undefined: implicit
  midpoint, `VPRK(Gauss(1))`, `PMVImidpoint`, and — 2D only — `CMDVI`. The `GaussVPRK` wrapper
  rebuilds `VPRK(Gauss(1))` at the run precision, its own `initmethod` otherwise baking in
  `Float64`.
- **A divergence guard** (2026-07-11). `integrate_bounded` replicates GeometricIntegrators' own
  stepping loop and stops as soon as the state goes non-finite or exceeds `bound` (default `1e3`),
  filling the tail with `NaN` so the metrics and plots break cleanly at the blow-up.
  `Run.diverged` records the step.
- **The type-purity gate.** `assert_precision` / `verify_precision` check that `datatype`,
  `timetype` and the stored `q`/`p` element types all equal the requested precision, for the
  problem and the solution, on every successful run.
- **Local-frame stepping** (2026-07-29). `reset_local!` relabels the solution step's clock after
  each history shift so every step looks like the second one — history at `0, Δt`, target at
  `nhistory · Δt`. Legitimate because the integrators take `Δt` from the problem and never from a
  clock difference, and all six problems are autonomous.
- **A tableau-driven initial guess** (2026-07-29). `src/initial_guess.jl` overrides `initial_guess!`
  for `IPRK` on `PODE`/`HODE` to pass the tableau node `c[i]` straight to `extrapolate!` with
  `NormalizedHermiteExtrapolation`, so no clock value enters the guess at all.
- **`capped_final_time`**, the diagnostic that measures where a `T`-typed grid stops advancing.
- **`src/bfloat16_compat.jl`** (2026-07-29) — the `Base` methods BFloat16s.jl lacks. `rem` is the
  one that cannot be done without: `fld`/`div`/`mod`/`cld`/`divrem` route through it and a float
  range's length is `fld(stop - start, step) + 1`, which is how `GeometricSolutions.TimeSeries`
  builds a time grid, so without it no `BFloat16` problem can be allocated at all. Plus
  `Integer(::BFloat16)` (the other half of that expression), `BFloat16(::BigInt)` (how RungeKutta's
  exact-rational tableau coefficients arrive), `sincos` (which otherwise recurses into itself and
  overflows the stack), `atan`, `fma`, and `NaNMath.log`.
- **Three plotting routines**, all on a shared grid layout of one panel per precision:
  `plot_energy_error`, `plot_solution_error` and `plot_solution` (phase space for a one-degree-of-
  freedom system, configuration space otherwise, overridable through `coords`/`xlabel`/`ylabel`).
  Each writes one file per method group.
- **`scripts/run_all.jl`** (2026-07-12, replacing a shell runner): discovers `scripts/*.jl`
  automatically and includes each into its own throwaway module in a *single* Julia session, so the
  shared packages compile once instead of once per script. Per-script failures are collected and
  reported at the end with a non-zero exit.
- **A Documenter site** (2026-07-11): Home, Methodology, a page per problem, and Findings.
- **CI**, split by concern (2026-07-11): `CI.yml` runs only the tests, on Julia LTS (`1.10`) and
  latest stable (`1`) across Linux, macOS and Windows, reporting coverage to Codecov;
  `Documenter.yml` regenerates `plots/` by running all twelve scripts, then builds and deploys the
  site. Figures are never committed.

### Changed
- **Dependency bounds over the development period**: `GeometricProblems` `0.6.24` → `0.6.25`
  (2026-07-11) → `0.7` (2026-07-13); `GeometricIntegrators` `0.16.4` → `0.16.5` and
  `GeometricIntegratorsBase` `0.3.2` → `0.4` (2026-07-14) → `0.4.2` (2026-07-29, the release adding
  `NormalizedHermiteExtrapolation`); `RungeKutta 0.5` and `SimpleSolvers 0.9` added 2026-07-11;
  `BFloat16s` and `NaNMath` added 2026-07-29.
- **From path dependencies to the registry** (2026-07-11). The initial `[sources]` block pointed at
  sibling checkouts, which meant CI resolved something different from local work. Replaced by
  `[compat]` bounds on registered versions.
- **A single problem form runs the whole method set.** The literal method list cannot: the special
  `ExplicitEuler`/`ImplicitEuler` are ODE-only while `SymplecticEulerA`/`B` are PODE/HODE-only. The
  studies use the partitioned form (`PODE` for oscillator and pendulum, `HODE` for double pendulum
  and Toda lattice) together with the numerically identical RK twins `ExplicitEulerRK` /
  `ImplicitEulerRK`, which auto-promote to partitioned RK.
- **The `other` plotting group became a 2 × 2** (explicit versus implicit at order 2, then at order
  4). Crank–Nicolson was dropped as not fitting it, and `ImplicitMidpoint()` was replaced by
  `Gauss(1)` so both implicit rules sit on the Runge–Kutta code path and share the tableau-driven
  initial guess. On a partitioned problem `Gauss(2)` reduces to `PRK Gauss(2)`: `Implicit
  Runge-Kutta 4` and `PRK Gauss(2)` are the *same* integrator, verified bit-identical at every
  precision.
- **`DogLeg` is the default nonlinear solver** (2026-07-11), the trust-region iteration being more
  robust in reduced precision than a line-search Newton. The double-pendulum scripts opt back into
  `Newton()` with a `Backtracking` line search and `max_iterations = 100`.
- **Figure filenames encode the timestep** (`…_dt_<Δt>_<group>.png`) rather than carrying a
  "longtime" label, so a problem's two scenarios are distinguished by `Δt`. Titles carry the run
  parameters. The `=` character was dropped from filenames to unbreak the Documenter build
  (2026-07-11).
- **Coarse-step scenarios reference a fine-step integration.** For the non-analytic problems the
  coarse run computes its `Gauss(8)` reference at the short scenario's `Δt` through a
  `make_reference(T)` closure, and `solution_error` subsamples it onto the coarse grid, so the
  reference is trustworthy independent of the coarse step. The refinement factor is taken from the
  two timesteps rather than from the length ratio, which is only equal when both grids span the same
  horizon.
- **Both scenarios of the double pendulum and the Toda lattice share a horizon**, differing only in
  `Δt` (2026-07-12).
- **`timevalues` uses the nominal grid `t₀ + n·Δt`** rather than the stored clock, which saturates
  at long horizons in half precision.

### Bug fixes
- **Pendulum energy error on CI but not locally** (2026-07-11). The registered
  `GeometricProblems ≤ 0.6.24` had only a parameter-free `hamiltonian(t, q, p)`, while
  `energy_error` passes `params` through `compute_invariant_error` — so the script failed against
  the registry and passed against the dev-linked checkout. Fixed upstream in `GeometricProblems`
  0.6.25, which restores the `(t, q, p, params)` method; the local workaround is dropped and the
  compat bound raised.
- **Overflow-safe shared y-limits.** Non-geometric energy errors reach ~`1e308` while still finite,
  which overflows Makie's log autolimit padding to `Inf` and dwarfs the scale. All panels of a figure
  share one range computed from the finite data, capped at `1e5`, with runaway lines clipped.
- **Exact x-limits** taken from the problem's `timespan` rather than the accumulated grid endpoint,
  which rounds slightly short or long at low precision and would drop the final tick.

### Findings
- **Precision purity holds throughout.** No library in the stack — GeometricIntegrators,
  GeometricIntegratorsBase, GeometricSolutions, GeometricEquations, GeometricBase, SimpleSolvers —
  silently promotes to `Float64`, for any method, problem or precision, including the hand-built
  half-precision constructions of the double pendulum and Toda lattice.
- **Most of what looked like a half-precision limit was a bookkeeping limit** (2026-07-29). Three
  quantities were being carried or compared at the working precision without needing to be, and each
  produced failures that read as hardware limits:
  1. **The clock.** A `T`-typed time variable stops advancing once `ulp(t) ≥ Δt`, after which the
     Hermite guess fails with `t₀ == t₁`. Onsets, at `Δt = 0.01 / 0.1 / 1.0`: `BFloat16` `t ≈ 2 /
     16 / 256`, `Float16` `t ≈ 16 / 128 / 2048`. Local-frame stepping removes the limit; explicit
     methods are bit-identical to global-clock runs and implicit ones differ only at solver
     tolerance.
  2. **The initial guess.** The stage node is a tableau constant `c[i]`, but upstream reconstructed
     it as `history[1].t + Δt·c[i]` and differenced it back down. Feeding the constant directly
     makes the partitioned RK methods *completely* clock-independent: past the `BFloat16` saturation
     point, `localclock = true` and `false` give bit-identical results.
  3. **The solver tolerance.** `f_abstol` was pinned to `8eps(Float64)` at every precision, which a
     half-precision residual can never satisfy — so every implicit solve burnt all 1000 iterations
     per step, silently. Scaling it to the run's own precision left `Float64` bit-identical and the
     rest unchanged to round-off while running 2.6–40× faster and taking a full `run_all.jl` log
     from **900 MB** to a few MB. The same floor was unreachable for the *`Float64`* `Gauss(8)`
     reference on the 4D Lotka–Volterra system — from conditioning, not precision — which capped out
     on every step for a solution the study treats as ground truth; loosening it by the stage size
     gave a **563× speedup** (36.9 s → 0.066 s over 500 steps, 270 capped solves → 0) with the two
     references agreeing to `1.4e-12`.
- **Consequently the `MidpointExtrapolation` workaround was removed** (2026-07-29). It had been
  introduced on 2026-07-14 to rescue the `Float16` double-pendulum solves from a `NaN` in the Newton
  direction. With the clock out of the way `HermiteExtrapolation` converges for every Gauss-family
  method at every precision on every Hamiltonian problem, and Midpoint is now the *worse* choice —
  five failures at `Float16` and one at `BFloat16` where Hermite has none. The option remains
  available through `run_study(…; initialguess = …)`.
- **The reduced-precision horizon caps were removed** (2026-07-29). `Float16` horizons had been
  capped at `t = 2000` (2026-07-14), then made `Δt`-aware (2026-07-14), to keep the time grid
  resolvable. Local-frame stepping makes them unnecessary; `capped_final_time` remains as the
  diagnostic that measures the onset.
- **What is genuinely precision-limited.** The round-off floor on the achievable energy error — the
  cost the study is actually about. For implicit midpoint at `Δt = 0.1`, across
  `BFloat16`/`Float16`/`Float32`/`Float64`: harmonic oscillator `1.1e-2` / `8.9e-3` / `2.1e-6` /
  `5.4e-15`, pendulum `2.5e-1` / `3.7e-2` / `6.2e-4` / `6.2e-4`. The pendulum's equal `Float32` and
  `Float64` figures show that rule is truncation-limited there, so which floor a curve sits on is
  method-dependent. `BFloat16` is a factor ≈ 8 worse than `Float16`, exactly the ratio of their
  `eps`: its wider exponent range buys nothing for these bounded Hamiltonian systems.
- **The degenerate-Lagrangian variational integrators break down in half precision, and genuinely
  so.** These are singular systems whose one-form involves `log(q)`, so the iterate must stay in the
  positive orthant, and 8–11 significand bits are not enough to keep it there. This is the one
  family whose half-precision failures are not a clock or tolerance artefact, and it is not a clean
  `BFloat16`-versus-`Float16` story: the two formats fail on overlapping but not identical sets.
  `CMDVI` integrates the 2D singular-gauge system but fails on the 4D Lagrangian — its Jacobian is
  singular under `A_quasicanonical_reduced` and it diverges under other gauge matrices — so the 4D
  comparison drops it.
- **Robustness is problem-dependent.** The Toda lattice, whose bump initial data keeps the state
  bounded, runs every method at every precision.
- **The partitioned `Gauss(2)` variants coincide on the linear oscillator and separate on the
  nonlinear problems**, where the symplectic-versus-duplicated tableau construction and the
  compensation coefficients `â, b̂, ĉ` leave a visible imprint on the energy-error fine structure —
  most so at `Float64`, since in half precision the round-off floor swamps the difference.

### Tests
- Unit tests covering the method registry and its plotting groups, `run_study` plus the purity gate
  at every precision, the `capped_final_time` onsets, the tableau-driven guess (including that
  `Gauss(2)` and `PRK Gauss(2)` agree exactly), long horizons past clock saturation, the error
  metrics, per-run failure capture, the divergence guard, and that the plotting routines write one
  figure per group including a script-supplied custom group set.

### Repository hygiene
- `.JuliaFormatter.toml`, `.gitignore` (excluding `Manifest.toml`, `plots/*.png` and
  `docs/src/figures/`), the MIT `LICENSE`, and `.github/dependabot.yml` keeping the GitHub Actions
  current.
- `scripts/experiments/` holds two one-off investigations kept for reference:
  `iguess_extrapolation.jl` (Hermite versus Midpoint initial guess) and `verify_reset_fix.jl`
  (local-frame versus global-clock stepping).

## Open Issues

- **Not every quoted figure has been re-measured since the 0.18 / 0.13 dependency bump.**
  `scripts/experiments/findings_energy_floor.jl` covers the harmonic oscillator and the pendulum. The
  energy- and solution-error levels quoted in `docs/src/double_pendulum.md`,
  `docs/src/toda_lattice.md`, the two Lotka–Volterra pages, and the solution-error levels on
  `docs/src/harmonic_oscillator.md` were not: the sweep re-runs and the figures regenerate, but those
  numbers are carried over. They are stated to one significant figure, and the movement measured on
  the two problems that *were* re-measured would leave most of them standing — but that is an
  argument, not a measurement. Extending the script to the remaining problems is the fix.
- The `1.4e-12` agreement and the `2.6–40×` speed-ups in `docs/src/findings.md` are **historical A/B
  figures**: they compare a fixed solver tolerance against a precision-scaled one, and the unscaled
  variant no longer exists in the stack. They are unaffected by a dependency bump and were not
  re-run.
- **A pre-bump baseline cannot be reproduced on this machine.** RungeKutta 0.5.23 pulls in
  GenericLinearAlgebra 0.4.0, which does not precompile on Julia 1.13, and setting up an older Julia
  to work around it was blocked here. Any A/B against a state before this bump needs a machine with
  Julia 1.11 or 1.12 available.
- A sweep emits some 3800 warnings, dominated by `DogLeg trust-region radius Δ underflowed` and
  Hermite extrapolation's `history[1] and history[2] are identical` — both pre-existing upstream
  messages, present in SimpleSolvers 0.10.1 and GeometricIntegratorsBase 0.5.3 as well. They are not
  suppressed, because `verbosity = 0` would also hide the genuine non-convergence reports the study
  reads; separating the two is open work.
