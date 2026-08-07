# Methodology

## Precisions

Every experiment is run at four precisions:

```julia
const PRECISIONS = (BFloat16, Float16, Float32, Float64)
```

The order is ascending significand precision, which is also the left-to-right panel order in every
figure. The two 16-bit formats differ in how they spend those bits:

| Type | Significand bits | `eps` | Max finite |
|:--|:--:|:--|:--|
| `BFloat16` ([BFloat16s.jl](https://github.com/JuliaMath/BFloat16s.jl)) | 8 | `2⁻⁷ ≈ 7.8e-3` | `≈ 3.4e38` |
| `Float16` | 11 | `2⁻¹⁰ ≈ 9.8e-4` | `65504` |
| `Float32` | 24 | `≈ 1.2e-7` | `≈ 3.4e38` |
| `Float64` | 53 | `≈ 2.2e-16` | `≈ 1.8e308` |

`BFloat16` keeps `Float32`'s exponent range at the cost of three significand bits, so it is the
*coarsest* of the four for accuracy while being the hardest to overflow. That trade-off is what makes
it interesting here: the exponent range buys nothing for these bounded Hamiltonian problems, while
the lost significand bits cost accuracy directly.

BFloat16s.jl does not yet define every `Base` method the stack needs (`rem` — and therefore
`fld`/`div`/`mod`, on which float-range construction depends — plus `Integer(::BFloat16)` and
construction from `BigInt`). `src/bfloat16_compat.jl` fills those in following BFloat16s' own
`Float32` round-trip pattern.

`NaNMath` has the same gap, and it is the larger one. Its domain-guarded functions — the variants that
return `NaN` where `Base` would throw a `DomainError`, which is what lets an implicit solver reject a
trial state outside a model's domain instead of aborting the run — are defined only for
`Float16`/`Float32`/`Float64`. Since GeometricProblems v0.8.0 passes `nanmath = true` to every
symbolically generated problem, an EulerLagrange-generated vector field reaches the NaNMath variant of
*every* elementary function it contains, so all eleven (`sin`, `cos`, `tan`, `asin`, `acos`, `acosh`,
`atanh`, `log`, `log2`, `log10`, `log1p`) are shimmed, mirroring NaNMath's own guards.

## Methods

The methods are split into a *geometric* (symplectic) group and a *non-geometric* group. The line
style in every plot follows this classification — **solid** for geometric, **dashed** for
non-geometric.

| Method | Group | Order | Notes |
|:--|:--|:--:|:--|
| Symplectic Euler A | geometric | 1 | symplectic partitioned Euler |
| Symplectic Euler B | geometric | 1 | symplectic partitioned Euler |
| Implicit Midpoint  | geometric | 2 | `Gauss(1)` — symplectic, symmetric |
| Implicit Runge-Kutta 4 | geometric | 4 | `Gauss(2)` — symplectic, symmetric |
| Explicit Euler     | non-geometric | 1 | energy-increasing |
| Implicit Euler     | non-geometric | 1 | energy-dissipating |
| Explicit Midpoint  | non-geometric | 2 | |
| Explicit Runge-Kutta 4 | non-geometric | 4 | classical `RK4` |

The last four form the **other methods** plotting group as a 2 × 2: explicit versus implicit at order
2, then at order 4. Within each order the only difference is explicit versus implicit; within each of
those, only the order. The implicit pair are the Gauss collocation rules — `Gauss(1)` *is* the implicit
midpoint rule and keeps that label, while `Gauss(2)` is labelled `Implicit Runge-Kutta 4` to pair it
with the classical explicit fourth-order rule.

Two consequences worth knowing. Using `Gauss(s)` rather than `ImplicitMidpoint()` puts both implicit
rules on the Runge–Kutta code path, which is what lets them share the tableau-driven initial guess
described under [Initial guess](@ref). And on a partitioned problem `Gauss(2)` reduces to the
duplicated partitioned Gauss(2) tableau, so `Implicit Runge-Kutta 4` and `PRK Gauss(2)` below are the
**same integrator** — verified bit-identical at every precision. It appears in both groups
deliberately: as the fourth-order implicit member of the 2 × 2, and as the baseline of the
implementation-variant comparison.

### Partitioned Gauss(2) variants

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

* **Short / fine step** — harmonic oscillator and pendulum use `Δt = 0.1`, `t ≤ 1000`; the Toda
  lattice uses `Δt = 0.1`, `t ≤ 100`; the double pendulum uses `Δt = 0.01`, `t ≤ 10` (its natural
  timescale is much shorter); the Lotka–Volterra problems use `Δt = 0.01`, `t ≤ 10`.
* **Coarse step** — harmonic oscillator and pendulum use `Δt = 1`, `t ≤ 10 000`; the Toda lattice
  uses `Δt = 1` over the same `t ≤ 100` as its short run; the double pendulum uses `Δt = 0.1` over
  the same `t ≤ 10` as its short run; the Lotka–Volterra problems use `Δt = 0.1`, `t ≤ 100`.

Every precision runs the same horizon. The output figure filenames encode the timestep (e.g.
`…_dt_0.1_…`), so the two scenarios of a problem are distinguished by `Δt` rather than by a
"longtime" label.

## Time stepping in a local frame

A reduced-precision *clock* saturates long before a reduced-precision *state* becomes useless. Once
`ulp(t) ≥ Δt` the stored time stops advancing (`t + Δt == t`), and the implicit methods'
`HermiteExtrapolation` initial guess — which divides by `t₁ - t₀` — fails outright with
`t₀ == t₁`. The onset depends on `Δt` and is much earlier for `BFloat16` than for `Float16`
(`capped_final_time` measures it):

| `Δt` | `BFloat16` | `Float16` | `Float32` / `Float64` |
|:--|:--|:--|:--|
| 0.01 | `t ≈ 2` | `t ≈ 16` | no limit at these horizons |
| 0.1 | `t ≈ 16` | `t ≈ 128` | no limit at these horizons |
| 1.0 | `t ≈ 256` | `t ≈ 2048` | no limit at these horizons |

Capping every horizon to those values would leave almost nothing to study. But the limit is an
artefact of *bookkeeping*, not a property of the integrator: none of the methods here derives its
step size from a difference of clock values — they all read `Δt` from the problem — so the absolute
time only ever reaches the vector field, and all six problems are **autonomous**, with `t` an unused
argument in every vector field and Hamiltonian.

`integrate_bounded` therefore advances the solution step in a **fixed local time frame**
(`reset_local!`): after each history shift the stored times are relabelled so that every step looks
like the *second* one — history at `0, Δt`, target at `2Δt`. Consecutive times then stay exactly `Δt`
apart at any precision and any horizon, and Hermite's extrapolation parameter is exactly `s = 2`.

This is a relabelling, not a change of dynamics. Explicit methods (which carry no initial guess)
produce **bit-identical** results either way; implicit methods differ only at solver tolerance
(~`1e-16` in `Float64`, ~`1e-7` in `Float32`). It is also strictly more accurate than a global clock
that has *drifted* without yet saturating: `Float16` resolves only about `0.004` near `t = 5`,
comparable to `Δt = 0.01`, so a differenced interval there is materially wrong. Two details it
depends on: the whole history must be relabelled after `reset!` (which copies stale times along with
the state), and the frame must be anchored at non-negative times, because
`MidpointExtrapolation`'s sub-stepping loop for partitioned problems tests `abs(t + extrap.Δt) <
abs(sol.t)`.

Pass `localclock = false` to `run_study` to step along the problem's own `T`-typed grid instead; it
is retained only to demonstrate the failure. A non-autonomous problem would need it, together with a
horizon capped by `capped_final_time`.

## Initial guess

The implicit solves are seeded by extrapolating from the two previous steps. That extrapolation is the
other place where clock precision leaks into the result, and the fix there is sharper than the local
time frame — it removes the clock from the guess altogether.

Upstream, the implicit partitioned Runge–Kutta integrator builds an *absolute* time for stage `i`,

```julia
t = history[1].t + timestep(int) * tableau(int).q.c[i]
```

and hands it to `HermiteExtrapolation`, which immediately differences it back down to a normalised
one (`Δt = t₁ - t₀`, `s = (tᵢ - t₀) / Δt`). Every quantity in that round trip is stored in the working
precision, so `c → t → c` loses accuracy in half precision — and once the clock saturates it fails
outright with `t₀ == t₁`.

But the normalised time is a **tableau constant**: stage `i` sits `c[i]` steps past the previous step
by construction. `NormalizedHermiteExtrapolation` accepts that constant directly and never looks at a
clock, so `src/initial_guess.jl` overrides `initial_guess!` to call `extrapolate!` with
`tableau(int).q.c[i]` instead of routing through `solutionstep!`. The guess then involves no time
arithmetic whatsoever, at any precision.

The effect is measurable: with this guess the partitioned Runge–Kutta methods become entirely
**clock-independent**. At `BFloat16`, `Δt = 0.1`, `t ≤ 100` — well past the `t ≈ 16` saturation point —
`localclock = true` and `localclock = false` produce *bit-identical* results.

The two mechanisms have complementary jobs, and it is worth being precise about which does what. On
the four **Hamiltonian** problems every method is either explicit (and so has no initial guess to
spoil) or a partitioned Runge–Kutta method whose guess is clock-free — so there the tableau-driven
guess alone suffices and the local frame is redundant. The local frame is what carries the
**Lotka–Volterra** (IODE/LODE) integrators, whose own `initial_guess!` methods difference absolute
times. It is applied unconditionally because it is free, and because it protects any method added
later that has not been given a clock-free guess.

Two implementation notes. `NormalizedHermiteExtrapolation` expects derivative samples scaled by `Δt`
and returns a `Δt`-scaled derivative, while the integrator cache holds unscaled vector-field values, so
the scaling is applied on the way in and undone on the way out. The two formulations agree exactly in
exact arithmetic, and the extrapolated state agrees to round-off at every precision; only the
extrapolated *derivative* rounds differently, because the Hermite derivative basis differences two
terms of size `6c(1+c)·x/Δt` and both formulations lose accuracy to that cancellation in half
precision. That is harmless — an initial guess can only change how many iterations a step needs, never
what it converges to. And the override necessarily applies to *every* implicit partitioned RK
method rather than only the two Gauss rules: `initmethod` collapses them all into the same `IPRK`
integrator type, so the originating method cannot be recovered at dispatch time. That is the coherent
choice regardless, since `Gauss(2)` and `PRK Gauss(2)` are the same integrator here and would otherwise
be seeded differently.

`NormalizedHermiteExtrapolation` requires `GeometricIntegratorsBase` ≥ v0.4.2; the package's
`[compat]` bound for that dependency is tighter (v0.5.1), set by the next section.

## Solver tolerances

The same "is this quantity really at the working precision?" question applies to the nonlinear
solve's *convergence criterion*, which is assessed as

```
rfₐ ≤ f_abstol + f_reltol · ‖F(x₀)‖
```

`run_study` and `integrate_bounded` pass no tolerances of their own; every term above has to be
scaled by the stack, and each is. `SimpleSolvers` builds its `Options` from the state's element type,
giving `x_abstol = 2eps(T)` and `f_reltol = √eps(T)`, and
`GeometricIntegratorsBase.default_options` supplies the absolute term as

```
f_abstol = max(8, solversize(method, problem)) · eps(datatype(problem))
```

scaled to the working precision *and* to the size of the stage system. Caller options are **merged
into** that bundle rather than replacing it, so overriding one keyword does not silently drop
`min_iterations = 1`.

Both factors are load-bearing, and for different reasons.

**The precision factor.** The better the initial guess, the smaller `‖F(x₀)‖` and the more the
*absolute* term alone decides the test — so an improved guess makes convergence *harder* to certify,
which is the wrong way round and is what makes a mis-scaled floor bite here in particular. A
half-precision residual bottoms out near `eps(T)`, so an absolute floor fixed at
`8eps(Float64) ≈ 1.8e-15` is simply unreachable: every implicit solve then runs to its iteration
limit, on every step, for a result that was already as good as the arithmetic allows. The partitioned
problems of this study report `solversize = 0` and so sit on the `max(8, …)` floor at `8eps(T)`,
which is the value the argument above asks for at each precision.

**The size factor.** The `Gauss(8)` **reference** integrations are `Float64`, so `8eps()` is the
"right" value for them by the same argument — and for four of the five it is: the pendulum (16
unknowns), the double pendulum (32), the Toda lattice (**256**) and the 2D Lotka–Volterra system (16)
all reach it with nothing capped. The 4D Lotka–Volterra system does not. It is the degenerate one,
posed with the quasi-canonical reduced gauge matrix, and its residual bottoms out just above that
threshold — so the obstruction is *conditioning*, not dimension. The size factor nevertheless covers
it: at 32 unknowns it gives `32eps ≈ 7.1e-15`, roughly a 4× loosening, and the reference converges in
a handful of iterations. Scaling by the stage size is dimensionally the right thing anyway, the
residual being an `l2` norm over `n` components whose own round-off grows with `n`.

The failure mode this avoids is worth naming, because it is silent: a non-convergent solve is only a
warning, so an unreachable floor lets a *reference* solution be used without having met its
convergence criterion.

`SimpleSolvers` supplies the complementary exit. `assess_convergence` cannot certify a solve that has
converged as far as its arithmetic allows but cannot say so, so a solve that stops making progress is
stopped instead: after `max_stalls = 2` stalled steps it gives up with one actionable message. A line
search is bounded by its own `linesearch_max_iterations` rather than by `max_iterations`, and shares
its solver's `Options` — which is what lets `verbosity = 0` reach it.

`solveropts` is the escape hatch: any `SimpleSolvers.Options` keywords, applied to every run or
resolved per precision by a callable `T -> NamedTuple`. `(verbosity = 0,)` is the useful one here,
silencing a configuration that is *expected* to fail.

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
* **Solution error** — the Euclidean norm of the state difference against the reference solution.
  For the problems without an analytic solution (all but the harmonic oscillator), the coarse-step
  scenarios compute the Float64 `Gauss(8)` reference at the **fine** step (the short scenario's
  `Δt`) and subsample it onto the coarse output grid, so the reference is accurate independent of
  the coarse step. (The short scenarios' reference already uses the fine step, being run at it.)

## Plotting conventions

* Every figure has **one panel per precision** (BFloat16 / Float16 / Float32 / Float64, coarsest
  first) and is produced once per method group. The four Hamiltonian problems use three groups —
  **Euler** (`_euler`), **other** (`_other`), and the **partitioned Gauss(2)** variants (`_gauss2`);
  the Lotka–Volterra problems use a single **variational** group (`_variational`). Scripts pass their
  group set to the plotting routines, which colour methods consistently within it.
* Error plots use a **shared logarithmic y-axis** across all four panels, with the upper limit
  **capped at `1e5`** so runaway (non-geometric) errors are clipped rather than dominating the
  scale; the x-axis is fitted **exactly** to the integration interval.
* Solution plots show the **2D trajectory** of each method (phase space `(q, p)`, or configuration
  space for multi-degree-of-freedom systems) with the **reference** drawn as a black backdrop.
* Colours are consistent per method; the legend is a horizontal row beneath the panels.
