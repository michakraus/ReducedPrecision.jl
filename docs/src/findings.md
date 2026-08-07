# Findings

## Geometric vs. non-geometric integrators

Across the four Hamiltonian problems and both scenarios the qualitative distinction is consistent
and matches the theory of geometric numerical integration (the two degenerate-Lagrangian
Lotka–Volterra problems use a separate variational-integrator comparison, discussed below):

* **Symplectic methods** (symplectic Euler A/B, implicit midpoint) keep the energy error
  **bounded** — it oscillates around a small value rather than growing — over arbitrarily long
  integrations. In phase space their trajectories stay on a closed level set of a nearby modified
  Hamiltonian.
* **Explicit Euler** systematically **increases** the energy (spirals outward in phase space);
  **implicit Euler** systematically **dissipates** it (spirals inward). Both are unusable for
  long-time dynamics regardless of precision.
* The **Gauss collocation rules** — implicit midpoint (`Gauss(1)`, order 2) and implicit RK4
  (`Gauss(2)`, order 4) — being symplectic *and* symmetric, nearly conserve energy, with a bounded
  error set by the precision.
* **Explicit RK4** is accurate over short times but, being non-symplectic, exhibits a slow energy
  drift over long horizons. Set against implicit RK4 at the same order, this isolates the effect of
  symplecticity from the effect of accuracy: the two track each other closely at first, then separate
  once the drift accumulates, while their *truncation* errors remain comparable.

## The role of precision

* For the energy-conserving methods the **precision sets the floor** of the (bounded) energy error.
  Measured for the implicit midpoint rule on the harmonic oscillator at `Δt = 0.1` over `t ≤ 1000`
  (mean over the second half of the run):

  | | BFloat16 | Float16 | Float32 | Float64 |
  |:--|:--|:--|:--|:--|
  | harmonic oscillator | `1.1e-2` | `8.9e-3` | `2.1e-6` | `5.4e-15` |
  | pendulum | `2.5e-1` | `3.7e-2` | `6.2e-4` | `6.2e-4` |

  The other problems show the same ordering.
* **The floor is not always round-off.** On the pendulum the implicit midpoint rule reaches the *same*
  `6.2e-4` at Float32 and Float64: there the error is set by the method's truncation error at
  `Δt = 0.1`, not by the arithmetic, so extra precision buys nothing. The `SPRK Gauss(2)` rule, being
  of higher order, does keep improving (`2.2e-5` → `2.7e-7`). Reading a precision study therefore
  requires knowing which of the two floors a given curve is sitting on.
* For the low-order methods the **solution error** is likewise often dominated by truncation rather
  than round-off, so Float32 and Float64 solution errors can be nearly identical while only the half
  precisions show a round-off floor.
* Reducing precision therefore mainly matters for (a) the achievable energy-conservation floor and
  (b) the robustness of the implicit solves and the time grid — not so much for the truncation-
  limited accuracy of the low-order methods.

## BFloat16 versus Float16: exponent range is the wrong trade

`BFloat16` and `Float16` are both 16 bits wide and differ only in how they split them. `BFloat16`
carries `Float32`'s 8-bit exponent, so it reaches `≈ 3.4e38` instead of `65504`, and pays for that
with three significand bits (`eps` `2⁻⁷ ≈ 7.8e-3` against `2⁻¹⁰ ≈ 9.8e-4`).

For the problems studied here that is a bad trade in both directions:

* **The extra range is never used.** These are bounded Hamiltonian systems; the state stays `O(1)`
  and the only thing that ever reaches `Float16`'s `65504` ceiling is a method that has already
  blown up (explicit Euler at a coarse step), which the divergence guard stops anyway. `BFloat16`
  merely lets such a run climb further before being cut off — visible as the explicit-Euler curves
  running higher up the clipped `1e5` ceiling.
* **The missing significand bits cost accuracy directly.** `BFloat16`'s bounded energy-error floor
  sits a factor of ≈ 8 above `Float16`'s — precisely the ratio of their `eps` — and a `BFloat16`
  time grid saturates eight times sooner.

So `BFloat16` is the weakest of the four precisions here, and predictably so. This is not an argument
against the format — its range is what makes it useful for the wide dynamic ranges of neural-network
training — but it does mean that a format chosen for machine learning is not automatically a good fit
for geometric integration, where accuracy per bit is what matters.

## Type purity

`verify_precision` passes for **every successful run across all six problems and all four
precisions**: `datatype`, `timetype`, and the element types of the stored `q`/`p` arrays all equal
the requested precision. No library in the stack (`GeometricIntegrators`, `GeometricIntegratorsBase`,
`GeometricSolutions`, `GeometricEquations`, `GeometricBase`, `SimpleSolvers`) silently promotes to
`Float64`, including for the hand-built half-precision constructions of the double pendulum and Toda
lattice. The key requirement is that **both** the initial conditions **and** the timespan/timestep
are created at the target precision, so that the data type and the time type agree.

`BFloat16` is the sharpest test of this, because its arithmetic is *implemented* by widening to
`Float32` and rounding back. The gate confirms that the widening never escapes into the stored state:
every `q`, `p` and time value comes back as `BFloat16`. The generic (non-BLAS) linear algebra the
implicit solvers need — `lu`, `\`, `norm`, `dot` — works on `BFloat16` unchanged.

## The reduced-precision clock, and why it is not a real limit

The most conspicuous half-precision failure mode in this study is an artefact of bookkeeping rather
than a property of the integrators.

A time variable stored in `T` stops advancing once `ulp(t) ≥ Δt`: successive stamps round to the same
value, and the implicit methods' `HermiteExtrapolation` — which divides by `t₁ - t₀` — fails with
`t₀ == t₁`. The onset arrives early, and eight times earlier for `BFloat16` than for `Float16`
(`t ≈ 16` vs `t ≈ 128` at `Δt = 0.1`). Capping every horizon there would leave almost nothing to
study.

But no method here derives its step size from a difference of clock values — they all read `Δt` from
the problem — so absolute time only ever reaches the vector field, and all six problems are
autonomous. Advancing each step in a [local time frame](@ref "Time stepping in a local frame") is
therefore exact, and it removes the limit entirely: every precision runs every horizon. Against
stepping along the problem's own grid (`localclock = false`) the explicit methods are bit-identical;
the implicit ones differ only at solver tolerance.

The same reasoning applied once more gives a sharper result still. The upstream *initial guess*
reaches its interpolation node by building an absolute stage time and differencing it back down —
even though that node is a tableau constant (`c[i]`) that needs no clock at all. Feeding the constant
directly to `NormalizedHermiteExtrapolation` (see [Initial guess](@ref)) makes the partitioned
Runge–Kutta methods entirely clock-independent: at `BFloat16` past the saturation point, stepping along
the saturating grid and stepping in the local frame give *bit-identical* results. Where the local frame
works around the coarse clock, this removes the dependence on it.

The third instance sits in the nonlinear solve, and this one is not confined to half precision. Its
convergence test is `rfₐ ≤ f_abstol + f_reltol·‖F(x₀)‖` (see [Solver tolerances](@ref)), and an
`f_abstol` pinned to `8eps(Float64)` regardless of the run's precision is unsatisfiable: a
half-precision residual bottoms out near `eps(T)`, so every implicit solve burns its full
1000-iteration budget on every step without ever being *reported* as non-convergent. Scaling the
tolerance to the run's own precision leaves `Float64` results bit-identical and the rest unchanged to
round-off, and makes the sweep 2.6–40× faster.

The `Float64` reference integrations have the same disease from a different cause. The 4D
Lotka–Volterra `Gauss(8)` reference cannot reach `8eps(Float64)` either — not because of precision but
because of *conditioning* — and so exhausts all 1000 iterations on every step, silently, for a solution
the study treats as ground truth. Scaling the floor with the stage size covers it: **two to three
orders of magnitude faster**, and the loosened and capped solutions agree to `1.4e-12`, so the capped
solve had converged all along without being able to certify it. Note the perverse coupling: the
relative term is measured at the *initial guess*, so the better the extrapolation, the smaller
`‖F(x₀)‖` and the more the test collapses onto the bare absolute floor — improving the guess makes
convergence harder to certify.

Both factors are supplied by the stack itself: `GeometricIntegratorsBase` defaults `f_abstol` to
`max(8, solversize(method, problem)) · eps(datatype(problem))`, precision- and size-scaled at once,
and `SimpleSolvers` gives a stalling solve the stagnation exit that the residual test cannot provide.
The study passes no tolerances of its own.

Five lessons generalise beyond this repository:

* **Distinguish the precision of the state from the precision of the bookkeeping.** They fail at very
  different horizons, and conflating them makes a solvable problem look like a hard limit.
* **Prefer the invariant quantity to the reconstructed one.** All three instances have the same
  shape: something already known exactly — the timestep, a tableau node, the machine epsilon of the
  working type — gets reconstructed from accumulated absolute values instead of used directly.
  In `Float64` the difference is invisible; in half precision it is the whole story.
* **Audit the tolerances too, not just the arrays.** A type-purity gate on the state will not catch a
  convergence threshold written for `Float64`. Asking a `BFloat16` solve for a `Float64` residual is
  not a stricter test, just an unsatisfiable one — and the `Float64` reference shows the same thing
  can happen from conditioning alone.
* **Treat "did not converge" as a result, not a warning.** Each of these degrades accuracy or wastes
  most of the compute while emitting nothing worse than a log line, so it can run unnoticed for a
  long time. A solve that exhausts its iteration budget on every step should be as visible as one
  that throws.
* **A "half precision breaks the solver" symptom may be an accumulated-time problem.** On the double
  pendulum a global clock costs the `Float16` implicit solves a `NaN` in the Newton direction, and
  swapping the initial-guess extrapolation only papers over it. The cause is that Hermite reads its
  interval off the stored history times, and `Float16` resolves only about `0.004` near `t = 5` —
  comparable to `Δt = 0.01` — so the extrapolation parameter goes materially wrong long before the
  clock saturates outright. In the local frame it is exact and the `NaN`s do not occur at any
  precision, so the default `HermiteExtrapolation` is the better guess everywhere.

## Genuine limits of half precision (not bugs)

What remains once the clock is out of the way is a much shorter list, and it is real:

* **The round-off floor on the achievable energy error.** This is the honest, unavoidable cost of
  reduced precision, and the one the study is really about: the bounded energy error of a symplectic
  method cannot go below what the working precision can represent.
* **One non-convergent solve on the Hamiltonian problems.** Across all eight Hamiltonian runs — every
  method at every precision — exactly one integration fails: `Implicit Euler` at `Float16` on the
  coarse-step pendulum (`Δt = 1`, `t ≤ 10 000`), with a `NaN` in the Newton direction. That is a
  first-order dissipative method at a very coarse step, i.e. the case with the least to recommend it
  anyway. Everything else completes, at all four precisions, over the full horizon.
* **The degenerate-Lagrangian variational integrators break down in half precision.** Across the four
  Lotka–Volterra runs, `VPRK Gauss(1)` fails with a `NaN` in the nonlinear-solver direction in *every*
  scenario at *both* 16-bit formats, and `PMVI Midpoint` in all but one; `CMDVI` and plain implicit
  midpoint each fail in one scenario. These are *singular* systems whose one-form involves `log(q)`, so
  the nonlinear iterate has to stay in the positive orthant — and 8–11 significand bits are not enough
  to keep it there. The iterate steps to a negative coordinate, the logarithm is evaluated there, and
  because GeometricProblems generates its vector fields with `nanmath = true` the model returns `NaN`
  rather than throwing a `DomainError` — so the breakdown surfaces one level up, in the solver.
  This is the one family where the failures are not a clock or tolerance artefact,
  and notably it is *not* a clean BFloat16-versus-Float16 story: the two 16-bit formats fail on
  overlapping but not identical sets, which is what one should expect when the breakdown turns on
  whether a particular iterate happens to step outside the domain.
* **Problem-dependent robustness.** The Toda lattice, whose bump initial data keeps the state
  bounded, runs every method at every precision — considerably more half-precision-friendly than the
  stiff, dimensional double pendulum.

All such failures are caught per run, reported as skips, and never abort the sweep.

## Implementation-detail and variational comparisons

Two further comparison groups isolate finer effects:

* **Partitioned Gauss(2) variants.** On the Hamiltonian problems, four algebraically
  equivalent forms of the 2-stage Gauss rule (symplectic-by-construction vs. by-duplication, with
  and without the rounding-compensation coefficients `â, b̂, ĉ`) coincide on the linear harmonic
  oscillator but separate on the nonlinear problems, where the tableau construction and the
  compensated-summation coefficients leave a visible imprint on the energy-error fine structure —
  most so at Float64, where the floor is not precision-limited.
* **Variational integrators on degenerate Lagrangians.** The Lotka–Volterra systems are compared
  with several flavours of the implicit midpoint rule that apply to degenerate (IODE/LODE) systems.
  They agree to their common order but differ in reduced-precision energy behaviour; `CMDVI`
  integrates the 2D system but not the 4D one (its iterate leaves the positive orthant and the solve
  breaks down there).

## Practical takeaways

* Use a **symplectic** (or at least symmetric) integrator whenever long-time energy behaviour
  matters; the qualitative advantage is independent of precision, and holds all the way down to
  `BFloat16`.
* **Float32** is a reasonable working precision for these problems: it preserves the geometric
  behaviour with a modest error floor and is far more robust than either 16-bit format.
* **Half precision is more usable than it first appears**, provided the clock, the initial-guess node
  and the solver tolerance are not themselves tied to it. With those separated out, both 16-bit
  formats run every method over every horizon studied on the four *Hamiltonian* problems; what they
  cost there is accuracy (an error floor of `1e-3`–`1e-1`) and, on chaotic systems, predictability
  horizon. The degenerate-Lagrangian problems are the exception, and `BFloat16` genuinely fails on
  them.
* **Between the two 16-bit formats, prefer `Float16` for this kind of work.** `BFloat16`'s wider
  exponent range is not useful for bounded Hamiltonian systems, and the three significand bits it
  costs show up directly in the error floor and in how early the time grid saturates.
* **Do not store accumulating quantities at the working precision by default.** Time is the example
  here, but the lesson is general: a monotonically growing counter in `BFloat16` runs out after
  ~256 increments, long before the state precision becomes the binding constraint.
