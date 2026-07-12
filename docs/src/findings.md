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
* **Implicit midpoint** and **Crank–Nicolson** (symmetric methods) nearly conserve energy, with a
  bounded error set by the precision.
* **RK4** is accurate over short times but, being non-symplectic, exhibits a slow energy drift over
  long horizons.

## The role of precision

* For the energy-conserving methods the **precision sets the floor** of the (bounded) energy error.
  For the harmonic oscillator this is roughly `1e-2` (Float16), `1e-6` (Float32) and `1e-15`
  (Float64); the other problems show the same ordering.
* For the low-order methods the **solution error** is often dominated by the *truncation* error of
  the method rather than round-off, so Float32 and Float64 solution errors can be nearly identical
  while only Float16 shows a round-off floor.
* Reducing precision therefore mainly matters for (a) the achievable energy-conservation floor and
  (b) the robustness of the implicit solves and the time grid — not so much for the truncation-
  limited accuracy of the low-order methods.

## Type purity

`verify_precision` passes for **every successful run across all six problems**: `datatype`,
`timetype`, and the element types of the stored `q`/`p` arrays all equal the requested precision.
No library in the stack (`GeometricIntegrators`, `GeometricIntegratorsBase`, `GeometricSolutions`,
`GeometricEquations`, `GeometricBase`, `SimpleSolvers`) silently promotes to `Float64`, including
for the hand-built `Float16`/`Float32` constructions of the double pendulum and Toda lattice. The
key requirement is that **both** the initial conditions **and** the timespan/timestep are created at
the target precision, so that the data type and the time type agree.

## Genuine limits of half precision (not bugs)

These are real properties of `Float16`, surfaced by the study rather than worked around:

* **Time-grid saturation at long horizons.** Once `t` exceeds the range where successive time
  stamps differ in `Float16` (e.g. `ulp(10000) ≈ 8 ≫ Δt = 1`), the implicit methods' Hermite
  initial guess sees `t₀ == t₁` and fails. This is why the harmonic-oscillator and pendulum coarse
  scenarios (which run out to `t ≤ 10 000`) drop some `Float16` implicit runs, whereas the shorter
  scenarios — including the double pendulum and Toda lattice, whose coarse runs share their short
  run's horizon — populate the full method × precision matrix.
* **Non-convergent implicit solves on stiff systems.** On the double pendulum the `Float16`
  Crank–Nicolson solve still produces a `NaN` direction and drops out — even after switching the
  nonlinear solve from the trust-region `DogLeg` to the line-search `Newton`/`Backtracking` solver
  (capped at 100 iterations). It is a genuinely half-precision-limited solve, not a solver-tuning
  issue.
* **Problem-dependent robustness.** The Toda lattice, whose bump initial data keeps the state
  bounded, runs every method at every precision in the short scenario — considerably more
  `Float16`-friendly than the stiff, dimensional double pendulum.

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
  matters; the qualitative advantage is independent of precision.
* **Float32** is a reasonable working precision for these problems: it preserves the geometric
  behaviour with a modest error floor and is far more robust than `Float16`.
* **Float16** is viable only for short horizons and well-conditioned (non-stiff, bounded) problems;
  its limited exponent/mantissa range breaks both the time grid and the implicit solves otherwise.
