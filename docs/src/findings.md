# Findings

## Geometric vs. non-geometric integrators

Across all four problems and both scenarios the qualitative distinction is consistent and matches
the theory of geometric numerical integration:

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

`verify_precision` passes for **every successful run across all four problems**: `datatype`,
`timetype`, and the element types of the stored `q`/`p` arrays all equal the requested precision.
No library in the stack (`GeometricIntegrators`, `GeometricIntegratorsBase`, `GeometricSolutions`,
`GeometricEquations`, `GeometricBase`, `SimpleSolvers`) silently promotes to `Float64`, including
for the hand-built `Float16`/`Float32` constructions of the double pendulum and Toda lattice. The
key requirement is that **both** the initial conditions **and** the timespan/timestep are created at
the target precision, so that the data type and the time type agree.

## Genuine limits of half precision (not bugs)

These are real properties of `Float16`, surfaced by the study rather than worked around:

* **Time-grid saturation at long horizons.** Once `t` exceeds the range where successive time
  stamps differ in `Float16` (e.g. `ulp(1000) ≈ 0.5 ≫ Δt = 0.1`), the implicit methods' Hermite
  initial guess sees `t₀ == t₁` and fails. This is why the short scenarios cap the horizon so that
  the full method × precision matrix is populated, while in the long scenarios some `Float16`
  implicit runs drop out.
* **Non-convergent implicit solves on stiff systems.** On the double pendulum the `Float16` Newton
  iteration produces `NaN` directions; at `Δt = 1` the double pendulum's implicit solves (and even
  the `Gauss(8)` reference) fail at *every* precision because the step is comparable to the natural
  period.
* **Problem-dependent robustness.** The Toda lattice, whose bump initial data keeps the state
  bounded, runs every method at every precision in the short scenario — considerably more
  `Float16`-friendly than the stiff, dimensional double pendulum.

All such failures are caught per run, reported as skips, and never abort the sweep.

## Practical takeaways

* Use a **symplectic** (or at least symmetric) integrator whenever long-time energy behaviour
  matters; the qualitative advantage is independent of precision.
* **Float32** is a reasonable working precision for these problems: it preserves the geometric
  behaviour with a modest error floor and is far more robust than `Float16`.
* **Float16** is viable only for short horizons and well-conditioned (non-stiff, bounded) problems;
  its limited exponent/mantissa range breaks both the time grid and the implicit solves otherwise.
