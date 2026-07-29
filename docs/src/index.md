# ReducedPrecision.jl

*Accuracy and long-time stability of geometric vs. non-geometric integrators in reduced
floating-point precision.*

## Overview

`ReducedPrecision` studies how numerical integrators behave when run in **BFloat16**, **Float16**,
**Float32** and **Float64**, and in particular how *geometric* (symplectic) integrators compare to
*non-geometric* ones with respect to

* **accuracy** — the error of the computed solution relative to a reference, and
* **long-time stability** — whether the energy error stays bounded over long integration times.

The study uses six example problems from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl), each integrated with a
range of symplectic and non-symplectic methods of
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl), and each run in two
scenarios (a fine short-horizon run and a coarser long-horizon run). The four Hamiltonian problems
(harmonic oscillator, pendulum, double pendulum, Toda lattice) compare Euler-type methods, higher-
order Runge–Kutta methods, and a group of partitioned Gauss(2) variants; the two
degenerate-Lagrangian problems (Lotka–Volterra 2D and 4D) compare several flavours of variational
implicit-midpoint integrator.

The implicit solves use the trust-region **`DogLeg`** nonlinear solver of
[SimpleSolvers.jl](https://github.com/JuliaGNI/SimpleSolvers.jl) by default, which is more robust in
reduced precision than a line-search Newton iteration.

A central goal of the implementation is **type purity**: every library in the stack
(`GeometricIntegrators`, `GeometricIntegratorsBase`, `GeometricSolutions`, `GeometricEquations`,
`GeometricBase`, `SimpleSolvers`) must honour the requested precision and never silently promote to
`Float64`. This is checked for every run (see [Methodology](@ref)).

## Key findings at a glance

* **Symplectic integrators conserve energy.** Symplectic Euler A/B and the implicit midpoint rule
  keep the energy error *bounded* (oscillating around a small value) over the entire integration,
  whereas explicit Euler grows without bound and implicit Euler dissipates. This is the expected
  qualitative difference between geometric and non-geometric integrators, and it holds at every
  precision.
* **Precision sets the error floor.** For the energy-conserving methods the size of the (bounded)
  energy error is set by the working precision — for the harmonic oscillator at `Δt = 0.1`, the
  implicit midpoint rule settles at roughly `1e-2` (BFloat16), `9e-3` (Float16), `2e-6` (Float32) and
  `5e-15` (Float64).
* **Exponent range buys nothing; significand bits do.** `BFloat16` and `Float16` are both 16-bit, but
  `BFloat16` spends three of its significand bits on `Float32`'s exponent range. On these bounded
  Hamiltonian problems that range is never needed, so `BFloat16` is simply the coarser of the two —
  its error floor sits a factor of ≈ 8 above `Float16`'s, exactly the ratio of their `eps`.
* **Reduced precision is type-pure.** No implicit promotion to `Float64` occurs in any library, for
  any method or problem, including the hand-built half-precision constructions of the double
  pendulum and Toda lattice.
* **Most of what looked like a half-precision limit was a bookkeeping limit.** Three quantities were
  being carried or compared at the working precision without needing to be — the *clock*, the initial
  guess's *interpolation node*, and the solver's *residual tolerance* — and each produced failures
  that read as hardware limits. A `T`-typed clock stops advancing once `ulp(t) ≥ Δt` (from `t ≈ 128`
  in `Float16`, `t ≈ 16` in `BFloat16` at `Δt = 0.1`), which used to cap the horizons studied and
  break the implicit solves. Stepping in a [local time frame](@ref "Time stepping in a local frame")
  and taking the [initial guess](@ref "Initial guess") from the tableau's normalised nodes remove that
  dependence exactly for these autonomous problems. Across all eight Hamiltonian runs — every method at
  every precision — exactly one integration now fails, against a whole column of them before.
* **What remains is genuinely precision-limited:** the round-off floor on the achievable energy error,
  and the degenerate-Lagrangian Lotka–Volterra variational integrators, which break down in half
  precision — those are singular systems with a `log(q)` one-form, so the nonlinear iterate must stay
  in the positive orthant, and 8–11 significand bits are not enough to keep it there.

See [Findings](@ref) for the full discussion.

## Running the experiments

```julia
# from the package root, with the project activated
julia --project=. scripts/harmonic_oscillator.jl
julia --project=. scripts/pendulum.jl
julia --project=. scripts/double_pendulum.jl
julia --project=. scripts/toda_lattice.jl
julia --project=. scripts/lotka_volterra_2d.jl
julia --project=. scripts/lotka_volterra_4d.jl
# coarser-step variants
julia --project=. scripts/harmonic_oscillator_longtime.jl
julia --project=. scripts/pendulum_longtime.jl
julia --project=. scripts/double_pendulum_longtime.jl
julia --project=. scripts/toda_lattice_longtime.jl
julia --project=. scripts/lotka_volterra_2d_longtime.jl
julia --project=. scripts/lotka_volterra_4d_longtime.jl
```

Alternatively, regenerate every figure at once with the runner (it discovers `scripts/*.jl`
automatically and runs them in a single Julia session, so the shared packages compile only once):

```bash
julia --project=. scripts/run_all.jl
```

Each script writes its figures to `plots/`. The documentation embeds those figures, so the scripts
must be run before building the docs.
