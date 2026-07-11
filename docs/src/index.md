# ReducedPrecision.jl

*Accuracy and long-time stability of geometric vs. non-geometric integrators in reduced
floating-point precision.*

## Overview

`ReducedPrecision` studies how numerical integrators behave when run in **Float16**, **Float32**,
and **Float64**, and in particular how *geometric* (symplectic) integrators compare to
*non-geometric* ones with respect to

* **accuracy** — the error of the computed solution relative to a reference, and
* **long-time stability** — whether the energy error stays bounded over long integration times.

The study uses six example problems from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl), each integrated with a
range of symplectic and non-symplectic methods of
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl), and each run in two
scenarios (a fine short-horizon run and a coarser long-horizon run). The four Hamiltonian problems
(harmonic oscillator, pendulum, double pendulum, Toda lattice) compare Euler-type methods, higher-
order Runge–Kutta methods, and a group of partitioned Gauss(2) midpoint variants; the two
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
  energy error is set by the working precision — e.g. for the harmonic oscillator roughly
  `1e-2` (Float16), `1e-6` (Float32), `1e-15` (Float64).
* **Reduced precision is type-pure.** No implicit promotion to `Float64` occurs in any library, for
  any method or problem, including the hand-built `Float16`/`Float32` constructions of the double
  pendulum and Toda lattice.
* **Float16 has hard limits at long horizons.** Once time exceeds the range where `t + Δt` is
  distinguishable in `Float16`, the implicit methods' initial guess breaks down; and on stiff
  systems the `Float16` implicit solves can produce `NaN` directions even with the robust `DogLeg`
  solver. These are genuine properties of half precision, surfaced (not hidden) by the study.

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

Each script writes its figures to `plots/`. The documentation embeds those figures, so the scripts
must be run before building the docs.
