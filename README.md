# ReducedPrecision

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://michakraus.github.io/ReducedPrecision.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://michakraus.github.io/ReducedPrecision.jl/dev/)
[![Build Status](https://github.com/michakraus/ReducedPrecision.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/michakraus/ReducedPrecision.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/michakraus/ReducedPrecision.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/michakraus/ReducedPrecision.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/ReducedPrecision.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/ReducedPrecision.html)

`ReducedPrecision` studies how numerical integrators behave in reduced floating-point precision
(`BFloat16`, `Float16`, `Float32`, `Float64`), and in particular how **geometric (symplectic)**
integrators compare to **non-geometric** ones with respect to accuracy and long-time stability.

Using six example problems from
[GeometricProblems.jl](https://github.com/JuliaGNI/GeometricProblems.jl) — the harmonic oscillator,
the pendulum, the double pendulum, the Toda lattice, and the 2D and 4D Lotka–Volterra systems —
each is integrated with
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl) at all four
precisions and in two scenarios (a fine short-horizon run and a coarser one), and the energy error
and solution error are compared across methods.

The two 16-bit formats spend their bits differently: `BFloat16`
([BFloat16s.jl](https://github.com/JuliaMath/BFloat16s.jl)) keeps `Float32`'s exponent range at the
cost of three significand bits, so it is the coarsest of the four (`eps` `2⁻⁷` vs `Float16`'s
`2⁻¹⁰`) while being the hardest to overflow.

A `T`-typed clock stops advancing once `ulp(t) ≥ Δt` — well before the horizons studied here, and much
earlier for `BFloat16` than for `Float16` — which used to break the implicit methods outright. Two
changes remove that limit, so every precision runs the full horizon: the solution step is advanced in a
**local time frame** rather than along the problem's own time grid, and the implicit methods' initial
guess is extrapolated from the **tableau's normalised nodes** instead of from differenced absolute
times, which makes it independent of the clock altogether. See the
[methodology](https://michakraus.github.io/ReducedPrecision.jl/stable/methodology/).

The methods compared are:

- **Euler methods:** symplectic Euler A/B (geometric) against explicit and implicit Euler;
- **midpoint / fourth-order methods**, as a 2 × 2 of explicit vs. implicit at order 2 and order 4:
  explicit midpoint, explicit RK4, implicit midpoint (`Gauss(1)`) and implicit RK4 (`Gauss(2)`);
- **partitioned Gauss(2) variants:** four algebraically-equivalent forms of the 2-stage
  Gauss rule (symplectic-by-construction vs. by-duplication, with/without the rounding-compensation
  coefficients `â, b̂, ĉ`), isolating implementation-detail effects on energy conservation.

The two **Lotka–Volterra** problems are degenerate Lagrangian systems (posed as LODEs), on which the
above methods do not apply; they are compared instead with several flavours of variational
implicit-midpoint integrator (`Implicit Midpoint`, `VPRK(Gauss(1))`, `PMVImidpoint`, and — 2D only —
`CMDVI`).

The implicit solves use the trust-region **`DogLeg`** nonlinear solver (`SimpleSolvers`) by default,
which is more robust in reduced precision than a line-search Newton iteration, with the solve's
absolute residual tolerance scaled to the working precision and to the size of the stage system. The
stack's default is a fixed `8eps(Float64)`, which neither a half-precision residual nor the degenerate
4D Lotka–Volterra `Gauss(8)` reference can reach — both then exhausted all 1000 iterations on every
step, silently.

A central design goal is **type purity**: every library in the stack (`GeometricIntegrators`,
`GeometricIntegratorsBase`, `GeometricSolutions`, `GeometricEquations`, `GeometricBase`,
`SimpleSolvers`) must honour the requested precision and never silently promote to `Float64`. This
is asserted for every run.

## Usage

Each experiment is a script in [`scripts/`](scripts) that runs the full sweep, verifies precision
purity, and writes its figures to `plots/`:

```julia
julia --project=. scripts/harmonic_oscillator.jl          # short: Δt = 0.1, t ≤ 1000
julia --project=. scripts/harmonic_oscillator_longtime.jl # coarse: Δt = 1,  t ≤ 10000
# ... likewise for pendulum, double_pendulum, toda_lattice,
#     lotka_volterra_2d, and lotka_volterra_4d
```

To regenerate every figure at once, run all experiment scripts via the runner (it discovers
`scripts/*.jl` automatically, so new examples are picked up without editing it, and runs them in a
single Julia session so the shared packages compile only once):

```bash
julia --project=. scripts/run_all.jl
```

Output figure filenames encode the timestep (e.g. `…_dt_0.1_…`), so a problem's two scenarios are
distinguished by `Δt`.

The reusable pipeline (`run_study`, `verify_precision`, `energy_error`, `solution_error`, and the
CairoMakie plotting routines) lives in [`src/ReducedPrecision.jl`](src/ReducedPrecision.jl).

## Documentation

A summary of all experiments and findings is provided as a
[Documenter](https://github.com/JuliaDocs/Documenter.jl) site under [`docs/`](docs). After running
the scripts (so the figures exist), build it with:

```julia
julia --project=docs docs/make.jl
```
