"""
    ReducedPrecision

Reusable pipeline for studying geometric (symplectic) vs non-geometric integrators
run in varying floating-point precision (Float16, Float32, Float64).

The study is driven by a per-problem `make_problem(T)` closure that must return a
partitioned problem (PODE/HODE) at precision `T` with **both** initial conditions and
timespan/timestep in `T` (so `datatype(prob) === timetype(prob) === T`). All methods in
the registry are then run across every precision, and energy-error / solution-error plots
are produced with CairoMakie.

Symplectic Euler (A/B) require a partitioned problem; the non-geometric explicit/implicit
Euler methods are represented by their Runge–Kutta tableau twins (`ExplicitEulerRK` /
`ImplicitEulerRK`), which auto-promote to partitioned RK on a PODE/HODE — so a single
problem form runs the whole method set for a fair comparison.

The implementation is split into logical units:

- `methods.jl`     — precisions, the `MethodSpec` registry and the plotting method groups;
- `study.jl`       — the `Run` type and `run_study` sweep;
- `diagnostics.jl` — precision-purity checks and the error metrics;
- `plotting.jl`    — the CairoMakie plotting routines.
"""
module ReducedPrecision

using GeometricIntegrators
using GeometricIntegratorsBase: Solution, solutionstep, current,
    HermiteExtrapolation, MidpointExtrapolation
import GeometricIntegratorsBase: initmethod, isimplicit, default_options
using GeometricSolutions
using GeometricBase
using GeometricEquations: parameters, GeometricProblem
using SimpleSolvers: DogLeg, Newton, StrongWolfe, Backtracking
using RungeKutta: Tableau, PartitionedTableau, SymplecticPartitionedTableau, TableauGauss
using CairoMakie

export PRECISIONS, MethodSpec, GEOMETRIC_METHODS, NONGEOMETRIC_METHODS, ALL_METHODS
export EULER_METHODS, OTHER_METHODS, GAUSS2_METHODS, METHOD_GROUPS
export LV2D_METHODS, LV4D_METHODS, LV2D_GROUPS, LV4D_GROUPS
export Run, run_study, integrate_bounded, assert_precision, verify_precision, capped_final_time
export DogLeg, Newton, StrongWolfe, Backtracking
export HermiteExtrapolation, MidpointExtrapolation
export energy_error, solution_error, timevalues
export plot_energy_error, plot_solution_error, plot_solution

include("methods.jl")
include("study.jl")
include("diagnostics.jl")
include("plotting.jl")

end # module
