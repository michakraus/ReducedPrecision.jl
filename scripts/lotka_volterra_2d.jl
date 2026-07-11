# Reduced-precision study: Lotka–Volterra 2D (singular / degenerate Lagrangian).
#
# The Lotka–Volterra system is a degenerate Lagrangian system, posed here as a LODE via the
# LotkaVolterra2dSingular module (the "singular" gauge is the one the CMDVI variational integrator
# needs). On such systems the explicit / symplectic-Euler / DIRK methods of the main study do not
# apply; the comparison (`LV2D_METHODS`) is instead between several flavours of the implicit midpoint
# rule that do: Implicit Midpoint, VPRK(Gauss(1)), PMVI midpoint, and CMDVI. No closed-form
# solution exists — the solution error is measured against a Float64 Gauss(8) reference.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.LotkaVolterra2dSingular: lodeproblem, hamiltonian
import GeometricProblems.LotkaVolterra2dSingular as LV

const Δt = 0.01
const nt = 1_000
const t₁ = nt * Δt

# Degenerate-Lagrangian (LODE) form at precision T, with initial conditions and time all in T.
make_problem(::Type{T}) where {T} = lodeproblem(T.(LV.q₀);
    timespan   = (T(LV.timespan[begin]), T(t₁)),
    timestep   = T(Δt),
    parameters = map(T, LV.default_parameters))

# Hamiltonian closure: the energy depends on q only, but energy_error evaluates it as (t,q,p,params).
ham(t, q, p, params) = hamiltonian(t, q, params)

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem; methods = LV2D_METHODS)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_energy_error(runs, ham; groups = LV2D_GROUPS,
    path  = joinpath(plotdir, "lotka_volterra_2d_energy_error_dt=$(Δt).png"),
    title = "Lotka–Volterra 2D — Relative Energy Error (Δt = 0.01, t ≤ 10)")

plot_solution_error(runs, reference; groups = LV2D_GROUPS,
    path  = joinpath(plotdir, "lotka_volterra_2d_solution_error_dt=$(Δt).png"),
    title = "Lotka–Volterra 2D — Solution Error (Δt = 0.01, t ≤ 10, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference, groups = LV2D_GROUPS,
    path   = joinpath(plotdir, "lotka_volterra_2d_solution_dt=$(Δt).png"),
    title  = "Lotka–Volterra 2D — Configuration-Space Trajectory (Δt = 0.01, t ≤ 10)",
    xlabel = "q₁", ylabel = "q₂")
