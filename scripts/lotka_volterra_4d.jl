# Reduced-precision study: Lotka–Volterra 4D (degenerate Lagrangian).
#
# Built from the LotkaVolterra4dLagrangian module as a LODE, using the quasi-canonical reduced
# gauge (A = A_quasicanonical_reduced, plus the exact one-form B) so the discrete system is
# non-singular. As for the 2D case, only variational / implicit-midpoint methods apply; here the
# comparison (`LV4D_METHODS`) is Implicit Midpoint, VPRK(Gauss(1)) and PMVI midpoint. CMDVI is
# omitted — it does not converge on this 4D degenerate Lagrangian (its iterate leaves the positive
# orthant and the solve breaks down). No closed-form solution — the error is measured against a
# Float64 Gauss(8) reference; the Hamiltonian is the (NaNMath-safe) LotkaVolterra4d energy.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.LotkaVolterra4d: hamiltonian
import GeometricProblems.LotkaVolterra4dLagrangian as LV

const Δt = 0.01
const nt = 1_000
const t₁ = nt * Δt

# Degenerate-Lagrangian (LODE) form at precision T. The gauge matrices A/B are exact rationals and
# stay unconverted (they promote to T when evaluated); only q₀, timespan, timestep, parameters are T.
make_problem(::Type{T}) where {T} = LV.lodeproblem(T.(LV.q₀), LV.A_quasicanonical_reduced, LV.B;
    timespan   = (T(0), T(t₁)),
    timestep   = T(Δt),
    parameters = LV.default_parameters(T))

# Hamiltonian closure (energy depends on q only; energy_error evaluates it as (t,q,p,params)).
ham(t, q, p, params) = hamiltonian(t, q, params)

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem; methods = LV4D_METHODS)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_energy_error(runs, ham; groups = LV4D_GROUPS,
    path  = joinpath(plotdir, "lotka_volterra_4d_energy_error_dt_$(Δt).png"),
    title = "Lotka–Volterra 4D — Relative Energy Error (Δt = 0.01, t ≤ 10)")

plot_solution_error(runs, reference; groups = LV4D_GROUPS,
    path  = joinpath(plotdir, "lotka_volterra_4d_solution_error_dt_$(Δt).png"),
    title = "Lotka–Volterra 4D — Solution Error (Δt = 0.01, t ≤ 10, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference, groups = LV4D_GROUPS,
    path   = joinpath(plotdir, "lotka_volterra_4d_solution_dt_$(Δt).png"),
    title  = "Lotka–Volterra 4D — Configuration-Space Trajectory (Δt = 0.01, t ≤ 10)",
    xlabel = "q₁", ylabel = "q₂")
