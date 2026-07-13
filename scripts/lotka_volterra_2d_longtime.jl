# Reduced-precision study: Lotka–Volterra 2D (singular / degenerate Lagrangian), coarse-step variant.
#
# Same pipeline as lotka_volterra_2d.jl but with a coarser timestep Δt = 0.1 and a longer horizon
# t = 100 (nt = 1000 steps), stressing the variational integrators' long-time behaviour. The
# solution error is measured against a Float64 Gauss(8) reference computed at the *fine* step
# (Δt_ref = 0.01, as in the short scenario) and subsampled onto the coarse grid.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.LotkaVolterra2dSingular: lodeproblem, hamiltonian
import GeometricProblems.LotkaVolterra2dSingular as LV

const Δt = 0.1
const nt = 1_000
const t₁ = nt * Δt
const Δt_ref = 0.01  # fine reference step (matches the short scenario)

_lv_problem(::Type{T}, dt) where {T} = lodeproblem(T.(LV.q₀);
    timespan = (T(0), T(t₁)), timestep = T(dt), parameters = LV.default_parameters(T))

make_problem(::Type{T})   where {T} = _lv_problem(T, Δt)
make_reference(::Type{T}) where {T} = _lv_problem(T, Δt_ref)

ham(t, q, p, params) = hamiltonian(t, q, params)

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem; methods = LV2D_METHODS)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, fine step, subsampled to the grid)
reference = try
    integrate(make_reference(Float64), Gauss(8))
catch e
    @warn "reference integration failed; skipping solution-error and trajectory plots" error = sprint(showerror, e)
    nothing
end

plot_energy_error(runs, ham; groups = LV2D_GROUPS,
    path  = joinpath(plotdir, "lotka_volterra_2d_energy_error_dt_$(Δt).png"),
    title = "Lotka–Volterra 2D — Relative Energy Error (Δt = 0.1, t ≤ 100)")

if reference !== nothing
    plot_solution_error(runs, reference; groups = LV2D_GROUPS,
        path  = joinpath(plotdir, "lotka_volterra_2d_solution_error_dt_$(Δt).png"),
        title = "Lotka–Volterra 2D — Solution Error (Δt = 0.1, t ≤ 100, vs. Float64 Gauss(8) at Δt = 0.01)")

    plot_solution(runs; reference = reference, groups = LV2D_GROUPS,
        path   = joinpath(plotdir, "lotka_volterra_2d_solution_dt_$(Δt).png"),
        title  = "Lotka–Volterra 2D — Configuration-Space Trajectory (Δt = 0.1, t ≤ 100)",
        xlabel = "q₁", ylabel = "q₂")
end
