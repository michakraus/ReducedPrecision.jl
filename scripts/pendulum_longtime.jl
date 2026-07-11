# Reduced-precision study: mathematical pendulum, long-time / coarse-step variant.
#
# Same pipeline as pendulum.jl but with a coarse timestep Δt = 1 and a long horizon
# t = 10_000 (nt = 10_000 steps). The solution error is measured against a Float64 Gauss(8)
# reference computed at the *fine* step (Δt_ref = 0.1, as in the short scenario) and subsampled
# onto the coarse grid, so the reference is trustworthy independent of the coarse step; the
# reference integration is guarded so the energy-error plot is still produced if it fails.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.Pendulum: podeproblem, hamiltonian

const t₀ = 0.0
const Δt = 1.0
const nt = 10_000
const t₁ = nt * Δt
const Δt_ref = 0.1   # fine reference step (matches the short scenario)

make_problem(::Type{T})   where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt))
make_reference(::Type{T}) where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt_ref))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "pendulum_energy_error_dt_$(Δt).png"),
    title = "Pendulum — Relative Energy Error (Δt = 1, t ≤ 10⁴)")

# high-precision reference (Float64, high-order symplectic, fine step, subsampled to the grid)
reference = try
    integrate(make_reference(Float64), Gauss(8))
catch e
    @warn "reference integration failed; skipping solution-error plot" error = sprint(showerror, e)
    nothing
end

if reference !== nothing
    plot_solution_error(runs, reference;
        path  = joinpath(plotdir, "pendulum_solution_error_dt_$(Δt).png"),
        title = "Pendulum — Solution Error (Δt = 1, t ≤ 10⁴, vs. Float64 Gauss(8) at Δt = 0.1)")

    plot_solution(runs; reference = reference,
        path   = joinpath(plotdir, "pendulum_solution_dt_$(Δt).png"),
        title  = "Pendulum — Phase-Space Trajectory (Δt = 1, t ≤ 10⁴)",
        xlabel = "q", ylabel = "p")
end
