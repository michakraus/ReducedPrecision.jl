# Reduced-precision study: mathematical pendulum, long-time / coarse-step variant.
#
# Same pipeline as pendulum.jl but with a coarse timestep Δt = 1 and a long horizon
# t = 10_000 (nt = 10_000 steps). Solution error is measured against a Float64 Gauss(8)
# reference on the same grid; the reference integration is guarded so the energy-error plot
# is still produced even if the high-order solve does not converge at this coarse step.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.Pendulum: podeproblem, hamiltonian

const t₀ = 0.0
const Δt = 1.0
const nt = 10_000
const t₁ = nt * Δt

make_problem(::Type{T}) where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "pendulum_longtime_energy_error.png"),
    title = "Pendulum — Relative Energy Error (Δt = 1, t ≤ 10⁴)")

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = try
    integrate(make_problem(Float64), Gauss(8))
catch e
    @warn "reference integration failed; skipping solution-error plot" error = sprint(showerror, e)
    nothing
end

if reference !== nothing
    plot_solution_error(runs, reference;
        path  = joinpath(plotdir, "pendulum_longtime_solution_error.png"),
        title = "Pendulum — Solution Error (Δt = 1, t ≤ 10⁴, vs. Float64 Gauss(8))")

    plot_solution(runs; reference = reference,
        path   = joinpath(plotdir, "pendulum_longtime_solution.png"),
        title  = "Pendulum — Phase-Space Trajectory (Δt = 1, t ≤ 10⁴)",
        xlabel = "q", ylabel = "p")
end
