# Reduced-precision study: mathematical pendulum.
#
# Same pipeline as the harmonic oscillator, using the partitioned (PODE) form. The pendulum
# has no closed-form solution, so the solution error is measured against a high-precision
# Float64 reference computed with the high-order symplectic Gauss(8) method on the same grid.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.Pendulum: podeproblem, hamiltonian
import GeometricProblems.Pendulum as PD

# Nominal horizon t = 1000 (nt = 10_000) at Δt = 0.1. The Float16 final time is capped to the last
# resolvable grid point via `capped_final_time` (the Float16 grid saturates at t ≈ 128 for Δt = 0.1,
# which breaks the implicit methods' Hermite initial guess); see harmonic_oscillator.jl. Float32/
# Float64 keep the full t ≤ 1000 horizon.
const t₀ = 0.0
const Δt = 0.1
const nt = 10_000
const t₁ = nt * Δt

# GeometricProblems v0.7.0 dropped the `podeproblem(::Type{T})` precision constructor, so the
# T-typed initial conditions are built here from the module defaults.
make_problem(::Type{T}) where {T} =
    podeproblem(T.(PD.q₀), T.(PD.p₀);
        timespan = (T(t₀), T(capped_final_time(T, t₁, Δt))), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "pendulum_energy_error_dt_$(Δt).png"),
    title = "Pendulum — Relative Energy Error (Δt = 0.1, t ≤ 10^4)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "pendulum_solution_error_dt_$(Δt).png"),
    title = "Pendulum — Solution Error (Δt = 0.1, t ≤ 10^4, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "pendulum_solution_dt_$(Δt).png"),
    title  = "Pendulum — Phase-Space Trajectory (Δt = 0.1, t ≤ 10^4)",
    xlabel = "q", ylabel = "p")
