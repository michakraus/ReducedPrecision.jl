# Reduced-precision study: mathematical pendulum.
#
# Same pipeline as the harmonic oscillator, using the partitioned (PODE) form. The pendulum
# has no closed-form solution, so the solution error is measured against a high-precision
# Float64 reference computed with the high-order symplectic Gauss(8) method on the same grid.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.Pendulum: podeproblem, hamiltonian

# Horizon capped at t = 100 (nt = 1000) so Float16 can resolve the time grid at Δt = 0.1
# (see harmonic_oscillator.jl for the rationale).
const t₀ = 0.0
const Δt = 0.1
const nt = 1_000
const t₁ = nt * Δt

make_problem(::Type{T}) where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt))

# The pendulum Hamiltonian is parameter-free (l = m = g = 1, so H = p²/2 + cos q); wrap it in the
# 4-argument (t, q, p, params) form expected by `energy_error` / `compute_invariant_error`.
ham(t, q, p, params) = hamiltonian(t, q, p)

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_energy_error(runs, ham;
    path  = joinpath(plotdir, "pendulum_energy_error.png"),
    title = "Pendulum — Relative Energy Error (Δt = 0.1, t ≤ 100)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "pendulum_solution_error.png"),
    title = "Pendulum — Solution Error (Δt = 0.1, t ≤ 100, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "pendulum_solution.png"),
    title  = "Pendulum — Phase-Space Trajectory (Δt = 0.1, t ≤ 100)",
    xlabel = "q", ylabel = "p")
