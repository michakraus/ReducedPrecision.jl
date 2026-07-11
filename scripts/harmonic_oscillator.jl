# Reduced-precision study: harmonic oscillator.
#
# This is the reference pipeline. It additionally runs the precision-verification gate
# (`assert_precision`) for every (method, precision) combination, proving that none of the
# involved libraries silently promotes to Float64. The harmonic oscillator has a closed-form
# exact solution, so the solution error is measured against the analytic reference.

using ReducedPrecision
using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution

# Horizon capped at t = 100 (nt = 1000): still many oscillation periods, but within the
# time range Float16 can resolve at Δt = 0.1 (ulp(100) ≈ 0.06 < Δt). A longer horizon makes
# successive Float16 time stamps indistinguishable, which breaks the implicit methods'
# Hermite initial guess (t₀ == t₁).
const t₀ = 0.0
const Δt = 0.1
const nt = 1_000
const t₁ = nt * Δt

# Partitioned form at precision T, with initial conditions *and* time all in T.
make_problem(::Type{T}) where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

# --- verification gate: no implicit type conversion -----------------------------------
verify_precision(runs)

# --- reference & plots ----------------------------------------------------------------
reference = exact_solution(make_problem(Float64))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "harmonic_oscillator_energy_error_dt=$(Δt).png"),
    title = "Harmonic Oscillator — Relative Energy Error (Δt = 0.1, t ≤ 100)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "harmonic_oscillator_solution_error_dt=$(Δt).png"),
    title = "Harmonic Oscillator — Solution Error (Δt = 0.1, t ≤ 100, vs. analytic)")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "harmonic_oscillator_solution_dt=$(Δt).png"),
    title  = "Harmonic Oscillator — Phase-Space Trajectory (Δt = 0.1, t ≤ 100)",
    xlabel = "q", ylabel = "p")
