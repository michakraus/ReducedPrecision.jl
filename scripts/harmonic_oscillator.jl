# Reduced-precision study: harmonic oscillator.
#
# This is the reference pipeline. It additionally runs the precision-verification gate
# (`assert_precision`) for every (method, precision) combination, proving that none of the
# involved libraries silently promotes to Float64. The harmonic oscillator has a closed-form
# exact solution, so the solution error is measured against the analytic reference.

using ReducedPrecision
using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution
import GeometricProblems.HarmonicOscillator as HO

# Horizon t = 1000 (nt = 10_000) at Δt = 0.1, at every precision. A T-typed global clock would stop
# advancing long before that in half precision (around t ≈ 128 in Float16 and t ≈ 16 in BFloat16;
# see `capped_final_time`), but `integrate_bounded` steps in a local time frame, so the horizon is
# not limited by the clock.
const t₀ = 0.0
const Δt = 0.1
const nt = 10_000
const t₁ = nt * Δt

# Partitioned form at precision T, with initial conditions *and* time all in T. GeometricProblems has
# no `podeproblem(::Type{T})` precision constructor, so the T-typed initial conditions are built here
# from the module defaults.
make_problem(::Type{T}) where {T} =
    podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(t₀), T(t₁)), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

# --- verification gate: no implicit type conversion -----------------------------------
verify_precision(runs)

# --- reference & plots ----------------------------------------------------------------
reference = exact_solution(make_problem(Float64))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "harmonic_oscillator_energy_error_dt_$(Δt).png"),
    title = "Harmonic Oscillator — Relative Energy Error (Δt = 0.1, t ≤ 1000)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "harmonic_oscillator_solution_error_dt_$(Δt).png"),
    title = "Harmonic Oscillator — Solution Error (Δt = 0.1, t ≤ 1000, vs. analytic)")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "harmonic_oscillator_solution_dt_$(Δt).png"),
    title  = "Harmonic Oscillator — Phase-Space Trajectory (Δt = 0.1, t ≤ 1000)",
    xlabel = "q", ylabel = "p")
