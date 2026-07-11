# Reduced-precision study: harmonic oscillator, long-time / coarse-step variant.
#
# Same pipeline as harmonic_oscillator.jl but with a coarse timestep Δt = 1 and a long
# horizon t = 10_000 (nt = 10_000 steps), stressing long-time stability. Note that at this
# horizon low precisions cannot resolve the time grid (in Float16, t + 1 == t once t exceeds
# ~2048), so implicit methods that rely on distinct successive times may fail there — those
# runs are caught and reported as skips.

using ReducedPrecision
using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution

const t₀ = 0.0
const Δt = 1.0
const nt = 10_000
const t₁ = nt * Δt

make_problem(::Type{T}) where {T} = podeproblem(T; timespan = (T(t₀), T(t₁)), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

# reference (analytic)
reference = exact_solution(make_problem(Float64))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "harmonic_oscillator_longtime_energy_error.png"),
    title = "Harmonic Oscillator — Relative Energy Error (Δt = 1, t ≤ 10⁴)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "harmonic_oscillator_longtime_solution_error.png"),
    title = "Harmonic Oscillator — Solution Error (Δt = 1, t ≤ 10⁴, vs. analytic)")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "harmonic_oscillator_longtime_solution.png"),
    title  = "Harmonic Oscillator — Phase-Space Trajectory (Δt = 1, t ≤ 10⁴)",
    xlabel = "q", ylabel = "p")
