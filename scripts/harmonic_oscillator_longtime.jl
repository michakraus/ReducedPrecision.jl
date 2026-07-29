# Reduced-precision study: harmonic oscillator, long-time / coarse-step variant.
#
# Same pipeline as harmonic_oscillator.jl but with a coarse timestep Δt = 1 and a long
# horizon t = 10_000 (nt = 10_000 steps), stressing long-time stability. A T-typed global clock
# would stop advancing well before the end at half precision (t + 1 == t past t ≈ 2048 in Float16,
# t ≈ 256 in BFloat16); the local-frame stepping of `integrate_bounded` carries the full horizon at
# every precision.

using ReducedPrecision
using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution
import GeometricProblems.HarmonicOscillator as HO

const t₀ = 0.0
const Δt = 1.0
const nt = 10_000
const t₁ = nt * Δt

make_problem(::Type{T}) where {T} =
    podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(t₀), T(t₁)), timestep = T(Δt))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

# reference (analytic)
reference = exact_solution(make_problem(Float64))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "harmonic_oscillator_energy_error_dt_$(Δt).png"),
    title = "Harmonic Oscillator — Relative Energy Error (Δt = 1, t ≤ 10⁴)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "harmonic_oscillator_solution_error_dt_$(Δt).png"),
    title = "Harmonic Oscillator — Solution Error (Δt = 1, t ≤ 10⁴, vs. analytic)")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "harmonic_oscillator_solution_dt_$(Δt).png"),
    title  = "Harmonic Oscillator — Phase-Space Trajectory (Δt = 1, t ≤ 10⁴)",
    xlabel = "q", ylabel = "p")
