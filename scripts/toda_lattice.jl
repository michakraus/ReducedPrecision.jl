# Reduced-precision study: Toda lattice (periodic N-site lattice, a completely integrable system).
#
# The TodaLattice module has NO `::Type{T}` constructor, so `make_problem(T)` hand-builds T-typed
# initial conditions, timespan, timestep and parameters (and passes the lattice size N). Its
# `hamiltonian` has a four-argument `(t, q, p, params)` method that recovers the lattice size from
# the state, so it can be handed to `plot_energy_error` directly. The one wrinkle compared to the
# other problems is that the state has N degrees of freedom, so the trajectory plot uses the phase
# space of the first lattice site. No closed-form solution; the solution error is measured against a
# Float64 Gauss(8) reference on the same grid.
#
# A modest lattice (N = 16) is used to keep the precision sweep tractable — the default N = 200
# soliton example would make the high-order reference and the implicit solves prohibitively slow.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.TodaLattice: hodeproblem, hamiltonian
import GeometricProblems.TodaLattice as TL

const N  = 16       # lattice size
const μ  = 0.3      # bump-width parameter of the initial condition
const t₀ = 0.0
const Δt = 0.1
const nt = 1_000
const t₁ = nt * Δt

function make_problem(::Type{T}) where {T}
    q₀     = T.(TL.compute_initial_q(μ, N))
    p₀     = zero(q₀)
    params = TL.default_parameters(T)
    hodeproblem(N, q₀, p₀; timespan = (T(t₀), T(t₁)), timestep = T(Δt), parameters = params)
end

# 2D projection for the trajectory plot: phase space of the first lattice site.
coords(sol) = (Float64.(vec(Array(sol.q)[1, :])), Float64.(vec(Array(sol.p)[1, :])))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "toda_lattice_energy_error_dt_$(Δt).png"),
    title = "Toda Lattice — Relative Energy Error (Δt = 0.1, t ≤ 100)")

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "toda_lattice_solution_error_dt_$(Δt).png"),
    title = "Toda Lattice — Solution Error (Δt = 0.1, t ≤ 100, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "toda_lattice_solution_dt_$(Δt).png"),
    title  = "Toda Lattice — Phase-Space Trajectory (Δt = 0.1, t ≤ 100)",
    coords = coords, xlabel = "q₁", ylabel = "p₁")
