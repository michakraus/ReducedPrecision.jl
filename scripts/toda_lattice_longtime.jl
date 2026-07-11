# Reduced-precision study: Toda lattice, long-time / coarse-step variant.
#
# Same pipeline as toda_lattice.jl (see there for why the initial conditions/parameters are
# hand-built at precision T and why a Hamiltonian closure is used), but with a coarse timestep
# Δt = 1 and a long horizon t = 10_000 (nt = 10_000 steps). The reference integration is guarded
# so the energy-error plot is still produced even if the high-order solve does not converge.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.TodaLattice: hodeproblem, hamiltonian
import GeometricProblems.TodaLattice as TL

const N  = 16       # lattice size
const μ  = 0.3      # bump-width parameter of the initial condition
const t₀ = 0.0
const Δt = 1.0
const nt = 10_000
const t₁ = nt * Δt

function make_problem(::Type{T}) where {T}
    q₀     = T.(TL.compute_initial_q(μ, N))
    p₀     = zero(q₀)
    params = map(T, TL.default_parameters)
    hodeproblem(N, q₀, p₀; timespan = (T(t₀), T(t₁)), timestep = T(Δt), parameters = params)
end

ham(t, q, p, params) = hamiltonian(t, q, p, params, N)
coords(sol) = (Float64.(vec(Array(sol.q)[1, :])), Float64.(vec(Array(sol.p)[1, :])))

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

plot_energy_error(runs, ham;
    path  = joinpath(plotdir, "toda_lattice_longtime_energy_error.png"),
    title = "Toda Lattice — Relative Energy Error (Δt = 1, t ≤ 10⁴)")

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = try
    integrate(make_problem(Float64), Gauss(8))
catch e
    @warn "reference integration failed; skipping solution-error and trajectory plots" error = sprint(showerror, e)
    nothing
end

if reference !== nothing
    plot_solution_error(runs, reference;
        path  = joinpath(plotdir, "toda_lattice_longtime_solution_error.png"),
        title = "Toda Lattice — Solution Error (Δt = 1, t ≤ 10⁴, vs. Float64 Gauss(8))")

    plot_solution(runs; reference = reference,
        path   = joinpath(plotdir, "toda_lattice_longtime_solution.png"),
        title  = "Toda Lattice — Phase-Space Trajectory (Δt = 1, t ≤ 10⁴)",
        coords = coords, xlabel = "q₁", ylabel = "p₁")
end
