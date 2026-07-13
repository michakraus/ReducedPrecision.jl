# Reduced-precision study: double pendulum.
#
# The DoublePendulum module (EulerLagrange-generated) exposes only `hodeproblem`/`lodeproblem`
# and has NO `::Type{T}` precision parameter, so `make_problem(T)` hand-builds T-typed initial
# conditions, timespan, timestep and parameters. Chaotic system with no closed-form solution;
# the solution error is measured against a Float64 Gauss(8) reference on the same grid.
#
# The timespan/timestep are set here (matching the DoublePendulum defaults, t ≤ 10 at Δt = 0.01)
# rather than read from the module, so the script does not depend on the module's timespan/timestep
# symbol names (only θ₀/p₀/default_parameters are pulled from the module).
#
# The implicit solves use the line-search `Newton` solver with a `Backtracking` linesearch and a
# 100-iteration cap (matching the coarse-step variant), instead of the default trust-region `DogLeg`.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.DoublePendulum: hodeproblem, hamiltonian
import GeometricProblems.DoublePendulum as DP

const Δt = 0.01
const t₀ = 0.0
const t₁ = 10.0

function make_problem(::Type{T}) where {T}
    q₀     = T.(DP.θ₀)
    p₀     = T.(DP.p₀)
    tspan  = (T(t₀), T(t₁))
    dt     = T(Δt)
    params = DP.default_parameters(T)
    hodeproblem(q₀, p₀; timespan = tspan, timestep = dt, parameters = params)
end

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem; solver = Newton(), linesearch = Backtracking(), max_iterations = 100)

verify_precision(runs)

# high-precision reference (Float64, high-order symplectic, same time grid)
reference = integrate(make_problem(Float64), Gauss(8))

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "double_pendulum_energy_error_dt_$(Δt).png"),
    title = "Double Pendulum — Relative Energy Error (Δt = 0.01, t ≤ 10)")

plot_solution_error(runs, reference;
    path  = joinpath(plotdir, "double_pendulum_solution_error_dt_$(Δt).png"),
    title = "Double Pendulum — Solution Error (Δt = 0.01, t ≤ 10, vs. Float64 Gauss(8))")

plot_solution(runs; reference = reference,
    path   = joinpath(plotdir, "double_pendulum_solution_dt_$(Δt).png"),
    title  = "Double Pendulum — Configuration-Space Trajectory (Δt = 0.01, t ≤ 10)",
    xlabel = "θ₁", ylabel = "θ₂")
