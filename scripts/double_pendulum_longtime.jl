# Reduced-precision study: double pendulum, coarse-step variant.
#
# Same pipeline as double_pendulum.jl (hand-built T-typed initial conditions/parameters,
# since DoublePendulum has no ::Type{T} constructor), but with a coarser timestep Δt = 0.1 and
# a longer horizon t = 1000 (nt = 10_000 steps). Δt = 0.1 is already coarse for this stiff,
# chaotic system, so several runs (and possibly the reference) may fail — all are guarded. The
# solution error is measured against a Float64 Gauss(8) reference computed at the *fine* step
# (Δt_ref = 0.01, as in the short scenario) and subsampled onto the coarse grid.

using ReducedPrecision
using GeometricIntegrators: Gauss, integrate
using GeometricProblems.DoublePendulum: hodeproblem, hamiltonian
import GeometricProblems.DoublePendulum as DP

const t₀ = 0.0
const Δt = 0.1
const nt = 10_000
const t₁ = nt * Δt
const Δt_ref = 0.01  # fine reference step (matches the short scenario)

function _dp_problem(::Type{T}, dt) where {T}
    hodeproblem(T.(DP.θ₀), T.(DP.p₀);
        timespan = (T(t₀), T(t₁)), timestep = T(dt), parameters = map(T, DP.default_parameters))
end

make_problem(::Type{T})   where {T} = _dp_problem(T, Δt)
make_reference(::Type{T}) where {T} = _dp_problem(T, Δt_ref)

const plotdir = normpath(joinpath(@__DIR__, "..", "plots"))

runs = run_study(make_problem)

verify_precision(runs)

plot_energy_error(runs, hamiltonian;
    path  = joinpath(plotdir, "double_pendulum_energy_error_dt_$(Δt).png"),
    title = "Double Pendulum — Relative Energy Error (Δt = 0.1, t ≤ 10³)")

# high-precision reference (Float64, high-order symplectic, fine step, subsampled to the grid)
reference = try
    integrate(make_reference(Float64), Gauss(8))
catch e
    @warn "reference integration failed; skipping solution-error plot" error = sprint(showerror, e)
    nothing
end

if reference !== nothing
    plot_solution_error(runs, reference;
        path  = joinpath(plotdir, "double_pendulum_solution_error_dt_$(Δt).png"),
        title = "Double Pendulum — Solution Error (Δt = 0.1, t ≤ 10³, vs. Float64 Gauss(8) at Δt = 0.01)")

    plot_solution(runs; reference = reference,
        path   = joinpath(plotdir, "double_pendulum_solution_dt_$(Δt).png"),
        title  = "Double Pendulum — Configuration-Space Trajectory (Δt = 0.1, t ≤ 10³)",
        xlabel = "θ₁", ylabel = "θ₂")
end
