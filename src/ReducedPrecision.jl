"""
    ReducedPrecision

Reusable pipeline for studying geometric (symplectic) vs non-geometric integrators
run in varying floating-point precision (Float16, Float32, Float64).

The study is driven by a per-problem `make_problem(T)` closure that must return a
partitioned problem (PODE/HODE) at precision `T` with **both** initial conditions and
timespan/timestep in `T` (so `datatype(prob) === timetype(prob) === T`). All methods in
the registry are then run across every precision, and energy-error / solution-error plots
are produced with CairoMakie.

Symplectic Euler (A/B) require a partitioned problem; the non-geometric explicit/implicit
Euler methods are represented by their Runge–Kutta tableau twins (`ExplicitEulerRK` /
`ImplicitEulerRK`), which auto-promote to partitioned RK on a PODE/HODE — so a single
problem form runs the whole method set for a fair comparison.
"""
module ReducedPrecision

using GeometricIntegrators
using GeometricSolutions
using GeometricBase
using GeometricEquations: parameters
using CairoMakie

export PRECISIONS, MethodSpec, GEOMETRIC_METHODS, NONGEOMETRIC_METHODS, ALL_METHODS
export EULER_METHODS, OTHER_METHODS, METHOD_GROUPS
export Run, run_study, assert_precision, verify_precision
export energy_error, solution_error, timevalues
export plot_energy_error, plot_solution_error, plot_solution

const PRECISIONS = (Float16, Float32, Float64)

"""
    MethodSpec(name, method, geometric)

A named integration method together with a flag marking it as geometric (symplectic) or not.
"""
struct MethodSpec
    name::String
    method::GeometricMethod
    geometric::Bool
end

const GEOMETRIC_METHODS = MethodSpec[
    MethodSpec("Symplectic Euler A", SymplecticEulerA(), true),
    MethodSpec("Symplectic Euler B", SymplecticEulerB(), true),
    MethodSpec("Implicit Midpoint",  ImplicitMidpoint(),  true),
]

const NONGEOMETRIC_METHODS = MethodSpec[
    MethodSpec("Explicit Euler",    ExplicitEulerRK(),  false),
    MethodSpec("Implicit Euler",    ImplicitEulerRK(),  false),
    MethodSpec("Explicit Midpoint", ExplicitMidpoint(), false),
    MethodSpec("Crank-Nicolson",    CrankNicolson(),    false),
    MethodSpec("RK4",               RK4(),              false),
]

const ALL_METHODS = vcat(GEOMETRIC_METHODS, NONGEOMETRIC_METHODS)

# Grouping used for plotting (orthogonal to the geometric/non-geometric classification, which
# still drives the line style). Plots are produced once per group so each figure stays readable.
const EULER_METHODS = filter(
    m -> m.name in ("Symplectic Euler A", "Symplectic Euler B", "Explicit Euler", "Implicit Euler"),
    ALL_METHODS)
const OTHER_METHODS = filter(
    m -> m.name in ("Implicit Midpoint", "Explicit Midpoint", "Crank-Nicolson", "RK4"),
    ALL_METHODS)
const METHOD_GROUPS = ["euler" => EULER_METHODS, "other" => OTHER_METHODS]

_group_title(label) = label == "euler" ? "Euler methods" : "other methods"

"""
    Run

Result of a single (method, precision) integration. `sol` is the `GeometricSolution` on
success, or `nothing` if the integration threw (with the message stored in `error`).
"""
struct Run
    method::MethodSpec
    precision::DataType
    prob::Any
    sol::Any
    error::Union{Nothing,String}
end

"""
    run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS)

Run every `method` at every `precision` on the problem produced by `make_problem(T)`.
Integration failures (e.g. Float16 blow-ups / non-convergent implicit solves) are caught
per run so a single failure does not abort the sweep. Returns a `Vector{Run}`.
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS)
    runs = Run[]
    for T in precisions
        for spec in methods
            prob = make_problem(T)
            sol = nothing
            err = nothing
            try
                sol = integrate(prob, spec.method)
            catch e
                err = sprint(showerror, e)
                @warn "integration failed" method=spec.name precision=T
            end
            push!(runs, Run(spec, T, prob, sol, err))
        end
    end
    return runs
end

"""
    assert_precision(prob, sol, T)

Verification gate: assert that neither the problem nor the resulting solution silently
promoted to another precision. Checks the problem's and solution's `datatype`/`timetype`
and the element types of the stored `q`, `p` arrays all equal `T`. Throws an
`AssertionError` on any mismatch; returns `true` otherwise.
"""
function assert_precision(prob, sol, ::Type{T}) where {T}
    @assert datatype(prob) === T "datatype(prob) = $(datatype(prob)) ≠ $T"
    @assert timetype(prob) === T "timetype(prob) = $(timetype(prob)) ≠ $T"
    @assert datatype(sol)  === T "datatype(sol) = $(datatype(sol)) ≠ $T"
    @assert timetype(sol)  === T "timetype(sol) = $(timetype(sol)) ≠ $T"
    @assert eltype(Array(sol.q)) === T "eltype(sol.q) = $(eltype(Array(sol.q))) ≠ $T"
    @assert eltype(Array(sol.p)) === T "eltype(sol.p) = $(eltype(Array(sol.p))) ≠ $T"
    return true
end

"""
    verify_precision(runs)

Run `assert_precision` on every successful run in `runs` and print a per-run report
(failed integrations are reported and skipped). Throws on the first purity violation.
This is the gate proving no library silently promotes the requested precision to Float64.
"""
function verify_precision(runs)
    println("Verifying precision purity (datatype == timetype == requested precision):")
    for run in runs
        if run.sol === nothing
            println("  [skip] $(run.method.name) @ $(run.precision): integration failed ($(run.error))")
            continue
        end
        assert_precision(run.prob, run.sol, run.precision)
        println("  [ ok ] $(run.method.name) @ $(run.precision)")
    end
    println("Precision verification passed for all successful runs.")
    return nothing
end

"""
    timevalues(sol) -> Vector{Float64}

Nominal time grid `t₀ + n·Δt` of the solution (including the initial time), as `Float64`
for plotting. The nominal grid is used rather than the stored clock `sol.t` because at low
precision and long horizons the stored clock saturates (e.g. in Float16, `t + Δt == t` once
`t` exceeds the representable integer range), which would corrupt the time axis.
"""
function timevalues(sol)
    t₀ = Float64(sol.t[0])
    Δt = Float64(timestep(sol))
    return Float64[t₀ + n * Δt for n in 0:ntime(sol)]
end

"""
    energy_error(sol, hamiltonian) -> Vector{Float64}

Relative energy error `|(H(qₙ,pₙ) − H₀) / H₀|` at each time step. The Hamiltonian is
evaluated in the solution's own precision (via `compute_invariant_error`); the absolute
value is returned as `Float64` for log-scale plotting.
"""
function energy_error(sol, hamiltonian)
    params = parameters(sol.problem)
    _, errds = compute_invariant_error(sol.t, sol.q, sol.p, params, hamiltonian)
    return Float64[abs(Float64(errds[n])) for n in 0:ntime(sol)]
end

"""
    solution_error(sol, reference) -> Vector{Float64}

Euclidean norm of the state error `‖(qₙ,pₙ) − (qₙ,pₙ)_ref‖` at each time step, comparing
index by index against `reference` (which must share the same time grid). Both solutions
are materialized to `Float64` before comparison, so `reference` may be at a different
precision than `sol` (e.g. a Float64 high-order reference).
"""
function solution_error(sol, reference)
    q  = Float64.(Array(sol.q));       p  = Float64.(Array(sol.p))
    qr = Float64.(Array(reference.q)); pr = Float64.(Array(reference.p))
    @assert size(q, 2) == size(qr, 2) "solution and reference have different lengths"
    n = size(q, 2)
    return Float64[
        sqrt(sum(abs2, @view(q[:, i]) .- @view(qr[:, i])) +
             sum(abs2, @view(p[:, i]) .- @view(pr[:, i]))) for i in 1:n
    ]
end


# --- plotting -------------------------------------------------------------------------

_poslog(v) = (v > 0 && isfinite(v)) ? v : NaN
_finite_or_nan(v) = isfinite(v) ? Float64(v) : NaN

# Wong 8-colour, colour-blind-safe palette (one colour per method in ALL_METHODS).
const _PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00",
                  "#56B4E9", "#CC79A7", "#F0E442", "#999999"]
const _REFERENCE_COLOR = :black

function _method_colors()
    Dict(spec.name => _PALETTE[mod1(i, length(_PALETTE))]
         for (i, spec) in enumerate(ALL_METHODS))
end

# Insert a suffix before the file extension, e.g. "foo.png" -> "foo_euler.png".
function _suffix_path(path, suffix)
    base, ext = splitext(path)
    return string(base, "_", suffix, ext)
end

# Legend as a horizontal row below all panels (spanning every column). `extra` optionally
# appends one more (element, label) entry, e.g. the reference trajectory.
function _legend_below!(fig, methods, colors, np; extra = nothing)
    elems  = LineElement[LineElement(color = colors[s.name],
                          linestyle = s.geometric ? :solid : :dash,
                          linewidth = 2) for s in methods]
    labels = String[s.name for s in methods]
    if extra !== nothing
        push!(elems, extra[1])
        push!(labels, extra[2])
    end
    Legend(fig[2, 1:np], elems, labels, "Method";
        orientation = :horizontal, framevisible = true, titleposition = :left)
    return nothing
end

# Figure width holds the per-panel width fixed (~430 px) and drops the old side-legend
# column; the extra height leaves room for the legend row below.
_grid_figure(np, ptitle) = begin
    fig = Figure(size = (430 * np, 500), fontsize = 14)
    Label(fig[0, 1:np], ptitle; fontsize = 18, font = :bold)
    fig
end

# Generic grid plot: one panel per precision, one line per method, log-scale y.
function _plot_grid(runs, seriesfun, ylabel, ptitle, path; methods = ALL_METHODS, precisions = PRECISIONS)
    colors = _method_colors()
    wanted = Set(m.name for m in methods)
    np = length(precisions)

    fig = _grid_figure(np, ptitle)

    for (j, T) in enumerate(precisions)
        ax = Axis(fig[1, j];
            yscale = log10,
            xlabel = "t",
            ylabel = j == 1 ? ylabel : "",
            title  = string(T),
        )
        finite_pos = Float64[]
        xr = nothing
        for run in runs
            run.precision === T || continue
            run.method.name in wanted || continue
            run.sol === nothing && continue
            y = _poslog.(seriesfun(run))
            all(isnan, y) && continue
            append!(finite_pos, filter(isfinite, y))
            # Exact integration interval from the problem itself (0 and t₁ are exactly
            # representable in every precision), not the accumulated grid endpoint which
            # rounds slightly short/long at low precision.
            xr === nothing && (xr = Float64.(timespan(run.prob)))
            lines!(ax, timevalues(run.sol), y;
                color = colors[run.method.name],
                linestyle = run.method.geometric ? :solid : :dash,
                linewidth = 2,
            )
        end
        # Explicit finite y-limits, with the upper limit capped at 1e5 so runaway
        # (non-geometric) errors are clipped at the top rather than dominating the scale.
        if isempty(finite_pos)
            ylims!(ax, 1e-16, 1e0)
        else
            lo, hi = extrema(finite_pos)
            loglo, loghi = log10(lo), log10(hi)
            pad = 0.05 * max(loghi - loglo, 1) + 0.1
            ylims!(ax, 10.0^(loglo - pad), 10.0^min(loghi + pad, 5))
        end
        # Fit the x-axis exactly to the time interval (no padding).
        xr !== nothing && xlims!(ax, xr[1], xr[2])
    end

    _legend_below!(fig, methods, colors, np)

    mkpath(dirname(path))
    save(path, fig)
    @info "saved figure" path
    return fig
end

# Default 2D coordinates for the solution plot: phase space (q, p) for a one-degree-of-freedom
# system, configuration space (q₁, q₂) for a multi-degree-of-freedom system.
function _default_coords(sol)
    q = Float64.(Array(sol.q))
    p = Float64.(Array(sol.p))
    xs = vec(q[1, :])
    ys = size(q, 1) == 1 ? vec(p[1, :]) : vec(q[2, :])
    return (xs, ys)
end

# Axis limits (with a little padding) taken from a bounded reference trajectory, so runaway
# methods are clipped to the region of the true solution rather than dominating the scale.
function _traj_limits(reference, coordsfun)
    xs, ys = coordsfun(reference)
    xf = filter(isfinite, xs)
    yf = filter(isfinite, ys)
    (isempty(xf) || isempty(yf)) && return nothing
    xlo, xhi = extrema(xf)
    ylo, yhi = extrema(yf)
    px = 0.08 * (xhi - xlo) + eps()
    py = 0.08 * (yhi - ylo) + eps()
    return ((xlo - px, xhi + px), (ylo - py, yhi + py))
end

# Grid of 2D trajectory plots: one panel per precision, one trajectory per method.
function _plot_trajectory_grid(runs, coordsfun, xlabel, ylabel, ptitle, path;
        methods = ALL_METHODS, precisions = PRECISIONS, reference = nothing)
    colors = _method_colors()
    wanted = Set(m.name for m in methods)
    np = length(precisions)
    lims = reference === nothing ? nothing : _traj_limits(reference, coordsfun)

    fig = _grid_figure(np, ptitle)

    for (j, T) in enumerate(precisions)
        ax = Axis(fig[1, j]; xlabel = xlabel, ylabel = j == 1 ? ylabel : "", title = string(T))
        # Reference trajectory as a black backdrop, so each method's deviation is visible.
        if reference !== nothing
            rx, ry = coordsfun(reference)
            lines!(ax, _finite_or_nan.(rx), _finite_or_nan.(ry);
                color = _REFERENCE_COLOR, linewidth = 2.5)
        end
        for run in runs
            run.precision === T || continue
            run.method.name in wanted || continue
            run.sol === nothing && continue
            xs, ys = coordsfun(run.sol)
            xs = _finite_or_nan.(xs)
            ys = _finite_or_nan.(ys)
            (all(isnan, xs) || all(isnan, ys)) && continue
            lines!(ax, xs, ys;
                color = colors[run.method.name],
                linestyle = run.method.geometric ? :solid : :dash,
                linewidth = 1.5,
            )
        end
        if lims !== nothing
            xlims!(ax, lims[1]...)
            ylims!(ax, lims[2]...)
        end
    end

    extra = reference === nothing ? nothing :
            (LineElement(color = _REFERENCE_COLOR, linewidth = 2.5), "reference")
    _legend_below!(fig, methods, colors, np; extra)

    mkpath(dirname(path))
    save(path, fig)
    @info "saved figure" path
    return fig
end

"""
    plot_energy_error(runs, hamiltonian; path, title)

Plot the relative energy error over time, one panel per precision, methods overlaid
(geometric = solid, non-geometric = dashed). Two figures are written — one for the Euler
methods and one for the remaining methods — with `_euler` / `_other` appended to `path`.
"""
function plot_energy_error(runs, hamiltonian; path, title)
    figs = Any[]
    for (label, methods) in METHOD_GROUPS
        push!(figs, _plot_grid(runs, run -> energy_error(run.sol, hamiltonian),
            "|ΔH / H₀|", "$title — $(_group_title(label))", _suffix_path(path, label); methods))
    end
    return figs
end

"""
    plot_solution_error(runs, reference; path, title)

Plot the state solution error over time versus `reference`, one panel per precision. As with
`plot_energy_error`, two figures are written (`_euler` / `_other`).
"""
function plot_solution_error(runs, reference; path, title)
    figs = Any[]
    for (label, methods) in METHOD_GROUPS
        push!(figs, _plot_grid(runs, run -> solution_error(run.sol, reference),
            "‖x − x_ref‖", "$title — $(_group_title(label))", _suffix_path(path, label); methods))
    end
    return figs
end

"""
    plot_solution(runs; path, title, reference = nothing, coords = _default_coords,
                  xlabel = "q", ylabel = "p")

Plot the 2D trajectory of each method's solution, one panel per precision (geometric = solid,
non-geometric = dashed). `coords(sol)` returns the `(xs, ys)` to plot — by default phase space
`(q, p)` for a one-degree-of-freedom system and configuration space `(q₁, q₂)` otherwise. When
`reference` is given, the axes are fitted to the reference trajectory so runaway methods are
clipped. Two figures are written (`_euler` / `_other`).
"""
function plot_solution(runs; path, title, reference = nothing,
        coords = _default_coords, xlabel = "q", ylabel = "p")
    figs = Any[]
    for (label, methods) in METHOD_GROUPS
        push!(figs, _plot_trajectory_grid(runs, coords, xlabel, ylabel,
            "$title — $(_group_title(label))", _suffix_path(path, label); methods, reference))
    end
    return figs
end

end # module
