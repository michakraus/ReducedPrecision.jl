# CairoMakie plotting of the study results.

_poslog(v) = (v > 0 && isfinite(v)) ? v : NaN
_finite_or_nan(v) = isfinite(v) ? Float64(v) : NaN

# Locate the run for a given method name and precision (nothing if it was not run).
function _find_run(runs, name, T)
    for run in runs
        run.precision === T && run.method.name == name && return run
    end
    return nothing
end

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

# Shared log-y limits for a whole figure: one range derived from every panel's finite,
# positive data, with the upper limit capped at 1e5 so runaway (non-geometric) errors are
# clipped at the top rather than dominating the scale.
function _shared_ylims(finite_pos)
    isempty(finite_pos) && return (1e-16, 1e0)
    lo, hi = extrema(finite_pos)
    loglo, loghi = log10(lo), log10(hi)
    pad = 0.05 * max(loghi - loglo, 1) + 0.1
    return (10.0^(loglo - pad), 10.0^min(loghi + pad, 5))
end

# Generic grid plot: one panel per precision, one line per method, log-scale y. Methods are
# drawn (and listed) in the order given by `methods`; every panel shares the same y-limits.
function _plot_grid(runs, seriesfun, ylabel, ptitle, path; methods = ALL_METHODS, precisions = PRECISIONS)
    colors = _method_colors()
    np = length(precisions)

    # First pass: collect each panel's (spec, t, y) series and the global finite range.
    finite_pos = Float64[]
    xr = nothing
    panels = Vector{Any}[]
    for T in precisions
        series = Any[]
        for spec in methods
            run = _find_run(runs, spec.name, T)
            (run === nothing || run.sol === nothing) && continue
            y = _poslog.(seriesfun(run))
            all(isnan, y) && continue
            append!(finite_pos, filter(isfinite, y))
            # Exact integration interval from the problem itself (0 and t₁ are exactly
            # representable in every precision), not the accumulated grid endpoint which
            # rounds slightly short/long at low precision.
            xr === nothing && (xr = Float64.(timespan(run.prob)))
            push!(series, (spec, timevalues(run.sol), y))
        end
        push!(panels, series)
    end

    yl = _shared_ylims(finite_pos)

    fig = _grid_figure(np, ptitle)
    for (j, T) in enumerate(precisions)
        ax = Axis(fig[1, j];
            yscale = log10,
            xlabel = "t",
            ylabel = j == 1 ? ylabel : "",
            title  = string(T),
        )
        for (spec, t, y) in panels[j]
            lines!(ax, t, y;
                color = colors[spec.name],
                linestyle = spec.geometric ? :solid : :dash,
                linewidth = 2,
            )
        end
        ylims!(ax, yl...)                                  # same y-limits for every panel
        xr !== nothing && xlims!(ax, xr[1], xr[2])         # exact time interval, no padding
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

# Grid of 2D trajectory plots: one panel per precision, one trajectory per method (drawn and
# listed in the order given by `methods`).
function _plot_trajectory_grid(runs, coordsfun, xlabel, ylabel, ptitle, path;
        methods = ALL_METHODS, precisions = PRECISIONS, reference = nothing)
    colors = _method_colors()
    np = length(precisions)
    lims = reference === nothing ? nothing : _traj_limits(reference, coordsfun)

    fig = _grid_figure(np, ptitle)

    for (j, T) in enumerate(precisions)
        ax = Axis(fig[1, j]; xlabel = xlabel, ylabel = j == 1 ? ylabel : "", title = string(T))
        # Reference trajectory drawn first, as a black backdrop, so every method lies on top.
        if reference !== nothing
            rx, ry = coordsfun(reference)
            lines!(ax, _finite_or_nan.(rx), _finite_or_nan.(ry);
                color = _REFERENCE_COLOR, linewidth = 2.5)
        end
        for spec in methods
            run = _find_run(runs, spec.name, T)
            (run === nothing || run.sol === nothing) && continue
            xs, ys = coordsfun(run.sol)
            xs = _finite_or_nan.(xs)
            ys = _finite_or_nan.(ys)
            (all(isnan, xs) || all(isnan, ys)) && continue
            lines!(ax, xs, ys;
                color = colors[spec.name],
                linestyle = spec.geometric ? :solid : :dash,
                linewidth = 1.5,
            )
        end
        if lims !== nothing
            xlims!(ax, lims[1]...)
            ylims!(ax, lims[2]...)
        end
    end

    extra = reference === nothing ? nothing :
            (LineElement(color = _REFERENCE_COLOR, linewidth = 2.5), "Reference")
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
