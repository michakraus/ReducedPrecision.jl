# Running the method × precision sweep.

"""
    capped_final_time(T, t₁, Δt)

Largest final integration time `≤ t₁` whose time grid `0, Δt, 2Δt, …` is still strictly increasing
in precision `T`. A reduced-precision grid *saturates* once the spacing `Δt` drops below the local
resolution (`ulp(t) ≥ Δt`): successive grid points then round to the same value (`t + Δt == t`),
which stalls the integration and makes the implicit methods' Hermite initial guess throw
(`t₀ == t₁`). The saturation onset is **Δt-dependent** — in `Float16`, `t ≈ 2048` for `Δt = 1` but
already `t ≈ 128` for `Δt = 0.1` — so the cap is derived from the actual `T`-grid rather than a fixed
constant: it returns the time just before the first collision, or `t₁` if the grid stays resolvable
over the whole horizon (always the case for `Float32`/`Float64` at these horizons). Use it in a
`make_problem(T)` closure:

    make_problem(::Type{T}) where {T} =
        podeproblem(...; timespan = (T(t₀), T(capped_final_time(T, t₁, Δt))), timestep = T(Δt))
"""
function capped_final_time(::Type{T}, t₁, Δt) where {T}
    # `GeometricSolutions.TimeSeries` backs the grid with the range `tbegin:Δt:tend`, so mirror that
    # construction. First find where the full-horizon grid stops advancing (`grid[i] == grid[i-1]`,
    # the saturation point).
    step = T(Δt)
    tbeg = T(zero(t₁))
    full = tbeg:step:T(t₁)
    icut = length(full)
    for i in 2:length(full)
        full[i] <= full[i-1] && (icut = i - 1; break)
    end
    # The solution rebuilds the grid as `tbeg:step:tend`; a `StepRangeLen` pins its endpoint, which
    # can re-collide with its predecessor even though the same value sits cleanly *inside* the
    # full-horizon range. Return the largest endpoint whose rebuilt grid is strictly increasing —
    # everything at or below `icut` already has a collision-free interior, so only the endpoint is
    # in question (at most a step or two back).
    for i in icut:-1:1
        g = tbeg:step:full[i]
        (length(g) < 2 || g[end] > g[end-1]) && return Float64(full[i])
    end
    return Float64(tbeg)
end

"""
    Run

Result of a single (method, precision) integration.

* `sol` is the `GeometricSolution` on success, or `nothing` if the integration threw (with the
  message in `error`).
* `diverged` is the step index at which the divergence guard tripped (the state magnitude
  exceeded the bound or became non-finite), or `nothing` if the run stayed bounded. Steps after
  `diverged` are filled with `NaN`, so the error metrics and plots truncate at the blow-up.
"""
struct Run
    method::MethodSpec
    precision::DataType
    prob::Any
    sol::Any
    error::Union{Nothing,String}
    diverged::Union{Nothing,Int}
end

"""
    integrate_bounded(problem, method; bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing) -> (sol, diverged)

Integrate `problem` with `method` step by step (replicating GeometricIntegrators' own stepping
loop, so results are identical to `integrate` for well-behaved runs) while guarding against
divergence: after each step the state `(q, p)` is checked, and if any component is non-finite or
exceeds `bound` in absolute value the integration stops. This avoids wasting steps on runs that
have already blown up — typically non-convergent implicit solves in low precision or at a coarse
timestep. The remaining steps are filled with `NaN` so downstream diagnostics/plots stop at the
divergence point. Returns the solution and the divergence step (`nothing` if the run stayed within
`bound`). Pass `bound = nothing` to disable the magnitude check (the non-finite check still fires).

For implicit methods the nonlinear solve uses `solver` (default the trust-region `DogLeg`, which
is more robust than the line-search `Newton` in reduced precision); explicit methods carry no
solver and ignore it. Pass `solver = Newton()` to reproduce the previous behaviour.

`linesearch` selects the line-search method for the nonlinear solve (e.g. `Backtracking()`); the
default `nothing` leaves the solver's own default in place (`Backtracking` for `Newton`; `DogLeg`
is a trust-region method and ignores it, as do explicit methods). `max_iterations` caps the
nonlinear iterations per step (default `nothing` → the solver default of 1000); lowering it makes
hopeless implicit solves at a coarse timestep bail out early instead of churning. Any override is
forwarded alongside the integrator's `default_options` so `min_iterations`/`f_abstol` are preserved
(passing any option kwarg otherwise replaces the whole default-option set).

`initialguess` overrides the per-step initial guess (the extrapolation seeding the nonlinear solve)
for implicit methods; the default `nothing` uses the method's own default (`HermiteExtrapolation()`
for the Runge–Kutta / variational methods here). Pass `MidpointExtrapolation()` to seed each step by
integrating one history point forward with the vector field instead of the two-point Hermite
polynomial — this is more robust when the low-precision time grid saturates (`t₀ == t₁`), which
breaks the Hermite guess. Ignored by explicit methods (they carry no solver).
"""
function integrate_bounded(problem, method; bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing)
    overrides = merge(
        linesearch     === nothing ? (;) : (; linesearch),
        max_iterations === nothing ? (;) : (; max_iterations),
    )
    iguesskw = initialguess === nothing ? (;) : (; initialguess)
    integrator = if isimplicit(method) === true
        isempty(overrides) ?
            GeometricIntegrator(problem, method; solver, iguesskw...) :
            GeometricIntegrator(problem, method; solver, iguesskw..., default_options(method)..., overrides...)
    else
        GeometricIntegrator(problem, method)
    end
    sol = Solution(problem)
    solstep = solutionstep(integrator, sol[0])
    curstate = current(solstep)
    nt = ntime(sol)
    T = datatype(sol)

    diverged = nothing
    for n in 1:nt
        reset!(solstep, timesteps(sol)[n])
        integrate!(solstep, integrator)     # advance one step
        copy!(sol, curstate, n)             # store it into the solution
        q = sol.q[n]
        p = sol.p[n]
        if !(all(isfinite, q) && all(isfinite, p)) ||
           (bound !== nothing && (maximum(abs, q) > bound || maximum(abs, p) > bound))
            diverged = n
            break
        end
    end

    # Blank out the (un-integrated) tail so energy_error / solution_error / the trajectory plots
    # break cleanly at the divergence point instead of reading leftover zeros.
    if diverged !== nothing
        for n in (diverged + 1):nt
            fill!(sol.q[n], T(NaN))
            fill!(sol.p[n], T(NaN))
        end
    end

    return sol, diverged
end

"""
    run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing)

Run every `method` at every `precision` on the problem produced by `make_problem(T)`.
The problem is built once per precision and reused across methods (it is immutable input to
`integrate`). Each integration is guarded by `integrate_bounded` (see there): divergent runs stop
early rather than producing runaway errors. Integration failures (e.g. a non-convergent implicit
solve that throws) are caught per run so a single failure does not abort the sweep. The implicit
methods use `solver` (and, if given, `linesearch` / `max_iterations` / `initialguess`) for the
nonlinear solve (see [`integrate_bounded`](@ref)). Returns a `Vector{Run}`.

`initialguess` may be an `Extrapolation` (or `nothing` for the method default) applied to every run,
or a **callable** `T -> Extrapolation | nothing` resolved per precision. The latter is how a study
uses a different initial guess only at a specific precision — e.g. `MidpointExtrapolation()` in
`Float16` (where it rescues stiff implicit solves that the default `HermiteExtrapolation` fails on)
and the method default elsewhere (where Midpoint would regress the multi-stage methods):

    initialguess = T -> T === Float16 ? MidpointExtrapolation() : nothing
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing)
    runs = Run[]
    for T in precisions
        prob = make_problem(T)
        iguess = initialguess isa Function ? initialguess(T) : initialguess
        for spec in methods
            sol = nothing
            err = nothing
            diverged = nothing
            try
                sol, diverged = integrate_bounded(prob, spec.method; bound, solver, linesearch, max_iterations, initialguess = iguess)
            catch e
                err = sprint(showerror, e)
                @warn "integration failed" method = spec.name precision = T
            end
            push!(runs, Run(spec, T, prob, sol, err, diverged))
        end
    end
    return runs
end
