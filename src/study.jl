# Running the method × precision sweep.
#
# ## Solver tolerances
#
# Nothing here scales the nonlinear solve's convergence criterion, because the framework already
# does — and the sweep depends on it. Convergence is `rfₐ ≤ f_abstol + f_reltol·‖F(x₀)‖`, and the
# better the initial guess, the smaller `‖F(x₀)‖` and the more the absolute term alone decides it.
# A floor fixed at `8eps(Float64) ≈ 1.8e-15` is then unreachable for a half-precision residual, which
# bottoms out near `eps(T)`. `GeometricIntegratorsBase.default_options` instead gives
#
#     f_abstol = max(8, solversize(method, problem)) · eps(datatype(problem))
#
# scaled to the working precision *and* to the stage-system size. This sweep's partitioned problems
# report `solversize = 0` and so sit on the `max(8, …)` floor at `8eps(T)`; the `Gauss(8)` references
# scale above it, which the degenerate 4D Lotka–Volterra one needs (32 unknowns, and a residual
# bottoming out just above `8eps(Float64)` for conditioning rather than precision reasons). Caller
# options are *merged into* `default_options`, so passing one does not drop `min_iterations`.
#
# `SimpleSolvers` supplies the complementary exit: a solve that stops making progress gives up after
# `max_stalls` steps rather than running to `max_iterations`, and a line search shares its solver's
# `Options` — which is what lets `verbosity = 0` (via `solveropts`) silence a configuration that is
# expected to fail.

"""
    capped_final_time(T, t₁, Δt)

Largest final integration time `≤ t₁` whose time grid `0, Δt, 2Δt, …` is still strictly increasing
in precision `T`. A reduced-precision grid *saturates* once the spacing `Δt` drops below the local
resolution (`ulp(t) ≥ Δt`): successive grid points then round to the same value (`t + Δt == t`).
This function returns the time just before the first such collision, or `t₁` if the grid stays
resolvable over the whole horizon (always the case for `Float32`/`Float64` at these horizons).

This is a **diagnostic**, not a constraint on the studies: it measures how far a `T`-typed *global*
clock can be carried before it stops advancing, and is therefore the motivation for the two things
that remove the limit — the local-frame stepping in [`integrate_bounded`](@ref) and the tableau-driven
initial guess in `initial_guess.jl`. The onset is Δt-dependent and much earlier for `BFloat16`
(8 significand bits) than for `Float16` (11):

| Δt   | `BFloat16` | `Float16` |
|:-----|:-----------|:----------|
| 0.01 | ≈ 2        | ≈ 16      |
| 0.1  | ≈ 16       | ≈ 128     |
| 1.0  | ≈ 256      | ≈ 2048    |
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
    reset_local!(solstep, Δt)

Prepare `solstep` for the next step in a **fixed local time frame**: shift the solution history back
by one slot (as [`GeometricBase.reset!`](@ref) does) and then label the times so that every step
looks like the *second* one — history at `0, Δt, …` and the target time at `nhistory · Δt`.

This is what lets the sweep carry arbitrarily long horizons in reduced precision. A `T`-typed global
clock saturates once `ulp(t) ≥ Δt` (see [`capped_final_time`](@ref)); the stored times then collide
and the implicit methods' `HermiteExtrapolation` throws `t₀ == t₁`. Re-anchoring near the origin
keeps consecutive times exactly `Δt` apart in every precision, because the step size the integrators
use comes from the problem (`timestep(problem)`) and never from a difference of clock values.

Two details the implementation depends on:

* `reset!` copies the stale time along with the state, so the whole history has to be relabelled
  after it — merely holding the target time fixed would collapse `t₀` onto `t₁`.
* The frame is anchored at **non-negative** times: `MidpointExtrapolation`'s sub-stepping loop for
  partitioned problems tests `abs(t + extrap.Δt) < abs(sol.t)`, which stays dormant for `t ≥ 0` at
  the step sizes used here but would silently engage around a negative history.

Scaling by an integer is exact in binary floating point, so `2Δt - Δt == Δt` and Hermite's
`s = (tᵢ - t₀)/Δt == 2` hold exactly — i.e. this reproduces what a non-saturating global clock does,
and is *more* accurate than one whose stored times have drifted.

Valid only for **autonomous** problems, where absolute time does not enter the vector field. That
holds for every problem in this study; `t` is an unused argument in all of their vector fields and
Hamiltonians.
"""
function reset_local!(solstep, Δt::T) where {T}
    N = nhistory(solstep)
    reset!(solstep, T(N) * Δt)
    for i in 1:N
        state(solstep, i).time .= T(N - i) * Δt
    end
    return solstep
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
is more robust than the line-search `Newton` in reduced precision; pass `solver = Newton()` to
compare); explicit methods carry no solver and ignore it.

`linesearch` selects the line-search method for the nonlinear solve (e.g. `Backtracking()`); the
default `nothing` leaves the solver's own default in place (`Backtracking` for `Newton`; `DogLeg`
is a trust-region method and ignores it, as do explicit methods). `max_iterations` caps the
nonlinear iterations per step (default `nothing` → the solver default of 1000); lowering it makes
hopeless implicit solves at a coarse timestep bail out early instead of churning, though it is
rarely the binding constraint, since a solve that stops making progress already gives up after
`max_stalls = 2` steps.

`solveropts` passes `SimpleSolvers.Options` keywords through to the nonlinear solve — e.g.
`(verbosity = 0,)` to silence a configuration that is *expected* to fail, or an `f_abstol` of your
own. It is empty by default, since the framework's tolerances are already precision- and size-scaled
(see the note at the top of this file), and overrides are merged into the method's `default_options`
rather than replacing them.

`initialguess` overrides the per-step initial guess (the extrapolation seeding the nonlinear solve)
for implicit methods; the default `nothing` uses the method's own default (`HermiteExtrapolation()`
for the Runge–Kutta / variational methods here). Pass `MidpointExtrapolation()` to seed each step by
integrating one history point forward with the vector field instead of the two-point Hermite
polynomial — this is more robust on stiff systems in low precision, where the Hermite guess can
produce `NaN` search directions. Ignored by explicit methods (they carry no solver).

`localclock` selects how the solution step's clock is advanced. The default `true` re-anchors it
every step via [`reset_local!`](@ref), which is what makes long horizons possible in reduced
precision. Pass `false` to advance it along the problem's own `T`-typed time grid, which saturates
at low precision and is retained only to demonstrate that failure.
"""
function integrate_bounded(problem, method; bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing, localclock = true, solveropts = (;))
    overrides = merge(
        solveropts,
        linesearch     === nothing ? (;) : (; linesearch),
        max_iterations === nothing ? (;) : (; max_iterations),
    )
    iguesskw = initialguess === nothing ? (;) : (; initialguess)
    integrator = if isimplicit(method) === true
        GeometricIntegrator(problem, method; solver, iguesskw..., overrides...)
    else
        GeometricIntegrator(problem, method)
    end
    sol = Solution(problem)
    solstep = solutionstep(integrator, sol[0])
    curstate = current(solstep)
    nt = ntime(sol)
    T = datatype(sol)
    Δt = timestep(problem)

    diverged = nothing
    for n in 1:nt
        # Advance the step's clock in a fixed local frame (see `reset_local!`) so it cannot
        # saturate; `localclock = false` walks the problem's own T-typed grid instead.
        localclock ? reset_local!(solstep, Δt) : reset!(solstep, timesteps(sol)[n])
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
    run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing, localclock = true)

Run every `method` at every `precision` on the problem produced by `make_problem(T)`.
The problem is built once per precision and reused across methods (it is immutable input to
`integrate`). Each integration is guarded by `integrate_bounded` (see there): divergent runs stop
early rather than producing runaway errors. Integration failures (e.g. a non-convergent implicit
solve that throws) are caught per run so a single failure does not abort the sweep. The implicit
methods use `solver` (and, if given, `linesearch` / `max_iterations` / `initialguess`) for the
nonlinear solve (see [`integrate_bounded`](@ref)). Returns a `Vector{Run}`.

`initialguess` may be an `Extrapolation` (or `nothing` for the method default) applied to every run,
or a **callable** `T -> Extrapolation | nothing` resolved per precision. The latter is how a study
uses a different initial guess only at the lowest precisions — e.g. `MidpointExtrapolation()` in
half precision (where it rescues stiff implicit solves that the default `HermiteExtrapolation` fails
on) and the method default elsewhere (where Midpoint would regress the multi-stage methods):

    initialguess = T -> T in (BFloat16, Float16) ? MidpointExtrapolation() : nothing

`solveropts` is resolved the same way — a `NamedTuple` of `SimpleSolvers.Options` keywords applied to
every run, or a callable `T -> NamedTuple` resolved per precision. See [`integrate_bounded`](@ref).
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing, localclock = true, solveropts = (;))
    runs = Run[]
    for T in precisions
        prob = make_problem(T)
        iguess = initialguess isa Function ? initialguess(T) : initialguess
        opts = solveropts isa Function ? solveropts(T) : solveropts
        for spec in methods
            sol = nothing
            err = nothing
            diverged = nothing
            try
                sol, diverged = integrate_bounded(prob, spec.method; bound, solver, linesearch, max_iterations, initialguess = iguess, localclock, solveropts = opts)
            catch e
                err = sprint(showerror, e)
                @warn "integration failed" method = spec.name precision = nameof(T)
            end
            push!(runs, Run(spec, T, prob, sol, err, diverged))
        end
    end
    return runs
end
