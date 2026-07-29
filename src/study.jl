# Running the method × precision sweep.

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
    solver_tolerances(T)

Nonlinear-solver convergence tolerance scaled to the working precision `T`, as a `NamedTuple` of
`SimpleSolvers.Options` keywords.

`SimpleSolvers` already builds its `Options` at the problem's precision, so `x_abstol = 2eps(T)`,
`f_reltol = √eps(T)` and friends come out correctly scaled on their own. The exception is the
*absolute* residual tolerance: `GeometricIntegratorsBase.default_options` hard-codes
`f_abstol = 8eps()`, i.e. `8eps(Float64) ≈ 1.8e-15`, whatever the precision of the run — overriding
SimpleSolvers' own `absolute_tolerance(T) = zero(T)`.

That matters because convergence is assessed as `rfₐ ≤ f_abstol + f_reltol · ‖F(x₀)‖`. The better the
initial guess, the smaller `‖F(x₀)‖` and the more the *absolute* term decides the test — so with a
good guess and a `1.8e-15` floor a half-precision solve, whose residual bottoms out near `eps(T)`,
can never satisfy it and exhausts `max_iterations` (1000) on every step. That is slow, floods the log
with trust-region warnings, and is meaningless as a criterion: it asks a `BFloat16` solve to deliver a
`Float64` residual.

Scaling it as `8eps(T)` reproduces the current value exactly at `Float64` and lets the lower
precisions stop once they have converged as far as their arithmetic allows.

The two-argument form additionally scales by the size `n` of the stage system, as `8·√n·eps(T)`: the
residual is measured as an `l2` norm over `n` components, whose own round-off grows like `√n`, so an
`n`-independent absolute floor is dimensionally wrong. This matters for the high-stage-count reference
integrations — see [`reference_solution`](@ref).
"""
solver_tolerances(::Type{T}) where {T} = (f_abstol = 8 * Float64(eps(T)),)
solver_tolerances(::Type{T}, n::Integer) where {T} = (f_abstol = 8 * sqrt(n) * Float64(eps(T)),)

"""
    reference_solution(problem, method; kwargs...)

Integrate `problem` with `method` and return the solution, for use as the high-accuracy reference that
a study measures its solution error against.

This is `integrate(problem, method)` with the nonlinear solve's absolute residual tolerance scaled to
the size of the stage system via [`solver_tolerances`](@ref), rather than left at the stack's fixed
`8eps()`.

The references here are high-stage-count rules (`Gauss(8)`), so their stage systems are large: 16
unknowns for the pendulum, 32 for the double pendulum and the 4D Lotka–Volterra system, 256 for the
16-site Toda lattice. Most of them still reach `8eps(Float64) ≈ 1.8e-15` — but the 4D Lotka–Volterra
system does not. It is degenerate, posed with the quasi-canonical reduced gauge matrix, and its
residual bottoms out just above that threshold, so the solve exhausted all 1000 iterations on *every*
step. Loosening the floor by the ~5.7× that `√32` supplies converges it in a handful of iterations
instead: measured over 20 steps, 4 capped solves and 0.59 s become 0 and 0.002 s.

Two things made this worth fixing rather than tolerating. It dominated the cost of a full
`run_all.jl`, and — because a non-convergent solve is only reported as a warning — a *reference*
solution was being used without having actually met its convergence criterion.
"""
function reference_solution(problem, method; kwargs...)
    isimplicit(method) === true || return integrate(problem, method; kwargs...)
    n = length(nlsolution(GeometricIntegrator(problem, method)))
    opts = solver_tolerances(datatype(problem), n)
    return integrate(problem, method; default_options(method)..., opts..., kwargs...)
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
hopeless implicit solves at a coarse timestep bail out early instead of churning.

The solve's absolute residual tolerance comes from [`solver_tolerances`](@ref) at the problem's
precision, which is what stops the half-precision solves from exhausting `max_iterations` on an
unreachable `Float64` residual. Pass `solveropts` to override it (or anything else
`SimpleSolvers.Options` accepts). All overrides are forwarded alongside the integrator's
`default_options`, so `min_iterations` is preserved — passing any option keyword otherwise replaces
the whole default set.

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
function integrate_bounded(problem, method; bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing, localclock = true, solveropts = solver_tolerances(datatype(problem)))
    overrides = merge(
        solveropts,
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

`solveropts` is resolved the same way: by default the callable [`solver_tolerances`](@ref), applied
per precision. Pass a `NamedTuple` to use one fixed set of `SimpleSolvers.Options` keywords for every
run instead.
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3, solver = DogLeg(), linesearch = nothing, max_iterations = nothing, initialguess = nothing, localclock = true, solveropts = solver_tolerances)
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
