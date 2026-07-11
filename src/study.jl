# Running the method × precision sweep.

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
    integrate_bounded(problem, method; bound = 1e3) -> (sol, diverged)

Integrate `problem` with `method` step by step (replicating GeometricIntegrators' own stepping
loop, so results are identical to `integrate` for well-behaved runs) while guarding against
divergence: after each step the state `(q, p)` is checked, and if any component is non-finite or
exceeds `bound` in absolute value the integration stops. This avoids wasting steps on runs that
have already blown up — typically non-convergent implicit solves in low precision or at a coarse
timestep. The remaining steps are filled with `NaN` so downstream diagnostics/plots stop at the
divergence point. Returns the solution and the divergence step (`nothing` if the run stayed within
`bound`). Pass `bound = nothing` to disable the magnitude check (the non-finite check still fires).
"""
function integrate_bounded(problem, method; bound = 1e3)
    integrator = GeometricIntegrator(problem, method)
    sol = Solution(problem)
    solstep = solutionstep(integrator, sol[0])
    curstate = current(solstep)
    nt = ntime(sol)
    T = datatype(sol)

    diverged = nothing
    for n in 1:nt
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
    run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3)

Run every `method` at every `precision` on the problem produced by `make_problem(T)`.
The problem is built once per precision and reused across methods (it is immutable input to
`integrate`). Each integration is guarded by `integrate_bounded` (see there): divergent runs stop
early rather than producing runaway errors. Integration failures (e.g. a non-convergent implicit
solve that throws) are caught per run so a single failure does not abort the sweep. Returns a
`Vector{Run}`.
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS, bound = 1e3)
    runs = Run[]
    for T in precisions
        prob = make_problem(T)
        for spec in methods
            sol = nothing
            err = nothing
            diverged = nothing
            try
                sol, diverged = integrate_bounded(prob, spec.method; bound)
            catch e
                err = sprint(showerror, e)
                @warn "integration failed" method = spec.name precision = T
            end
            push!(runs, Run(spec, T, prob, sol, err, diverged))
        end
    end
    return runs
end
