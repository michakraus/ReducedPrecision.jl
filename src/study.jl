# Running the method × precision sweep.

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
The problem is built once per precision and reused across methods (it is immutable input to
`integrate`). Integration failures (e.g. Float16 blow-ups / non-convergent implicit solves)
are caught per run so a single failure does not abort the sweep. Returns a `Vector{Run}`.
"""
function run_study(make_problem; methods = ALL_METHODS, precisions = PRECISIONS)
    runs = Run[]
    for T in precisions
        prob = make_problem(T)
        for spec in methods
            sol = nothing
            err = nothing
            try
                sol = integrate(prob, spec.method)
            catch e
                err = sprint(showerror, e)
                @warn "integration failed" method = spec.name precision = T
            end
            push!(runs, Run(spec, T, prob, sol, err))
        end
    end
    return runs
end
