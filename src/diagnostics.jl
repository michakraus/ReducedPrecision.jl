# Precision-purity checks and error metrics.

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
        if run.diverged === nothing
            println("  [ ok ] $(run.method.name) @ $(run.precision)")
        else
            println("  [ ok ] $(run.method.name) @ $(run.precision) (diverged at step $(run.diverged); stopped)")
        end
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

Euclidean norm of the state error `‖(qₙ,pₙ) − (qₙ,pₙ)_ref‖` at each time step. Both solutions
are materialized to `Float64` before comparison, so `reference` may be at a different precision
than `sol` (e.g. a Float64 high-order reference). If `reference` was integrated on a **finer**
time grid that refines `sol`'s grid (an integer number of reference steps per solution step,
sharing `t₀`), it is subsampled onto `sol`'s grid before the index-by-index comparison — this is
how the coarse-step scenarios compare against a reference computed at the fine step.
"""
function solution_error(sol, reference)
    q  = Float64.(Array(sol.q));       p  = Float64.(Array(sol.p))
    qr = Float64.(Array(reference.q)); pr = Float64.(Array(reference.p))
    nsol = size(q, 2); nref = size(qr, 2)
    if nref > nsol
        # reference on a finer grid: pick the nsol evenly-spaced points (endpoints included) that
        # line up with the solution's output times. Robust to the ±1-step ambiguity from a
        # non-exact fine Δt (e.g. 0.1 in Float64).
        idx = round.(Int, range(1, nref; length = nsol))
        qr = qr[:, idx]; pr = pr[:, idx]
    end
    @assert size(q, 2) == size(qr, 2) "solution and reference have different lengths"
    n = size(q, 2)
    return Float64[
        sqrt(sum(abs2, @view(q[:, i]) .- @view(qr[:, i])) +
             sum(abs2, @view(p[:, i]) .- @view(pr[:, i]))) for i in 1:n
    ]
end
