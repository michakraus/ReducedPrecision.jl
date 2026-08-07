# A tableau-driven initial guess for the implicit partitioned Runge-Kutta integrators.
#
# The upstream `initial_guess!` for `IPRK` on a partitioned problem builds an absolute stage time
# and hands it to `solutionstep!`, which differences it back down to a normalised one:
#
#     soltmp = (t = history[1].t + timestep(int) * tableau(int).q.c[i], …)
#     solutionstep!(soltmp, history, problem(int), iguess(int))   # HermiteExtrapolation
#
# and inside `HermiteExtrapolation` that becomes `Δt = t₁ - t₀`, `s = (tᵢ - t₀) / Δt`. Every one of
# those quantities is stored in the working precision, so in half precision the round trip
# `c → t → c` loses accuracy — and once the clock saturates it fails outright with `t₀ == t₁`.
#
# But the normalised time is a *tableau constant*: stage `i` sits at `c[i]` steps past the previous
# step by construction. `NormalizedHermiteExtrapolation` takes that constant directly and never looks
# at a clock, so passing `tableau(int).q.c[i]` to `extrapolate!` removes the round trip entirely —
# the initial guess then involves no time arithmetic at all, at any precision.
#
# ## Scope
#
# This replaces the upstream method for **every** `IPRK` method on a partitioned problem, not just the
# two Gauss rules that motivated it. That is forced rather than chosen: `initmethod` collapses every
# `IPRKMethod` (and every `IRKMethod`, on a partitioned problem) into the same `IPRK{PartitionedTableau{T}}`
# integrator type, so the originating method is not recoverable at dispatch time. It is also the
# coherent choice — `Gauss(2)` and `PRK Gauss(2)` reduce to the identical integrator here, so giving
# them different initial guesses would produce two differently-behaving copies of one method. The
# substitution is an improvement for all of them: same interpolant, same nodes, less round-off.
#
# The signature narrows the problem type to `Union{PODEProblem, HODEProblem}` — the two forms the
# study actually uses — rather than the upstream `AbstractProblemPODE` (which also covers the DAE
# forms). That keeps this an *added* method rather than one overwriting the upstream definition, which
# Julia rejects during precompilation.
#
# This file needs `NormalizedHermiteExtrapolation`, i.e. GeometricIntegratorsBase ≥ v0.4.2. The
# package's `[compat]` bound is tighter than that (v0.5.1) for the reasons given in `study.jl`.

"""
    initial_guess!(sol, history, params, int::GeometricIntegrator{<:IPRK,<:Union{PODEProblem,HODEProblem}})

Initial guess for the internal stages of an implicit partitioned Runge–Kutta step, extrapolated from
the two previous steps with [`NormalizedHermiteExtrapolation`](@ref).

Replaces the upstream method, which reaches the same interpolant through an absolute stage time that
is then differenced back into a normalised one. Here the normalised time is read straight off the
tableau (`c[i]`), so no clock value enters the guess — which is what keeps it accurate in half
precision, where the differenced version degrades and eventually fails with `t₀ == t₁`.

`NormalizedHermiteExtrapolation` expects derivative samples scaled by `Δt` and returns a `Δt`-scaled
derivative, whereas the `IPRK` cache holds unscaled vector-field values; the scaling is applied on the
way in and undone on the way out. The two formulations agree exactly in exact arithmetic, and the
extrapolated *state* agrees to round-off in every precision. The extrapolated *derivative* rounds
differently: the Hermite derivative basis differences two terms of size `6c(1+c)·x/Δt`, which for
`Δt = 0.1` is a cancellation of quantities ~150× the result, so both formulations lose accuracy there
in half precision (the normalised one cancels at `1/Δt` smaller magnitude before rescaling). This is
harmless — an initial guess can only affect how many solver iterations a step needs, never what the
step converges to — and empirically every method converges at every precision.
"""
function GeometricIntegratorsBase.initial_guess!(sol, history, params,
        int::GeometricIntegrator{<:IPRK,<:Union{PODEProblem,HODEProblem}})
    local x = nlsolution(int)
    local C = cache(int)
    local Δt = timestep(int)

    # The two Hermite samples: history[2] is the earlier one (normalised time -1), history[1] the
    # previous step (normalised time 0). The extrapolation is normalised, so its derivative arguments
    # are Δt·v rather than v.
    q₀, q₁ = history[2].q, history[1].q
    p₀, p₁ = history[2].p, history[1].p
    Δv₀, Δv₁ = Δt .* history[2].q̇, Δt .* history[1].q̇
    Δf₀, Δf₁ = Δt .* history[2].ṗ, Δt .* history[1].ṗ

    for i in eachstage(int)
        # Stage i is c[i] steps past history[1] by definition of the tableau — no time arithmetic.
        extrapolate!(q₀, Δv₀, q₁, Δv₁, tableau(int).q.c[i], C.Q[i], C.V[i],
            NormalizedHermiteExtrapolation())
        extrapolate!(p₀, Δf₀, p₁, Δf₁, tableau(int).p.c[i], C.P[i], C.F[i],
            NormalizedHermiteExtrapolation())
        # undo the normalisation: the cache holds vector-field values, not Δt-scaled ones
        C.V[i] ./= Δt
        C.F[i] ./= Δt
    end

    # assemble the initial guess for the nonlinear solver's solution vector (as upstream)
    D = length(C.V[1])
    for i in eachstage(int)
        for k in 1:D
            x[2*(D*(i-1)+k-1)+1] = 0
            x[2*(D*(i-1)+k-1)+2] = 0
            for j in eachstage(int)
                x[2*(D*(i-1)+k-1)+1] += tableau(int).q.a[i, j] * C.V[j][k]
                x[2*(D*(i-1)+k-1)+2] += tableau(int).p.a[i, j] * C.F[j][k]
            end
        end
    end

    return nothing
end
