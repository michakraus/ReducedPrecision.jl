# Gaps in BFloat16s.jl (v0.6.1) and NaNMath that the integration stack needs.
#
# BFloat16s defines arithmetic and the elementary functions by round-tripping through `Float32`
# (see its `bfloat16.jl`), but not a handful of two-argument ones. `Base` has no generic
# `AbstractFloat` fallback for those — it dispatches only on `IEEEFloat` (`Float16/32/64`) and errors
# out with `no_op_err` otherwise — so they are filled in here, following BFloat16s' own pattern.
# `Float32` covers `BFloat16` exactly, so each result is correctly rounded.
#
# The one we cannot do without is `rem`: `Base.:(:)` computes the length of a float range as
# `fld(stop - start, step) + 1`, and `fld`/`div`/`mod`/`cld`/`divrem` all route through `rem`. That
# range is how `GeometricSolutions.TimeSeries(tbegin, tend, Δt)` — and hence `Solution(problem)` —
# builds the time grid, so without it no `BFloat16` problem can even be allocated.
#
# Everything here extends a `Base` or `NaNMath` function on a type owned by BFloat16s, i.e. formally
# piracy. All of it is additive (no existing method is shadowed) and all of it belongs upstream —
# the `Base` methods in BFloat16s, the `NaNMath` ones in NaNMath alongside its `Float16` branch.
# Remove each once its owner ships it.

Base.rem(x::BFloat16, y::BFloat16) = BFloat16(rem(Float32(x), Float32(y)))
Base.atan(x::BFloat16, y::BFloat16) = BFloat16(atan(Float32(x), Float32(y)))
Base.fma(x::BFloat16, y::BFloat16, z::BFloat16) = BFloat16(fma(Float32(x), Float32(y), Float32(z)))

# `Base.sincos` recurses into itself for BFloat16 (it has no `IEEEFloat` guard), overflowing the
# stack rather than erroring; route it through the scalar `sin`/`cos` that BFloat16s does define.
Base.sincos(x::BFloat16) = (sin(x), cos(x))

# `Integer(x)` is likewise `IEEEFloat`-only (`boot.jl`: `Integer(x::Union{Float16,Float32,Float64})`).
# The float-range length is `convert(Integer, fld(stop - start, step)) + 1`, so this is the second
# half of what `Solution(problem)` needs.
Base.Integer(x::BFloat16) = Integer(Float32(x))

# BFloat16s constructs from the fixed-width integer types but not from `BigInt`, which is how the
# tableau coefficients of the fixed-coefficient rules (`ExplicitMidpoint`, `RK4`) arrive — `RungeKutta`
# builds them from exact rationals. This also covers `Rational{BigInt}`, since
# `convert(::Type{<:AbstractFloat}, ::Rational)` divides the converted numerator by the converted
# denominator.
BFloat16(x::Integer) = BFloat16(Float32(x))

# `NaNMath` returns `NaN` where `Base` would throw a `DomainError`, which is what lets an implicit
# solver reject a trial state outside a model's domain instead of aborting the run. Its
# domain-guarded variants are defined only for `Union{Float16,Float32,Float64}`; every other `Real`
# falls through to a generic method that throws a `MethodError` when `float(x) === x`, which
# `BFloat16` satisfies.
#
# The whole guarded set is needed, because GeometricProblems passes `nanmath = true` to every
# symbolic generation: an EulerLagrange-generated vector field uses the NaNMath variant of every
# elementary function it contains — `cos` for the double pendulum, `log` for the Lotka–Volterra
# one-form ϑ, and so on. NaNMath's `lgamma` is the one omission (`Base` has none to delegate to);
# its `sqrt`, `pow`, `max` and `min` need no shim, being generic over `AbstractFloat`/`Real`.
#
# The guards below mirror NaNMath's own. `Base` supplies each function for `BFloat16`, so the domain
# check is all that is missing.
for f in (:sin, :cos, :tan)                      # NaN at ±Inf
    @eval NaNMath.$f(x::BFloat16) = isinf(x) ? BFloat16(NaN) : Base.$f(x)
end
for f in (:asin, :acos, :atanh)                  # NaN outside [-1, 1]
    @eval NaNMath.$f(x::BFloat16) = abs(x) > BFloat16(1) ? BFloat16(NaN) : Base.$f(x)
end
for f in (:log, :log2, :log10)                   # NaN below 0
    @eval NaNMath.$f(x::BFloat16) = x < 0 ? BFloat16(NaN) : Base.$f(x)
end
NaNMath.acosh(x::BFloat16) = x < BFloat16(1) ? BFloat16(NaN) : Base.acosh(x)
NaNMath.log1p(x::BFloat16) = x < BFloat16(-1) ? BFloat16(NaN) : Base.log1p(x)
