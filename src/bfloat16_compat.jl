# Gaps in BFloat16s.jl (v0.6.1) that the integration stack needs.
#
# BFloat16s defines arithmetic and the elementary functions by round-tripping through `Float32`
# (see its `bfloat16.jl`), but a few two-argument functions were never added. `Base` has no generic
# `AbstractFloat` fallback for them — it dispatches only on `IEEEFloat` (`Float16/32/64`) and errors
# out with `no_op_err` otherwise — so they have to be filled in here. The definitions below follow
# BFloat16s' own pattern; `Float32` covers `BFloat16` exactly, so each result is correctly rounded.
#
# The one we cannot do without is `rem`: `Base.:(:)` computes the length of a float range as
# `fld(stop - start, step) + 1`, and `fld`/`div`/`mod`/`cld`/`divrem` all route through `rem`. That
# range is how `GeometricSolutions.TimeSeries(tbegin, tend, Δt)` — and hence `Solution(problem)` —
# builds the time grid, so without it no `BFloat16` problem can even be allocated.
#
# These are extensions of `Base` functions on a type owned by BFloat16s, i.e. formally piracy. They
# are additive (no existing method is shadowed) and belong upstream; remove them once BFloat16s
# ships its own.

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

# `NaNMath` returns `NaN` where `Base` would throw a `DomainError`, and the Lotka–Volterra problems
# use `NaNMath.log` in their one-form ϑ. It defines those variants only for
# `Union{Float16,Float32,Float64}`; every other `Real` falls through to a generic method that throws a
# `MethodError` when `float(x) === x`, which `BFloat16` satisfies. Mirror its `log` semantics — the
# only NaNMath function GeometricProblems calls.
NaNMath.log(x::BFloat16) = x < 0 ? BFloat16(NaN) : Base.log(x)
