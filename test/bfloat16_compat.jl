using ReducedPrecision
using Test
import NaNMath

@testset "BFloat16 compatibility shims" begin
    # `src/bfloat16_compat.jl` fills the gaps BFloat16s.jl and NaNMath leave. The NaNMath ones
    # are load-bearing for the whole BFloat16 column: GeometricProblems passes `nanmath = true`
    # to every symbolic generation, so an EulerLagrange vector field reaches the NaNMath variant
    # of *every* elementary function it contains — `cos` for the double pendulum, `log` for the
    # Lotka–Volterra one-form ϑ — and a missing one is a `MethodError` per run.
    guarded = (
        :sin, :cos, :tan, :asin, :acos, :atanh, :log, :log2, :log10, :log1p, :acosh)
    for f in guarded
        g = getfield(NaNMath, f)
        x = f === :acosh ? BFloat16(2.0) : BFloat16(0.5)      # in-domain for all of them
        @test g(x) isa BFloat16                                # method exists, stays in type
        @test g(x) === BFloat16(g(Float32(x)))                 # agrees with the Float32 value
    end
    # the domain guards return NaN rather than throwing, which is the whole point
    @test isnan(NaNMath.log(BFloat16(-1)))
    @test isnan(NaNMath.asin(BFloat16(2)))
    @test isnan(NaNMath.acosh(BFloat16(0.5)))
    @test isnan(NaNMath.cos(BFloat16(Inf)))
    # NaNMath's own `sqrt`/`pow`/`max`/`min` are generic enough to need no shim
    @test NaNMath.sqrt(BFloat16(0.25)) === BFloat16(0.5)
    @test isnan(NaNMath.sqrt(BFloat16(-1)))
    @test isnan(NaNMath.pow(BFloat16(-2), BFloat16(0.5)))

    # the Base gaps: `rem` (hence the float range behind `Solution`), `Integer`, `sincos`
    @test rem(BFloat16(5), BFloat16(3)) === BFloat16(2)
    @test Integer(BFloat16(7)) == 7
    @test sincos(BFloat16(0.5)) === (sin(BFloat16(0.5)), cos(BFloat16(0.5)))
    @test BFloat16(big(3)) === BFloat16(3.0)
end
