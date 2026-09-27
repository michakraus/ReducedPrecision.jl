# The BFloat16 shims in `src/bfloat16_compat.jl`, run against one `BFloat16s` / `NaNMath` version
# pair without the rest of the package. This is the measurement behind the `NaNMath = "1.1.4"`
# floor: 1.1.2 returns `NaN` for `NaNMath.acosh(2.0f0)`, and 1.1.2 and 1.1.3 recurse without end in
# `NaNMath.pow` for `BFloat16`. Usage, from the repository root:
#
#   julia +1.10 scripts/experiments/nanmath_floor.jl <BFloat16s version> <NaNMath version>
using Pkg
Pkg.activate(; temp = true, io = devnull)
Pkg.add(
    [PackageSpec(name = "BFloat16s", version = ARGS[1]),
        PackageSpec(name = "NaNMath", version = ARGS[2])];
    io = devnull)

using BFloat16s, NaNMath, Test
include(joinpath(@__DIR__, "..", "..", "src", "bfloat16_compat.jl"))

println("BFloat16s ", pkgversion(BFloat16s), "  NaNMath ", pkgversion(NaNMath))
@testset "BFloat16 shims" begin
    for f in (:sin, :cos, :tan, :asin, :acos, :atanh, :log, :log2, :log10, :log1p, :acosh)
        g = getfield(NaNMath, f)
        x = f === :acosh ? BFloat16(2.0) : BFloat16(0.5)
        @test g(x) isa BFloat16
        @test g(x) === BFloat16(g(Float32(x)))
    end
    @test isnan(NaNMath.log(BFloat16(-1)))
    @test isnan(NaNMath.sqrt(BFloat16(-1)))
    @test isnan(NaNMath.pow(BFloat16(-2), BFloat16(0.5)))
    @test Integer(BFloat16(7)) == 7
    @test BFloat16(big(3)) === BFloat16(3.0)
end
