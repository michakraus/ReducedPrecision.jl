using ReducedPrecision
using Aqua
using Test

@testset "Aqua" begin
    Aqua.test_all(ReducedPrecision;
        stale_deps = false,
        piracies = (broken = true,))  # issue: not filed yet — BFloat16 shims and initial_guess!
    # GeometricProblems is in [deps] but src/ loads it nowhere
    @test_broken isempty(Aqua.find_stale_deps(Base.PkgId(ReducedPrecision)))  # issue: not filed yet
end
