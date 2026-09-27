using ReducedPrecision
using Aqua
using Test

@testset "Aqua" begin
    Aqua.test_all(ReducedPrecision;
        stale_deps = false,
        persistent_tasks = false,  # issue #28
        piracies = (broken = true,))  # issue #25
    # GeometricProblems is in [deps] but src/ loads it nowhere
    @test_broken isempty(Aqua.find_stale_deps(Base.PkgId(ReducedPrecision)))  # issue #26
end
