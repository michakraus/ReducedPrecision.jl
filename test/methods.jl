using ReducedPrecision
using Test

@testset "method registry" begin
    @test length(ALL_METHODS) == 12
    @test length(GEOMETRIC_METHODS) == 4
    @test length(NONGEOMETRIC_METHODS) == 4
    @test length(GAUSS2_METHODS) == 4
    @test all(m.geometric for m in GEOMETRIC_METHODS)
    @test all(!m.geometric for m in NONGEOMETRIC_METHODS)
    @test all(m.geometric for m in GAUSS2_METHODS)          # partitioned Gauss(2): all symplectic

    # the three plotting groups partition all methods (no overlap, nothing missing)
    @test length(EULER_METHODS) == 4
    @test length(OTHER_METHODS) == 4
    groupnames = [Set(m.name for m in EULER_METHODS),
        Set(m.name for m in OTHER_METHODS),
        Set(m.name for m in GAUSS2_METHODS)]
    @test union(groupnames...) == Set(m.name for m in ALL_METHODS)
    @test sum(length, groupnames) == length(ALL_METHODS)      # pairwise disjoint

    # the "other" group is a 2x2: explicit/implicit at order 2, then at order 4
    @test [m.name for m in OTHER_METHODS] ==
          ["Explicit Midpoint", "Explicit Runge-Kutta 4",
        "Implicit Midpoint", "Implicit Runge-Kutta 4"]
    @test [m.geometric for m in OTHER_METHODS] == [false, false, true, true]
    @test [g.first for g in METHOD_GROUPS] == ["euler", "other", "gauss2"]
end
