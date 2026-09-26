using ReducedPrecision
using Test

using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution
import GeometricProblems.HarmonicOscillator as HO

include("helpers/problems.jl")

@testset "diagnostics" begin
    runs = run_study(make_ho; precisions = (Float64,))
    sea = runof(runs, "Symplectic Euler A")

    ee = energy_error(sea.sol, hamiltonian)
    tv = timevalues(sea.sol)
    ref = exact_solution(make_ho(Float64))
    se = solution_error(sea.sol, ref)

    # consistent lengths across the three metrics
    @test length(ee) == length(tv) == length(se)
    @test length(tv) ≥ 2

    # time grid: starts at 0, ends at t₁, strictly increasing
    @test tv[1] == 0.0
    @test tv[end] ≈ 1.0
    @test issorted(tv)

    # energy error: 0 at t₀, finite and non-negative everywhere
    @test ee[1] == 0.0
    @test all(isfinite, ee)
    @test all(≥(0), ee)

    # solution error vs the analytic reference: 0 at t₀, finite
    @test se[1] ≈ 0.0 atol = 1e-12
    @test all(isfinite, se)

    # a finer-grid reference (here Δt/2) is subsampled onto the solution's coarser grid
    make_ho_fine(::Type{T}) where {T} = podeproblem(
        T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(1.0)), timestep = T(0.05))
    se_fine = solution_error(sea.sol, exact_solution(make_ho_fine(Float64)))
    @test length(se_fine) == length(tv)          # subsampled to the solution grid
    @test se_fine[1] ≈ 0.0 atol = 1e-12
    @test all(isfinite, se_fine)

    # explicit Euler drifts more in energy than symplectic Euler by the final step
    ee_exp = energy_error(runof(runs, "Explicit Euler").sol, hamiltonian)
    @test ee_exp[end] > ee[end]
end
