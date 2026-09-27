using ReducedPrecision
using Test

using GeometricProblems.HarmonicOscillator: hamiltonian

include("helpers/problems.jl")

@testset "tableau-driven initial guess" begin
    # The Gauss rules must still reproduce the methods they are equivalent to. Gauss(1) *is* the
    # implicit midpoint rule, and on a partitioned problem Gauss(2) reduces to the duplicated
    # partitioned Gauss(2) tableau — so those two must agree exactly, at every precision.
    for T in PRECISIONS
        runs = run_study(make_ho; precisions = (T,))
        a, b = runof(runs, "Implicit Runge-Kutta 4"), runof(runs, "PRK Gauss(2)")
        @test Array(a.sol.q) == Array(b.sol.q)
        @test Array(a.sol.p) == Array(b.sol.p)
    end

    # Both Gauss rules are symplectic and should sit at the round-off floor on the oscillator.
    runs = run_study(make_ho; methods = OTHER_METHODS, precisions = (Float64,))
    for name in ("Implicit Midpoint", "Implicit Runge-Kutta 4")
        @test energy_error(runof(runs, name).sol, hamiltonian)[end] < 1e-13
    end
end
