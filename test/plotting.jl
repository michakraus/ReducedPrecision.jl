using ReducedPrecision
using Test

using GeometricProblems.HarmonicOscillator: hamiltonian, exact_solution

include("helpers/problems.jl")

@testset "plotting writes both group figures" begin
    runs = run_study(make_ho)                    # all precisions
    ref = exact_solution(make_ho(Float64))
    dir = mktempdir()

    plot_energy_error(runs, hamiltonian; path = joinpath(dir, "energy.png"), title = "t")
    @test isfile(joinpath(dir, "energy_euler.png"))
    @test isfile(joinpath(dir, "energy_other.png"))
    @test isfile(joinpath(dir, "energy_gauss2.png"))

    plot_solution_error(runs, ref; path = joinpath(dir, "solerr.png"), title = "t")
    @test isfile(joinpath(dir, "solerr_euler.png"))
    @test isfile(joinpath(dir, "solerr_other.png"))
    @test isfile(joinpath(dir, "solerr_gauss2.png"))

    plot_solution(runs; reference = ref, path = joinpath(dir, "traj.png"), title = "t")
    @test isfile(joinpath(dir, "traj_euler.png"))
    @test isfile(joinpath(dir, "traj_other.png"))
    @test isfile(joinpath(dir, "traj_gauss2.png"))

    # a script-supplied custom group set (as the Lotka–Volterra examples use)
    plot_energy_error(runs, hamiltonian; path = joinpath(dir, "grp.png"), title = "t",
        groups = ["mid" => GAUSS2_METHODS])
    @test isfile(joinpath(dir, "grp_mid.png"))
end
