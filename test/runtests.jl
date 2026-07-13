using ReducedPrecision
using Test

using GeometricBase: datatype, timetype
using GeometricIntegratorsBase: ExplicitEuler
using GeometricProblems.HarmonicOscillator: podeproblem, hamiltonian, exact_solution
import GeometricProblems.HarmonicOscillator as HO

# A small, fast harmonic-oscillator problem (10 steps) at precision T. GeometricProblems v0.7.0
# dropped the `podeproblem(::Type{T})` precision constructor, so build the T-typed initial
# conditions from the module defaults.
make_ho(::Type{T}) where {T} =
    podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(1.0)), timestep = T(0.1))

# Fetch the run for a given method name.
runof(runs, name) = only(filter(r -> r.method.name == name, runs))

@testset "ReducedPrecision.jl" begin

    @testset "method registry" begin
        @test length(ALL_METHODS) == 12
        @test length(GEOMETRIC_METHODS) == 3
        @test length(NONGEOMETRIC_METHODS) == 5
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

        # the requested "other" plotting order
        @test [m.name for m in OTHER_METHODS] ==
              ["RK4", "Explicit Midpoint", "Implicit Midpoint", "Crank-Nicolson"]
        @test [g.first for g in METHOD_GROUPS] == ["euler", "other", "gauss2"]
    end

    @testset "run_study + precision purity ($T)" for T in PRECISIONS
        runs = run_study(make_ho; precisions = (T,))
        @test length(runs) == length(ALL_METHODS)
        @test all(r.precision === T for r in runs)
        # the harmonic oscillator over this short horizon runs for every method/precision
        @test all(r.sol !== nothing for r in runs)
        @test all(r.error === nothing for r in runs)
        for r in runs
            @test assert_precision(r.prob, r.sol, T)   # no implicit promotion to Float64
            @test datatype(r.sol) === T
            @test timetype(r.sol) === T
        end
    end

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
        make_ho_fine(::Type{T}) where {T} =
            podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(1.0)), timestep = T(0.05))
        se_fine = solution_error(sea.sol, exact_solution(make_ho_fine(Float64)))
        @test length(se_fine) == length(tv)          # subsampled to the solution grid
        @test se_fine[1] ≈ 0.0 atol = 1e-12
        @test all(isfinite, se_fine)

        # explicit Euler drifts more in energy than symplectic Euler by the final step
        ee_exp = energy_error(runof(runs, "Explicit Euler").sol, hamiltonian)
        @test ee_exp[end] > ee[end]
    end

    @testset "failed integrations are captured" begin
        # The special ExplicitEuler is ODE-only and errors on a partitioned problem;
        # run_study must catch it per run rather than aborting the sweep.
        bad = MethodSpec("ODE-only Euler", ExplicitEuler(), false)
        runs = run_study(make_ho; methods = [bad], precisions = (Float64,))
        @test length(runs) == 1
        @test runs[1].sol === nothing
        @test runs[1].error !== nothing
    end

    @testset "divergence guard" begin
        # a coarse oscillator (Δt = 1) makes explicit Euler blow up while the symplectic
        # methods stay bounded
        make_coarse(::Type{T}) where {T} =
            podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(100.0)), timestep = T(1.0))
        runs = run_study(make_coarse; precisions = (Float64,), bound = 1e3)

        ee = runof(runs, "Explicit Euler")
        @test ee.sol !== nothing
        @test ee.diverged !== nothing          # guard tripped
        @test 0 < ee.diverged < 100            # stopped before the end
        @test all(isnan, Array(ee.sol.q)[:, end])   # tail blanked with NaN after divergence

        sea = runof(runs, "Symplectic Euler A")
        @test sea.diverged === nothing         # symplectic stays bounded
        @test all(isfinite, Array(sea.sol.q))

        # disabling the magnitude bound lets the same explicit-Euler run proceed to the end
        runs_nobound = run_study(make_coarse;
            methods = [ee.method], precisions = (Float64,), bound = nothing)
        @test runs_nobound[1].diverged === nothing
    end

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

end
