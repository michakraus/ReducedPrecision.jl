using ReducedPrecision
using Test

using GeometricBase: datatype, timetype, ntime
using GeometricIntegrators: Gauss
using GeometricIntegratorsBase: ExplicitEuler, GeometricIntegrator, default_options,
                                initmethod
import GeometricIntegratorsBase
import SimpleSolvers
using GeometricProblems.HarmonicOscillator: podeproblem, odeproblem
import GeometricProblems.HarmonicOscillator as HO

include("helpers/problems.jl")

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

@testset "capped_final_time (clock-saturation diagnostic)" begin
    # Where a T-typed *global* clock stops advancing. Float32/Float64 resolve these horizons
    # outright; the half precisions saturate, BFloat16 far earlier than Float16 (8 vs 11
    # significand bits). These are the numbers the local time frame and the tableau-driven initial
    # guess exist to defeat.
    @test capped_final_time(Float64, 1000.0, 0.1) == 1000.0
    @test capped_final_time(Float32, 1000.0, 0.1) == 1000.0
    @test capped_final_time(Float16, 1000.0, 0.1) == 128.0
    @test capped_final_time(BFloat16, 1000.0, 0.1) == 16.25
    @test capped_final_time(BFloat16, 10.0, 0.01) ≈ 2.0 atol = 0.1
    # a coarser step stays resolvable for longer
    @test capped_final_time(BFloat16, 10_000.0, 1.0) == 256.0
    @test capped_final_time(Float16, 10_000.0, 1.0) == 2048.0
end

@testset "solver tolerance is precision- and size-scaled upstream" begin
    # `run_study` passes no tolerances of its own, so the sweep depends on
    # `GeometricIntegratorsBase.default_options` scaling `f_abstol` as
    # `max(8, solversize(method, problem)) * eps(datatype(problem))`. These assertions pin that
    # property: they are the tripwire that fires if a release reverts to an absolute floor fixed
    # at `8eps(Float64)`, which a half-precision residual can never reach.
    fabstol(prob, method) = Float64(default_options(initmethod(method, prob), prob).f_abstol)

    # The precision scaling: an implicit solve is never asked for a residual its own arithmetic
    # cannot express.
    for T in PRECISIONS
        @test fabstol(make_ho(T), Gauss(2)) == 8 * Float64(eps(T))
    end
    for (lo, hi) in ((Float64, Float32), (Float32, Float16), (Float16, BFloat16))
        @test fabstol(make_ho(hi), Gauss(2)) > fabstol(make_ho(lo), Gauss(2))
    end

    # The size scaling, which is what keeps the high-stage-count `Gauss(8)` references from
    # stalling. It only bites above the `max(8, …)` floor, so it needs a stage system of more
    # than 8 unknowns: the 2-dof oscillator in ODE form under `Gauss(8)` has 16. (The partitioned
    # form the sweep itself uses reports `solversize = 0` and so sits on the floor at every stage
    # count.)
    ho_ode(::Type{T}) where {T} = odeproblem(
        T.([HO.q₀[1], HO.p₀[1]]); timespan = (T(0.0), T(1.0)), timestep = T(0.1))
    @test fabstol(ho_ode(Float64), Gauss(2)) == 8eps(Float64)     # 4 unknowns: on the floor
    @test fabstol(ho_ode(Float64), Gauss(8)) == 16eps(Float64)    # 16 unknowns: above it
    @test fabstol(ho_ode(Float16), Gauss(8)) > fabstol(ho_ode(Float64), Gauss(8))

    # A `solveropts` override must reach the nonlinear solve, and must be *merged* into the
    # method's `default_options` rather than replacing them — `min_iterations = 1` decides
    # whether a step is taken at all, so losing it is silent and severe. Built the way
    # `integrate_bounded` builds it.
    int = GeometricIntegrator(make_ho(Float64), Gauss(2);
        solver = DogLeg(), f_abstol = 1e-2, verbosity = 0)
    opts = SimpleSolvers.config(GeometricIntegratorsBase.solver(int))
    @test opts.f_abstol == 1e-2          # the override arrived
    @test opts.verbosity == 0            # …including the one that silences the solve
    @test opts.min_iterations == 1       # …without dropping the framework default

    # and the same options survive a whole sweep
    runs = run_study(make_ho; methods = GAUSS2_METHODS, precisions = (Float64,),
        solveropts = (verbosity = 0,))
    @test all(r.sol !== nothing && r.error === nothing for r in runs)
end

@testset "half precision carries a saturating horizon" begin
    # t₁ = 100 at Δt = 0.1 is far past where a BFloat16 global clock stops advancing (t ≈ 16).
    # Every method must nevertheless run the full horizon and stay type-pure.
    make_long(::Type{T}) where {T} = podeproblem(
        T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(100.0)), timestep = T(0.1))

    runs = run_study(make_long; precisions = (BFloat16,))
    @test all(r.sol !== nothing for r in runs)
    @test all(r.error === nothing for r in runs)
    @test all(assert_precision(r.prob, r.sol, BFloat16) for r in runs)
    @test all(ntime(r.sol) > capped_final_time(BFloat16, 100.0, 0.1) / 0.1
    for r in runs)

    # The partitioned RK methods take their initial guess from the tableau (see
    # `initial_guess.jl`), so no clock value enters it: their results must be *identical* whether
    # the step clock is re-anchored each step or walked along the problem's saturating grid.
    old = run_study(make_long; precisions = (BFloat16,), localclock = false)
    for name in ("Implicit Midpoint", "Implicit Runge-Kutta 4", "Implicit Euler", "PRK Gauss(2)")
        a, b = runof(runs, name), runof(old, name)
        @test b.sol !== nothing
        @test Array(a.sol.q) == Array(b.sol.q)
        @test Array(a.sol.p) == Array(b.sol.p)
    end

    # The local frame is a relabelling, not a change of dynamics: it must not perturb a
    # precision whose global clock never saturated in the first place.
    for T in (Float32, Float64)
        a = run_study(make_long; methods = EULER_METHODS, precisions = (T,))
        b = run_study(make_long; methods = EULER_METHODS, precisions = (T,), localclock = false)
        for (ra, rb) in zip(a, b)
            @test maximum(abs, Float64.(Array(ra.sol.q)) .- Float64.(Array(rb.sol.q))) <
                  100 * eps(T)
        end
    end
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
    make_coarse(::Type{T}) where {T} = podeproblem(
        T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(100.0)), timestep = T(1.0))
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
