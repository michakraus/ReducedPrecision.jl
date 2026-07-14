# Experiment: does MidpointExtrapolation beat HermiteExtrapolation as the implicit-solve initial
# guess, especially in Float16?
#
# HermiteExtrapolation (the default `iguess` for every implicit RK / variational method used in this
# study) seeds each nonlinear solve from a two-point cubic Hermite polynomial through the last two
# history points, dividing by Δt = t₁ − t₀. When the low-precision time grid saturates (t₀ == t₁,
# which happens in Float16 once t is large enough that ulp(t) ≫ Δt) that guess degenerates.
# MidpointExtrapolation instead advances a *single* history point with the vector field (a
# Gragg/midpoint scheme + Aitken–Neville extrapolation), so it never needs two distinct times.
#
# This script integrates the implicit methods with each extrapolation, under identical solver
# settings, and records per-step nonlinear-iteration counts plus the run outcome (completed /
# diverged / threw) and final relative energy error, across all three precisions.
#
#     julia --project=. scripts/experiments/iguess_extrapolation.jl

using ReducedPrecision
using Printf
using GeometricIntegrators
using GeometricIntegratorsBase: Solution, solutionstep, current, solverstate,
    HermiteExtrapolation, MidpointExtrapolation
import GeometricIntegratorsBase as GIB
using GeometricBase: ntime, datatype
import SimpleSolvers

using GeometricProblems.HarmonicOscillator: podeproblem as ho_pode, hamiltonian as ho_ham
import GeometricProblems.HarmonicOscillator as HO
using GeometricProblems.Pendulum: podeproblem as pd_pode, hamiltonian as pd_ham
import GeometricProblems.Pendulum as PD
using GeometricProblems.DoublePendulum: hodeproblem as dp_hode, hamiltonian as dp_ham
import GeometricProblems.DoublePendulum as DP

# number of nonlinear iterations of the most recent solve, or -1 if the method's solver state does
# not expose one (keeps the loop robust across method families)
safe_iters(int) = try
    SimpleSolvers.iteration_number(solverstate(int))
catch
    -1
end

# Mirrors ReducedPrecision.integrate_bounded, but records the nonlinear-solve iteration count after
# every step so we can compare how hard the solver works under each initial guess.
function instrumented_run(problem, method; solver, initialguess, bound = 1e3, max_iterations = 100)
    integrator = GeometricIntegrator(problem, method; solver, initialguess,
        GIB.default_options(method)..., max_iterations)
    sol = Solution(problem)
    solstep = solutionstep(integrator, sol[0])
    curstate = current(solstep)
    nt = ntime(sol)

    iters = Int[]
    diverged = nothing
    for n in 1:nt
        integrate!(solstep, integrator)
        push!(iters, safe_iters(integrator))
        copy!(sol, curstate, n)
        q = sol.q[n]; p = sol.p[n]
        if !(all(isfinite, q) && all(isfinite, p)) ||
           (maximum(abs, q) > bound || maximum(abs, p) > bound)
            diverged = n
            break
        end
    end
    return (; sol, iters, diverged, nt)
end

# one (problem, method, precision, iguess) cell
function measure(make_problem, hamiltonian, method, ::Type{T}, iguess; solver, max_iterations, cap) where {T}
    prob = make_problem(T)
    try
        r = instrumented_run(prob, method; solver, initialguess = iguess, max_iterations)
        ee = ReducedPrecision.energy_error(r.sol, hamiltonian)
        lastok = r.diverged === nothing ? length(ee) : min(r.diverged, length(ee))
        finite_ee = filter(isfinite, ee[1:lastok])
        final_energy = isempty(finite_ee) ? NaN : finite_ee[end]
        its = r.iters
        nfail = count(>=(cap), its)          # steps hitting the iteration cap → non-converged
        return (; ok = true, threw = "", diverged = r.diverged, nsteps = length(its),
                meaniter = isempty(its) ? NaN : sum(its) / length(its),
                maxiter = isempty(its) ? 0 : maximum(its), nfail, final_energy)
    catch e
        return (; ok = false, threw = first(split(sprint(showerror, e), '\n')),
                diverged = nothing, nsteps = 0, meaniter = NaN, maxiter = 0,
                nfail = 0, final_energy = NaN)
    end
end

# problem definitions matching the study scripts. HO/Pendulum use the long-time / coarse-step
# horizon (Δt = 1, t ≤ 10_000) of the `*_longtime.jl` scripts, where the Float16 time grid saturates
# hardest (t + 1 == t once t ≳ 2048). The double pendulum uses its short-step stiff scenario.
ho_make(::Type{T}) where {T} =
    ho_pode(T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(10_000.0)), timestep = T(1.0))  # nt = 10000
pd_make(::Type{T}) where {T} =
    pd_pode(T.(PD.q₀), T.(PD.p₀); timespan = (T(0.0), T(10_000.0)), timestep = T(1.0))  # nt = 10000
dp_make(::Type{T}) where {T} =
    dp_hode(T.(DP.θ₀), T.(DP.p₀); timespan = (T(0.0), T(10.0)), timestep = T(0.01),
        parameters = DP.default_parameters(T))                                          # nt = 1000

const IMPLICIT_METHODS = [
    ("Implicit Midpoint", ImplicitMidpoint()),
    ("Implicit Euler",    ImplicitEulerRK()),
    ("Crank-Nicolson",    CrankNicolson()),
    ("PRK Gauss(2)",      ReducedPrecision.GaussPRK()),
]

const CAP = 100   # iteration cap; a step reaching it counts as a non-converged solve

const PROBLEMS = [
    ("Harmonic osc. (Δt=1, t≤10000)",  ho_make, ho_ham, DogLeg()),
    ("Pendulum (Δt=1, t≤10000)",       pd_make, pd_ham, DogLeg()),
    ("Double pendulum (Δt=0.01,t≤10)", dp_make, dp_ham, Newton()),  # matches its study script
]

fresh_iguess(name) = name == "Hermite" ? HermiteExtrapolation() : MidpointExtrapolation()

function fmt(m)
    m.ok || return @sprintf("THREW: %.40s", m.threw)
    outcome = m.diverged === nothing ? @sprintf("done  %5d steps", m.nsteps) :
              @sprintf("DIVERGED @ step %-5d", m.diverged)
    @sprintf("%-22s  mean %6.2f  max %4d  ncap %4d  E=%.2e",
        outcome, m.meaniter, m.maxiter, m.nfail, m.final_energy)
end

for (plabel, make, ham, solver) in PROBLEMS
    println("\n", "="^112)
    println(plabel, "   [solver: ", nameof(typeof(solver)), ", iteration cap: ", CAP, "]")
    println("="^112)
    for (mname, method) in IMPLICIT_METHODS
        for T in (Float16, Float32, Float64)
            mh = measure(make, ham, method, T, fresh_iguess("Hermite");  solver, max_iterations = CAP, cap = CAP)
            mm = measure(make, ham, method, T, fresh_iguess("Midpoint"); solver, max_iterations = CAP, cap = CAP)
            @printf("%-17s %-8s  Hermite : %s\n", mname, string(T), fmt(mh))
            @printf("%-17s %-8s  Midpoint: %s\n", "", "", fmt(mm))
        end
    end
end

println("\nLegend: mean/max = nonlinear iterations per step; ncap = #steps hitting the iteration ",
        "cap (non-converged); E = final relative energy error at last valid step.")
