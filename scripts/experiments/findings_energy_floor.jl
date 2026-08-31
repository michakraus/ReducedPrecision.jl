# Reproduce the energy-error floors quoted in `docs/src/findings.md`, "The role of precision".
#
#     julia --project=. scripts/experiments/findings_energy_floor.jl
#
# Two claims are measured, both on the short scenario (Δt = 0.1, t ≤ 1000) that
# `scripts/harmonic_oscillator.jl` and `scripts/pendulum.jl` run:
#
#   1. the implicit midpoint rule's bounded energy error, per precision, on the harmonic
#      oscillator and on the pendulum — the table in that section;
#   2. the fourth-order pair quoted immediately after it, to show that the pendulum's equal
#      Float32 and Float64 midpoint figures are a *truncation* floor rather than a round-off one.
#      `docs/src/findings.md` names `SPRK Gauss(2)` there and `docs/src/pendulum.md` names implicit
#      RK4 (`Gauss(2)`), so both are measured — they are different integrators on a nonlinear
#      problem, symplectic-by-construction against by-duplication.
#
# The statistic is the mean of the relative energy error over the second half of the run, which is
# what "bounded error" means quantitatively: the error oscillates rather than growing, so a mean
# taken past the transient is the level it oscillates about.
#
# This script exists because a dependency bump can move these numbers. RungeKutta 0.6 takes the
# Gauss nodes and weights from QuadratureRules rather than computing them itself, and SimpleSolvers
# 0.13 makes `LapackLU` the default linear solver — both change results in the last floating-point
# digits, which is enough to move a figure sitting on a round-off floor. Run it after any such bump
# and compare against `docs/src/findings.md` before trusting the prose.

using ReducedPrecision
import GeometricProblems.HarmonicOscillator as HO
import GeometricProblems.Pendulum as PD

const Δt = 0.1
const nt = 10_000

# The two problems, each built at precision T with initial conditions from its own module defaults
# (GeometricProblems has no `podeproblem(::Type{T})` precision constructor).
const PROBLEMS = [
    ("harmonic oscillator",
        T -> HO.podeproblem(T.(HO.q₀), T.(HO.p₀);
            timespan = (T(0.0), T(nt * Δt)), timestep = T(Δt)),
        HO.hamiltonian),
    ("pendulum",
        T -> PD.podeproblem(T.(PD.q₀), T.(PD.p₀);
            timespan = (T(0.0), T(nt * Δt)), timestep = T(Δt)),
        PD.hamiltonian)
]

const REPORTED = ["Implicit Midpoint", "Implicit Runge-Kutta 4", "SPRK Gauss(2)"]
const SPECS = filter(m -> m.name in REPORTED, ALL_METHODS)

# Statistic behind every number quoted in the section: the mean over the second half of the run.
function tail_mean(run, hamiltonian)
    run.sol === nothing && return NaN
    e = energy_error(run.sol, hamiltonian)
    half = length(e) ÷ 2
    return sum(@view e[(half + 1):end]) / (length(e) - half)
end

# One sweep per problem, both methods and all four precisions in it.
results = Dict{Tuple{String, String, DataType}, Float64}()
for (label, make_problem, ham) in PROBLEMS
    println("==================== $label ====================")
    runs = run_study(make_problem; methods = SPECS)
    verify_precision(runs)
    for run in runs
        results[(label, run.method.name, run.precision)] = tail_mean(run, ham)
    end
    flush(stdout)
end

fmt(x) = isnan(x) ? "failed" : string(round(x; sigdigits = 2))

for name in REPORTED
    println("\n$name — mean relative energy error over the second half " *
            "(Δt = $Δt, t ≤ $(nt * Δt))\n")
    println("| | ", join(nameof.(PRECISIONS), " | "), " |")
    println("|:--|", repeat(":--|", length(PRECISIONS)))
    for (label, _, _) in PROBLEMS
        row = [fmt(results[(label, name, T)]) for T in PRECISIONS]
        println("| ", label, " | ", join(row, " | "), " |")
    end
end
flush(stdout)
