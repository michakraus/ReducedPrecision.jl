# HISTORICAL. Asked whether the upstream GeometricIntegratorsBase 0.4.0 `reset!` fix (set the step
# time from timesteps(sol)[n] instead of accumulating `t += Δt`) resolved the Float16
# time-accumulation issue that made the long-time HO/pendulum implicit runs throw
# "t₀ and t₁ … identical". It did not — the answer was a *representability* limit of the T-typed
# clock, which is why `capped_final_time` existed.
#
# Both the fix and the cap have since been superseded: `integrate_bounded` advances the step in a
# local time frame (`reset_local!`) and the partitioned RK initial guess is taken from the tableau's
# normalised nodes (`src/initial_guess.jl`), so the clock no longer constrains the horizon at all and
# no script caps its final time. Kept because the Float16 grid-collision census below still documents
# the underlying limit; the "uncapped vs capped" comparison it runs is now expected to succeed in
# both columns.
#
#     julia --project=. scripts/experiments/verify_reset_fix.jl

using ReducedPrecision
using Printf
using GeometricIntegrators
using GeometricIntegratorsBase: Solution, solutionstep, current
using GeometricBase: timesteps, timespan
import GeometricIntegratorsBase as GIB

using GeometricProblems.HarmonicOscillator: podeproblem as ho_pode, hamiltonian as ho_ham
import GeometricProblems.HarmonicOscillator as HO

println("GeometricIntegratorsBase version in use: ",
    pkgversion(GIB))

# --- 1. the stored Float16 time grid: does it still collide past ~2048? ------------------------
println("\n[1] Float16 stored time grid (Δt = 1), collisions = consecutive equal timestamps:")
for (t₁, label) in ((2000.0, "capped t≤2000"), (10_000.0, "full t≤10000"))
    prob = ho_pode(Float16.(HO.q₀), Float16.(HO.p₀);
        timespan = (Float16(0), Float16(t₁)), timestep = Float16(1))
    sol = Solution(prob)
    t = timesteps(sol)                       # OffsetVector of Float16 times
    tv = Float64.(collect(t))
    ncol = count(i -> tv[i] == tv[i-1], 2:length(tv))
    firstcol = findfirst(i -> tv[i] == tv[i-1], 2:length(tv))
    @printf("    %-14s  n=%5d  distinct=%5d  collisions=%5d  first collision at t≈%s\n",
        label, length(tv), length(unique(tv)), ncol,
        firstcol === nothing ? "none" : string(tv[firstcol+1]))
end

# --- 2. actual integrations: full (uncapped) vs capped Float16 horizon -------------------------
# Default Hermite initial guess (the one that previously threw t₀==t₁). Implicit methods only.
methods = [("Implicit Midpoint",      Gauss(1)),
           ("Implicit Euler",         ImplicitEulerRK()),
           ("Implicit Runge-Kutta 4", Gauss(2))]

function try_run(t₁, mname, method)
    prob = ho_pode(Float16.(HO.q₀), Float16.(HO.p₀);
        timespan = (Float16(0), Float16(t₁)), timestep = Float16(1))
    try
        sol, diverged = integrate_bounded(prob, method)
        ee = ReducedPrecision.energy_error(sol, ho_ham)
        last = diverged === nothing ? length(ee) : diverged
        fe = filter(isfinite, ee[1:min(last, length(ee))])
        status = diverged === nothing ? @sprintf("done %5d steps", last-1) :
                 @sprintf("diverged @ %d", diverged)
        @sprintf("%-18s  E_final=%.2e", status, isempty(fe) ? NaN : fe[end])
    catch e
        @sprintf("THREW: %.48s", first(split(sprint(showerror, e), '\n')))
    end
end

for (t₁, label) in ((2000.0, "Float16 capped t≤2000"), (10_000.0, "Float16 full  t≤10000"))
    println("\n[2] $label — implicit methods, default Hermite guess:")
    for (mname, method) in methods
        @printf("    %-17s : %s\n", mname, try_run(t₁, mname, method))
    end
end
