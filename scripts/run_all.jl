# Run every reduced-precision experiment script, regenerating all figures in plots/.
#
#     julia --project=. scripts/run_all.jl   (or: julia scripts/run_all.jl)
#
# Every `scripts/*.jl` is discovered automatically (so new examples are picked up without editing
# this script) and run in a *single* Julia session. The shared packages (ReducedPrecision,
# CairoMakie, GeometricIntegrators, …) are then compiled only once, instead of once per script as
# when a fresh `julia` is spawned per file — much faster over the full suite. Each script runs the
# full method × precision sweep for one problem/scenario, verifies precision purity, and writes its
# energy-error, solution-error and trajectory figures to plots/ (filenames encode the timestep,
# e.g. `…_dt_0.1_euler.png`). This regenerates the figures independently of the documentation build.
#
# Each script is `include`d into its own throwaway module so their top-level names (`Δt`,
# `make_problem`, `plotdir`, …) do not clash across scripts. A script that throws is reported and
# makes the run exit non-zero, but does not stop the remaining scripts (per-run integration
# failures are already caught inside `run_study` and do not fail a script).

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# every scripts/*.jl except this runner, in a stable (alphabetical) order
scripts = sort(filter(readdir(@__DIR__; join = true)) do path
    endswith(path, ".jl") && basename(path) != basename(@__FILE__)
end)

failed = String[]

for (i, path) in enumerate(scripts)
    name = basename(path)
    println("==================== $name ====================")
    try
        # fresh module per script isolates its top-level constants/definitions
        Base.include(Module(Symbol("Script_", i)), path)
    catch e
        push!(failed, name)
        @error "$name FAILED" exception = (e, catch_backtrace())
    end
end

println()
if isempty(failed)
    println("All $(length(scripts)) experiment scripts completed; figures written to plots/.")
else
    println("Completed with failures: ", join(failed, " "))
    exit(1)
end
