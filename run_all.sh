#!/usr/bin/env bash
#
# Run every reduced-precision experiment script, regenerating all figures in plots/.
#
# Every `scripts/*.jl` is run (discovered automatically, so new examples are picked up without
# editing this script). Each runs the full method × precision sweep for one problem/scenario,
# verifies precision purity, and writes its energy-error, solution-error and trajectory figures
# to plots/ (filenames encode the timestep, e.g. `…_dt_0.1_euler.png`). This regenerates the
# figures independently of the documentation build.
#
# Usage:  bash run_all.sh   (or ./run_all.sh)
#
# A single script that throws is reported at the end and makes the run exit non-zero, but does
# not stop the remaining scripts (per-run integration failures are already caught inside
# run_study and do not fail a script).

set -uo pipefail
shopt -s nullglob

cd "$(dirname "${BASH_SOURCE[0]}")"

# Make sure the environment is instantiated before the timing loop.
julia --project=. -e 'using Pkg; Pkg.instantiate()'

scripts=(scripts/*.jl)

failed=()
for s in "${scripts[@]}"; do
    echo "==================== $s ===================="
    if ! julia --project=. "$s"; then
        failed+=("$s")
        echo "!!! $s FAILED"
    fi
done

echo
if ((${#failed[@]})); then
    echo "Completed with failures: ${failed[*]}"
    exit 1
fi
echo "All ${#scripts[@]} experiment scripts completed; figures written to plots/."
