#!/usr/bin/env bash
#
# Run every reduced-precision experiment script, regenerating all figures in plots/.
#
# Each script runs the full method × precision sweep for one problem/scenario, verifies
# precision purity, and writes its energy-error, solution-error and trajectory figures to
# plots/ (filenames encode the timestep, e.g. `…_dt=0.1_euler.png`). This regenerates the
# figures independently of the documentation build.
#
# Usage:  ./run_all.sh
#
# A single script that throws is reported at the end and makes the run exit non-zero, but does
# not stop the remaining scripts (per-run integration failures are already caught inside
# run_study and do not fail a script).

set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

# Make sure the environment is instantiated before the timing loop.
julia --project=. -e 'using Pkg; Pkg.instantiate()'

scripts=(
    harmonic_oscillator
    harmonic_oscillator_longtime
    pendulum
    pendulum_longtime
    double_pendulum
    double_pendulum_longtime
    toda_lattice
    toda_lattice_longtime
    lotka_volterra_2d
    lotka_volterra_2d_longtime
    lotka_volterra_4d
    lotka_volterra_4d_longtime
)

failed=()
for s in "${scripts[@]}"; do
    echo "==================== $s ===================="
    if ! julia --project=. "scripts/$s.jl"; then
        failed+=("$s")
        echo "!!! $s FAILED"
    fi
done

echo
if ((${#failed[@]})); then
    echo "Completed with failures: ${failed[*]}"
    exit 1
fi
echo "All ${#scripts[@]} experiments completed; figures written to plots/."
