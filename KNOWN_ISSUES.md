# Known issues

What is known to be wrong in this experiment and is not fixed. An entry leaves this file when its
fix merges, and the CHANGELOG entry of the fix names its ID.

### K1 · Not every quoted figure has been re-measured since the 0.18 / 0.13 dependency bump.

- **location:** `scripts/experiments/findings_energy_floor.jl`
- **evidence:** `scripts/experiments/findings_energy_floor.jl` covers the harmonic oscillator and the pendulum. The
  energy- and solution-error levels quoted in `docs/src/double_pendulum.md`,
  `docs/src/toda_lattice.md`, the two Lotka–Volterra pages, and the solution-error levels on
  `docs/src/harmonic_oscillator.md` were not: the sweep re-runs and the figures regenerate, but those
  numbers are carried over. They are stated to one significant figure, and the movement measured on
  the two problems that *were* re-measured would leave most of them standing — but that is an
  argument, not a measurement. Extending the script to the remaining problems is the fix.
- **kind:** not verified
- **found:** 2026-08-31

### K2 · The `1.4e-12` agreement and the `2.6–40×` speed-ups in `docs/src/findings.md` are **historical A/B figures**: they compare a fixed solver tolerance against a precision-scaled one, and the unscaled variant no longer exists in the stack.

- **location:** `docs/src/findings.md`
- **evidence:** They are unaffected by a dependency bump and were not
  re-run.
- **kind:** not verified
- **found:** 2026-08-31

### K3 · A pre-bump baseline cannot be reproduced on this machine.

- **location:** —
- **evidence:** RungeKutta 0.5.23 pulls in
  GenericLinearAlgebra 0.4.0, which does not precompile on Julia 1.13, and setting up an older Julia
  to work around it was blocked here. Any A/B against a state before this bump needs a machine with
  Julia 1.11 or 1.12 available.
- **kind:** upstream
- **found:** 2026-08-31

### K4 · A sweep emits some 3800 warnings, dominated by `DogLeg trust-region radius Δ underflowed` and Hermite extrapolation's `history[1] and history[2] are identical` — both pre-existing upstream messages, present in SimpleSolvers 0.10.1 and GeometricIntegratorsBase 0.5.3 as well.

- **location:** —
- **evidence:** They are not
  suppressed, because `verbosity = 0` would also hide the genuine non-convergence reports the study
  reads; separating the two is open work.
- **kind:** upstream
- **found:** 2026-08-31
