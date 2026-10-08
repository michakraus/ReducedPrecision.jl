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

### K5 · No test sees the scaling of the initial guess.

- **location:** `src/initial_guess.jl:78–79`
- **evidence:** the mutant `C.V[i] ./= Δt` → `C.V[i] ./= 2Δt` (or `.*= 2Δt`) survives
  `test/initial_guess.jl` and the whole suite. Both methods that the file compares take the same
  guess, and the nonlinear solver converges from any finite guess, so the results do not change.
  A test that sees the guess would check it directly, or count the solver iterations.
  `julia --startup-file=no mutate.jl <repository> src/initial_guess.jl 'C.V[i] ./= Δt' 'C.V[i] ./= 2Δt' initial_guess.jl`
  gives `SURVIVED`.
- **kind:** missing test
- **found:** #27

### K6 · The label of `test/initial_guess.jl` does not match what it tests.

- **location:** `test/initial_guess.jl:8`
- **evidence:** the testset "tableau-driven initial guess" compares the Gauss(2) method with the
  partitioned Gauss(2) tableau, and checks the energy error. It does not call `initial_guess!`
  directly.
- **kind:** docs
- **found:** #27

### K7 · Two test files import names that `test/helpers/problems.jl` already imports.

- **location:** `test/study.jl:10–11`
- **evidence:** `test/study.jl:10–11` and `test/diagnostics.jl:4–5` import `podeproblem` and
  `GeometricProblems.HarmonicOscillator as HO`, which `test/helpers/problems.jl` already imports.
  The duplicates do no harm.
- **kind:** dead code
- **found:** #27

### K8 · `BFloat16(::Integer)` overwrites the BFloat16s method on aarch64 with Julia 1.11 and 1.12, so the package does not precompile there.

- **location:** `src/bfloat16_compat.jl:39`
- **evidence:** GitHub issue #28 reports `Method overwriting is not permitted during Module
  precompilation` on Julia 1.10. The same error shows with BFloat16s 0.6.2 on aarch64 with Julia
  1.10 to 1.12, 1.11.9 among them. BFloat16s 0.5.1 to 0.6.2 define `BFloat16(x::Integer)` only
  where `llvm_arithmetic` in their `src/bfloat16.jl` is false: without `Core.BFloat16`, which is
  Julia 1.10, and on aarch64 below LLVM 19, which is Julia 1.11 (LLVM 16) and 1.12 (LLVM 18). On
  x86_64 from Julia 1.11, and on aarch64 from Julia 1.13 (LLVM 20), BFloat16s does not define it,
  and nothing is overwritten. The package still loads, and the suite passes.
  `test/quality/aqua.jl:8` keeps `persistent_tasks = false`.
- **kind:** found late
- **found:** 2026-10-08

### K9 · The **Findings** energy-error figures are not measured at the 0.3.0 floors.

- **location:** `scripts/experiments/findings_energy_floor.jl`
- **evidence:** The `[Unreleased]` figures in `CHANGELOG.md` and `docs/src/findings.md` were measured
  on GeometricIntegrators 0.18.4, GeometricIntegratorsBase 0.6.4, RungeKutta 0.6.1 and
  SimpleSolvers 0.13.2. The `[compat]` floors are now 0.18.6, 0.6.9, 0.6.4 and 0.14.1, which exclude
  all four. A run of `scripts/experiments/findings_energy_floor.jl` in an environment resolved at
  the new floors answers it.
- **kind:** not verified
- **found:** 2026-10-08
