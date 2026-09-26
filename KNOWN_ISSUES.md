# Known issues

Defects that are known and not yet fixed. Each entry gives its kind and its evidence.

### K1 · No test sees the scaling of the initial guess.

- **location:** `src/initial_guess.jl:78–79`
- **evidence:** the mutant `C.V[i] ./= Δt` → `C.V[i] ./= 2Δt` (or `.*= 2Δt`) survives
  `test/initial_guess.jl` and the whole suite. Both methods that the file compares take the same
  guess, and the nonlinear solver converges from any finite guess, so the results do not change.
  A test that sees the guess would check it directly, or count the solver iterations.
  `julia --startup-file=no mutate.jl <repository> src/initial_guess.jl 'C.V[i] ./= Δt' 'C.V[i] ./= 2Δt' initial_guess.jl`
  gives `SURVIVED`.
- **kind:** missing test
- **found:** 2026-09-27

### K2 · The label of `test/initial_guess.jl` does not match what it tests.

- **location:** `test/initial_guess.jl:8`
- **evidence:** the testset "tableau-driven initial guess" compares the Gauss(2) method with the
  partitioned Gauss(2) tableau, and checks the energy error. It does not call `initial_guess!`
  directly.
- **kind:** docs
- **found:** 2026-09-27

### K3 · Two test files import names that `test/helpers/problems.jl` already imports.

- **location:** `test/study.jl:10–11`
- **evidence:** `test/study.jl:10–11` and `test/diagnostics.jl:4–5` import `podeproblem` and
  `GeometricProblems.HarmonicOscillator as HO`, which `test/helpers/problems.jl` already imports.
  The duplicates do no harm.
- **kind:** dead code
- **found:** 2026-09-27
