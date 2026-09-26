# Known issues

Defects that are known and not yet fixed. Each entry gives its kind and its evidence.

## KI-1 · missing test — no test sees the scaling of the initial guess

`src/initial_guess.jl:78–79`: the mutant `C.V[i] ./= Δt` → `C.V[i] ./= 2Δt` (or `.*= 2Δt`) survives
`test/initial_guess.jl` and the whole suite. Both methods that the file compares take the same
guess, and the nonlinear solver converges from any finite guess, so the results do not change.
A test that sees the guess would check it directly, or count the solver iterations.

Reproduce: `julia --startup-file=no mutate.jl <repository> src/initial_guess.jl 'C.V[i] ./= Δt' 'C.V[i] ./= 2Δt' initial_guess.jl`
gives `SURVIVED`.

## KI-2 · docs — the label of `test/initial_guess.jl` does not match what it tests

The testset "tableau-driven initial guess" compares the Gauss(2) method with the partitioned Gauss(2)
tableau, and checks the energy error. It does not call `initial_guess!` directly. The label is the
one the testset had in the single `test/runtests.jl`.

## KI-3 · dead code — duplicate imports in two test files

`test/study.jl:10–11` and `test/diagnostics.jl:4–5` import `podeproblem` and
`GeometricProblems.HarmonicOscillator as HO`, which `test/helpers/problems.jl` already imports.
The duplicates do no harm.
