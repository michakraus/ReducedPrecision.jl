using GeometricProblems.HarmonicOscillator: podeproblem
import GeometricProblems.HarmonicOscillator as HO

# A small, fast harmonic-oscillator problem (10 steps) at precision T. GeometricProblems has no
# `podeproblem(::Type{T})` precision constructor, so build the T-typed initial conditions from the
# module defaults.
function make_ho(::Type{T}) where {T}
    podeproblem(T.(HO.q₀), T.(HO.p₀); timespan = (T(0.0), T(1.0)), timestep = T(0.1))
end

# Fetch the run for a given method name.
runof(runs, name) = only(filter(r -> r.method.name == name, runs))
