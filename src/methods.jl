# Precisions and the integration-method registry.

"""
The floating-point precisions swept over in every study.
"""
const PRECISIONS = (Float16, Float32, Float64)

"""
    MethodSpec(name, method, geometric)

A named integration method together with a flag marking it as geometric (symplectic) or not.
The flag drives the line style in plots (geometric = solid, non-geometric = dashed).
"""
struct MethodSpec
    name::String
    method::GeometricMethod
    geometric::Bool
end

"""
Geometric (symplectic) methods.
"""
const GEOMETRIC_METHODS = MethodSpec[
    MethodSpec("Symplectic Euler A", SymplecticEulerA(), true),
    MethodSpec("Symplectic Euler B", SymplecticEulerB(), true),
    MethodSpec("Implicit Midpoint",  ImplicitMidpoint(),  true),
]

"""
Non-geometric methods. Explicit/implicit Euler are represented by their Runge–Kutta tableau
twins (`ExplicitEulerRK` / `ImplicitEulerRK`), which auto-promote to partitioned RK on a
PODE/HODE, so the whole set runs on a single partitioned problem form.
"""
const NONGEOMETRIC_METHODS = MethodSpec[
    MethodSpec("Explicit Euler",    ExplicitEulerRK(),  false),
    MethodSpec("Implicit Euler",    ImplicitEulerRK(),  false),
    MethodSpec("Explicit Midpoint", ExplicitMidpoint(), false),
    MethodSpec("Crank-Nicolson",    CrankNicolson(),    false),
    MethodSpec("RK4",               RK4(),              false),
]

# Partitioned Gauss(2) midpoint-rule variants (the third comparison group). All four are
# implicit partitioned Runge–Kutta (IPRK) methods built from the 2-stage Gauss tableau. They
# are all symplectic and differ only in implementation details:
#   * `PartitionedTableau(Gauss(2))`            — the partitioned tableau formed by duplicating
#                                                  the Gauss tableau for q and p;
#   * `SymplecticPartitionedTableau(Gauss(2))`  — the p-tableau is the symplectic conjugate of
#                                                  the q-tableau, so the symplecticity condition
#                                                  holds to floating-point accuracy by construction;
#   * the `…0` variants additionally zero the rounding-error compensation coefficients â, b̂, ĉ.
# Each defines only a `tableau(method, T)` accessor, so `initmethod` rebuilds the tableau at the
# run precision T (keeping every precision pure); everything downstream uses the concrete IPRK.
struct GaussPRK   <: IPRKMethod end   # PartitionedTableau(Gauss(2))
struct GaussSPRK  <: IPRKMethod end   # SymplecticPartitionedTableau(Gauss(2))
struct GaussPRK0  <: IPRKMethod end   # GaussPRK,  with â = b̂ = ĉ = 0
struct GaussSPRK0 <: IPRKMethod end   # GaussSPRK, with â = b̂ = ĉ = 0

# Reconstruct a tableau from its (already precision-T) a/b/c coefficients. The inner
# constructor derives the compensation terms as `â = a .- convert(T, a)` etc.; since a/b/c
# are already type T, the recomputed â/b̂/ĉ come out identically zero.
_zero_hats(t::Tableau{T}) where {T} = Tableau{T}(t.name, t.o, t.s, t.a, t.b, t.c; R∞ = t.R∞)
_zero_hats(pt::PartitionedTableau{T}) where {T} =
    PartitionedTableau{T}(pt.name, pt.o, _zero_hats(pt.q), _zero_hats(pt.p); R∞ = pt.R∞)

GeometricBase.tableau(::GaussPRK,   ::Type{T} = Float64) where {T} = PartitionedTableau(TableauGauss(T, 2))
GeometricBase.tableau(::GaussSPRK,  ::Type{T} = Float64) where {T} = SymplecticPartitionedTableau(TableauGauss(T, 2))
GeometricBase.tableau(::GaussPRK0,  ::Type{T} = Float64) where {T} = _zero_hats(PartitionedTableau(TableauGauss(T, 2)))
GeometricBase.tableau(::GaussSPRK0, ::Type{T} = Float64) where {T} = _zero_hats(SymplecticPartitionedTableau(TableauGauss(T, 2)))

"""
Partitioned Gauss(2) midpoint-rule variants (all symplectic). They differ only in
implementation details: symplectic-by-construction (`SPRK`) versus by-duplication (`PRK`),
and whether the rounding-error compensation coefficients â/b̂/ĉ are retained or zeroed.
"""
const MIDPOINT_METHODS = MethodSpec[
    MethodSpec("PRK Gauss(2)",            GaussPRK(),   true),
    MethodSpec("SPRK Gauss(2)",           GaussSPRK(),  true),
    MethodSpec("PRK Gauss(2), â=b̂=ĉ=0",  GaussPRK0(),  true),
    MethodSpec("SPRK Gauss(2), â=b̂=ĉ=0", GaussSPRK0(), true),
]

"""
All methods (geometric, then non-geometric, then the partitioned Gauss(2) variants). The index
in this vector fixes each method's colour, so grouping/reordering for plots never changes a
method's colour.
"""
const ALL_METHODS = vcat(GEOMETRIC_METHODS, NONGEOMETRIC_METHODS, MIDPOINT_METHODS)

_method_by_name(name) = ALL_METHODS[findfirst(m -> m.name == name, ALL_METHODS)]
_methods_in_order(names...) = MethodSpec[_method_by_name(n) for n in names]

# Grouping used for plotting (orthogonal to the geometric/non-geometric classification). Plots
# are produced once per group so each figure stays readable. The listed order is the plotting
# (draw and legend) order within each group.
"""
The Euler-type methods (both symplectic and non-symplectic).
"""
const EULER_METHODS = _methods_in_order(
    "Symplectic Euler A", "Symplectic Euler B", "Explicit Euler", "Implicit Euler")

"""
The remaining (midpoint / trapezoidal / higher-order) methods.
"""
const OTHER_METHODS = _methods_in_order(
    "RK4", "Explicit Midpoint", "Implicit Midpoint", "Crank-Nicolson")

"""
The default method groups, as `label => methods` pairs; each plotting routine produces one figure
per group. Scripts may pass their own `groups` to the plotting routines (e.g. the Lotka–Volterra
examples use a single variational-integrator group).
"""
const METHOD_GROUPS = [
    "euler"    => EULER_METHODS,
    "other"    => OTHER_METHODS,
    "midpoint" => MIDPOINT_METHODS,
]

# --- Lotka–Volterra (degenerate Lagrangian) variational-integrator comparison ----------------
# The Lotka–Volterra systems are degenerate Lagrangian systems posed as IODE/LODE problems, on
# which the explicit / symplectic-Euler / DIRK methods of the main registry are undefined. The
# comparison there is instead between several flavours of the implicit midpoint rule that DO
# apply to such systems and differ only in implementation details.
#
# `VPRK(Gauss(1))` does not rebuild its tableau at the problem's precision (its `initmethod`
# discards the problem), so it would bake in a Float64 tableau and break precision purity. The
# `GaussVPRK` wrapper mirrors the Runge–Kutta `initmethod(::RKMethod, ::GeometricProblem{…,TT})`
# pattern to construct `VPRK(Gauss(1))` at the run precision `TT`.
struct GaussVPRK <: VPRKMethod end
initmethod(::GaussVPRK, ::GeometricProblem{ST,DT,TT}) where {ST,DT,TT} = VPRK(TableauGauss(TT, 1))
# `isimplicit(::VPRKMethod)` inspects `tableau(method)`, which `GaussVPRK` only provides after
# `initmethod` (at the run precision); declare it directly so the DogLeg solver gate applies.
isimplicit(::GaussVPRK) = true

"""
Variational-integrator comparison for the Lotka–Volterra 2D example: several flavours of the
implicit midpoint rule (all symplectic) that apply to degenerate Lagrangian (IODE/LODE) systems
and differ only in implementation details.
"""
const LV2D_METHODS = MethodSpec[
    MethodSpec("Implicit Midpoint", ImplicitMidpoint(), true),
    MethodSpec("VPRK Gauss(1)",     GaussVPRK(),        true),
    MethodSpec("PMVI Midpoint",     PMVImidpoint(),     true),
    MethodSpec("CMDVI",             CMDVI(),            true),
]

"""
Variational-integrator comparison for the Lotka–Volterra 4D (Lagrangian) example. Same set as
[`LV2D_METHODS`](@ref) but without `CMDVI`, which does not converge on the 4D degenerate Lagrangian
(its Newton/DogLeg iterate leaves the positive orthant and the solve breaks down).
"""
const LV4D_METHODS = filter(m -> m.name != "CMDVI", LV2D_METHODS)

"""
Plotting groups (one `variational` group each) for the Lotka–Volterra examples.
"""
const LV2D_GROUPS = ["variational" => LV2D_METHODS]
const LV4D_GROUPS = ["variational" => LV4D_METHODS]

# Human-readable figure subtitle for a group label; unknown labels fall back to the label itself.
const _GROUP_TITLES = Dict(
    "euler"       => "Euler methods",
    "other"       => "other methods",
    "midpoint"    => "Gauss(2) midpoint variants",
    "variational" => "implicit-midpoint variational integrators",
)
_group_title(label) = get(_GROUP_TITLES, label, label)
