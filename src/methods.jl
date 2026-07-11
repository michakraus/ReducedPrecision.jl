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

"""
All methods (geometric followed by non-geometric). The index in this vector fixes each
method's colour, so grouping/reordering for plots never changes a method's colour.
"""
const ALL_METHODS = vcat(GEOMETRIC_METHODS, NONGEOMETRIC_METHODS)

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
The two method groups, as `label => methods` pairs; each plotting routine produces one figure
per group.
"""
const METHOD_GROUPS = ["euler" => EULER_METHODS, "other" => OTHER_METHODS]

_group_title(label) = label == "euler" ? "Euler methods" : "other methods"
