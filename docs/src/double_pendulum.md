# Double Pendulum

The double pendulum is a stiff, chaotic, two-degree-of-freedom system with dimensional parameters
(``g = 9.80665``). It is generated symbolically by `EulerLagrange` and has no `::Type{T}`
constructor, so the initial conditions, timespan, timestep and parameters are built at precision `T`
by hand. There is no closed-form solution; the reference is a Float64 `Gauss(8)` run. The trajectory
plots use the **configuration space** ``(\theta_1, \theta_2)``.

## Short scenario (Δt = 0.01, t ≤ 10)

### Energy error

![Energy error, Euler methods](figures/double_pendulum_energy_error_dt_0.01_euler.png)

![Energy error, other methods](figures/double_pendulum_energy_error_dt_0.01_other.png)

![Energy error, Gauss(2) variants](figures/double_pendulum_energy_error_dt_0.01_gauss2.png)

At Float32 and Float64 all methods run and the usual ordering holds — implicit midpoint and
Crank–Nicolson keep the energy error far below the Euler methods. At **Float16 the three implicit
methods fail** with a `NaN` in the nonlinear-solver direction: half precision is simply inadequate
for the implicit solves on this stiff system (the explicit and symplectic Euler methods still run).
Switching the nonlinear solver from `Newton` to the trust-region `DogLeg` (the default here) does
not rescue these Float16 solves — the breakdown is a genuine property of half precision, not of the
solver. The four partitioned-Gauss(2) variants share the fate of the other implicit methods at
Float16 but are informative at Float32/Float64, where the symplectic-vs-duplicated tableau choice
and the rounding-compensation coefficients ``â, b̂, ĉ`` produce visibly different energy-error fine
structure.

### Solution error

![Solution error, Euler methods](figures/double_pendulum_solution_error_dt_0.01_euler.png)

![Solution error, other methods](figures/double_pendulum_solution_error_dt_0.01_other.png)

![Solution error, Gauss(2) variants](figures/double_pendulum_solution_error_dt_0.01_gauss2.png)

Against the `Gauss(8)` reference every method tracks closely until the chaotic divergence sets in,
after which all trajectory errors saturate at order one; at Float16 the surviving methods diverge
noticeably earlier than at Float32/Float64.

### Configuration-space trajectory

![Configuration-space trajectory, Euler methods](figures/double_pendulum_solution_dt_0.01_euler.png)

![Configuration-space trajectory, other methods](figures/double_pendulum_solution_dt_0.01_other.png)

![Configuration-space trajectory, Gauss(2) variants](figures/double_pendulum_solution_dt_0.01_gauss2.png)

The methods track the reference until the chaotic divergence sets in; at Float16 the surviving
methods depart from the reference noticeably earlier.

## Coarse scenario (Δt = 0.1, t ≤ 10)

### Energy error

![Energy error, Euler methods](figures/double_pendulum_energy_error_dt_0.1_euler.png)

![Energy error, other methods](figures/double_pendulum_energy_error_dt_0.1_other.png)

![Energy error, Gauss(2) variants](figures/double_pendulum_energy_error_dt_0.1_gauss2.png)

At the ten-times-coarser step `Δt = 0.1` (over the same `t ≤ 10` horizon as the fine run) every
solve now stays stable — with the line-search `Newton`/`Backtracking` solver capped at 100
iterations, nothing trips the divergence guard. The geometric methods keep the smallest energy
error: symplectic Euler A/B, implicit midpoint and the partitioned Gauss(2) variants stay bounded
around `1e-2`–`1e-1`, while the non-geometric methods drift up toward order one (explicit Euler
worst, then explicit midpoint, RK4 and implicit Euler). The only outright failure is Crank–Nicolson
in `Float16` (a NaN in the Newton direction, guarded and skipped — hence its absence from the
`Float16` panel). Reduced precision mainly raises the error floor; the qualitative ranking is the
same across `Float16`/`Float32`/`Float64`. As always for this chaotic system the fine `Δt = 0.01`
run is the more informative one.

### Solution error

![Solution error, Euler methods](figures/double_pendulum_solution_error_dt_0.1_euler.png)

![Solution error, other methods](figures/double_pendulum_solution_error_dt_0.1_other.png)

![Solution error, Gauss(2) variants](figures/double_pendulum_solution_error_dt_0.1_gauss2.png)

At the coarse step the trajectory error against the `Gauss(8)` reference grows quickly for every
method as the chaotic orbits separate, with the geometric and Gauss(2) methods retaining the
smallest error before the divergence dominates.

### Configuration-space trajectory

![Configuration-space trajectory, Euler methods](figures/double_pendulum_solution_dt_0.1_euler.png)

![Configuration-space trajectory, other methods](figures/double_pendulum_solution_dt_0.1_other.png)

![Configuration-space trajectory, Gauss(2) variants](figures/double_pendulum_solution_dt_0.1_gauss2.png)
