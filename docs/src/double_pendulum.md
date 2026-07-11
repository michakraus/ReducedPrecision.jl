# Double Pendulum

The double pendulum is a stiff, chaotic, two-degree-of-freedom system with dimensional parameters
(``g = 9.80665``). It is generated symbolically by `EulerLagrange` and has no `::Type{T}`
constructor, so the initial conditions, timespan, timestep and parameters are built at precision `T`
by hand. There is no closed-form solution; the reference is a Float64 `Gauss(8)` run. The trajectory
plots use the **configuration space** ``(\theta_1, \theta_2)``.

## Short scenario (Δt = 0.01, t ≤ 10)

### Energy error

![Energy error, Euler methods](figures/double_pendulum_energy_error_dt=0.01_euler.png)

![Energy error, other methods](figures/double_pendulum_energy_error_dt=0.01_other.png)

At Float32 and Float64 all methods run and the usual ordering holds — implicit midpoint and
Crank–Nicolson keep the energy error far below the Euler methods. At **Float16 the three implicit
methods fail** with a `NaN` in the nonlinear-solver direction: half precision is simply inadequate
for the implicit solves on this stiff system (the explicit and symplectic Euler methods still run).
Switching the nonlinear solver from `Newton` to the trust-region `DogLeg` (the default here) does
not rescue these Float16 solves — the breakdown is a genuine property of half precision, not of the
solver.

### Partitioned Gauss(2) midpoint variants

![Energy error, Gauss(2) midpoint variants](figures/double_pendulum_energy_error_dt=0.01_midpoint.png)

The four partitioned-Gauss(2) variants share the fate of the other implicit methods at Float16 but
are informative at Float32/Float64, where the symplectic-vs-duplicated tableau choice and the
rounding-compensation coefficients ``â, b̂, ĉ`` produce visibly different energy-error fine structure.

### Configuration-space trajectory

![Configuration-space trajectory, other methods](figures/double_pendulum_solution_dt=0.01_other.png)

The methods track the reference until the chaotic divergence sets in; at Float16 the surviving
methods depart from the reference noticeably earlier.

## Coarse scenario (Δt = 0.1, t ≤ 1000)

![Energy error, Euler methods](figures/double_pendulum_energy_error_dt=0.1_euler.png)

![Energy error, other methods](figures/double_pendulum_energy_error_dt=0.1_other.png)

![Energy error, Gauss(2) midpoint variants](figures/double_pendulum_energy_error_dt=0.1_midpoint.png)

At the ten-times-coarser step `Δt = 0.1` the chaotic double pendulum is badly under-resolved and the
relative energy error is of order one for essentially every method. Counter-intuitively it is the
**implicit** solves (implicit midpoint, Crank–Nicolson, and the partitioned Gauss(2) variants) that
spike and diverge early — their nonlinear iterations blow up on the stiff coarse-step problem and
the guard truncates them — while the explicit RK methods merely hover around `1e0` for the whole
horizon (tracking a wrong-but-bounded trajectory). At `Float16` the implicit solves fail outright.
Reduced precision makes little qualitative difference at this step; the short-step run is the
informative one for the double pendulum.
