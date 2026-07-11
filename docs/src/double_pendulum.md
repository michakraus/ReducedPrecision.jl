# Double Pendulum

The double pendulum is a stiff, chaotic, two-degree-of-freedom system with dimensional parameters
(``g = 9.80665``). It is generated symbolically by `EulerLagrange` and has no `::Type{T}`
constructor, so the initial conditions, timespan, timestep and parameters are built at precision `T`
by hand. There is no closed-form solution; the reference is a Float64 `Gauss(8)` run. The trajectory
plots use the **configuration space** ``(\theta_1, \theta_2)``.

## Short scenario (Δt = 0.01, t ≤ 10)

### Energy error

![Energy error, Euler methods](figures/double_pendulum_energy_error_euler.png)

![Energy error, other methods](figures/double_pendulum_energy_error_other.png)

At Float32 and Float64 all methods run and the usual ordering holds — implicit midpoint and
Crank–Nicolson keep the energy error far below the Euler methods. At **Float16 the three implicit
methods fail** with a `NaN` in the Newton direction: half precision is simply inadequate for the
implicit solves on this stiff system (the explicit and symplectic Euler methods still run).

### Configuration-space trajectory

![Configuration-space trajectory, other methods](figures/double_pendulum_solution_other.png)

The methods track the reference until the chaotic divergence sets in; at Float16 the surviving
methods depart from the reference noticeably earlier.

## Long scenario (Δt = 1, t ≤ 10 000)

![Energy error, Euler methods](figures/double_pendulum_longtime_energy_error_euler.png)

![Energy error, other methods](figures/double_pendulum_longtime_energy_error_other.png)

`Δt = 1` is roughly one natural period per step — far too coarse for this system. **Every implicit
method fails at every precision**, and even the high-order `Gauss(8)` reference does not converge,
so the solution-error and trajectory plots are skipped (only the energy-error plots are produced).
The explicit and symplectic methods blow up almost immediately. This scenario documents the
breakdown; the short-step run is the informative one for the double pendulum.
