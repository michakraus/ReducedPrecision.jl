# Pendulum

The mathematical pendulum ``H(q,p) = p^2/2 + \cos q`` is nonlinear and has no closed-form solution,
so the solution error is measured against a high-precision **Float64 `Gauss(8)`** reference on the
same time grid. All methods are type-pure at all precisions.

## Short scenario (Δt = 0.1, t ≤ 100)

### Energy error

![Energy error, Euler methods](figures/pendulum_energy_error_euler.png)

![Energy error, other methods](figures/pendulum_energy_error_other.png)

As for the harmonic oscillator, the symplectic methods keep the energy error bounded while explicit
Euler grows and implicit Euler dissipates. Implicit midpoint and Crank–Nicolson track the precision
floor (≈ `1e-2` / `1e-5` / `1e-11` for Float16 / Float32 / Float64).

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/pendulum_solution_euler.png)

## Long scenario (Δt = 1, t ≤ 10 000)

![Energy error, Euler methods](figures/pendulum_longtime_energy_error_euler.png)

![Energy error, other methods](figures/pendulum_longtime_energy_error_other.png)

The symplectic and midpoint/trapezoidal methods remain bounded over the full horizon; the
non-symplectic Euler and explicit-midpoint methods drift. In Float16 the implicit methods fail on
the long-horizon time-grid saturation.

![Phase-space trajectory, Euler methods](figures/pendulum_longtime_solution_euler.png)
