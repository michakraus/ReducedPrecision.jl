# Pendulum

The mathematical pendulum ``H(q,p) = p^2/2 + \cos q`` is nonlinear and has no closed-form solution,
so the solution error is measured against a high-precision **Float64 `Gauss(8)`** reference (in the
coarse scenario computed at the fine step and subsampled to the output grid). All methods are
type-pure at all precisions.

## Short scenario (Δt = 0.1, t ≤ 100)

### Energy error

![Energy error, Euler methods](figures/pendulum_energy_error_dt_0.1_euler.png)

![Energy error, other methods](figures/pendulum_energy_error_dt_0.1_other.png)

As for the harmonic oscillator, the symplectic methods keep the energy error bounded while explicit
Euler grows and implicit Euler dissipates. Implicit midpoint and Crank–Nicolson track the precision
floor (≈ `1e-2` / `1e-5` / `1e-11` for Float16 / Float32 / Float64).

### Partitioned Gauss(2) variants

![Energy error, Gauss(2) variants](figures/pendulum_energy_error_dt_0.1_gauss2.png)

The four partitioned-Gauss(2) variants (symplectic-by-construction vs. by-duplication, with and
without the rounding-compensation coefficients ``â, b̂, ĉ``) are compared here; on this nonlinear
system the implementation differences become visible in the energy-error fine structure.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/pendulum_solution_dt_0.1_euler.png)

## Long scenario (Δt = 1, t ≤ 10 000)

![Energy error, Euler methods](figures/pendulum_energy_error_dt_1.0_euler.png)

![Energy error, other methods](figures/pendulum_energy_error_dt_1.0_other.png)

The symplectic and midpoint/trapezoidal methods remain bounded over the full horizon; the
non-symplectic Euler and explicit-midpoint methods drift. In Float16 the implicit methods fail on
the long-horizon time-grid saturation.

![Phase-space trajectory, Euler methods](figures/pendulum_solution_dt_1.0_euler.png)
