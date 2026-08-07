# Pendulum

The mathematical pendulum ``H(q,p) = p^2/2 + \cos q`` is nonlinear and has no closed-form solution,
so the solution error is measured against a high-precision **Float64 `Gauss(8)`** reference (in the
coarse scenario computed at the fine step and subsampled to the output grid). All methods are
type-pure at all four precisions and run the full horizon.

## Short scenario (Δt = 0.1, t ≤ 1000)

### Energy error

![Energy error, Euler methods](figures/pendulum_energy_error_dt_0.1_euler.png)

![Energy error, other methods](figures/pendulum_energy_error_dt_0.1_other.png)

![Energy error, Gauss(2) variants](figures/pendulum_energy_error_dt_0.1_gauss2.png)

As for the harmonic oscillator, the symplectic methods keep the energy error bounded while explicit
Euler grows and implicit Euler dissipates. The two orders separate instructively here. **Implicit
midpoint** sits at ≈ `2e-1` / `4e-2` / `6e-4` / `6e-4` for BFloat16 / Float16 / Float32 / Float64 —
the last two are *equal*, because at `Δt = 0.1` this second-order rule is already truncation-limited
on the nonlinear pendulum, so precision beyond Float32 buys it nothing. **Implicit RK4**, being
fourth-order, is still round-off-limited at those precisions and keeps improving (≈ `2e-5` at Float32,
`3e-7` at Float64). Which of the two floors a curve is sitting on — arithmetic or truncation — is
therefore method-dependent, and reading a precision study requires knowing which. The four
partitioned-Gauss(2) variants (symplectic-by-construction vs. by-duplication, with and without the
rounding-compensation coefficients ``â, b̂, ĉ``) become distinguishable here: on this nonlinear system
the implementation differences are visible in the energy-error fine structure.

### Solution error

![Solution error, Euler methods](figures/pendulum_solution_error_dt_0.1_euler.png)

![Solution error, other methods](figures/pendulum_solution_error_dt_0.1_other.png)

![Solution error, Gauss(2) variants](figures/pendulum_solution_error_dt_0.1_gauss2.png)

Against the `Gauss(8)` reference the trajectory error accumulates fastest for the non-symplectic
Euler and explicit-midpoint methods, while the symplectic, implicit-midpoint and Gauss(2) methods
keep the error small; reduced precision raises the floor without reordering the methods.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/pendulum_solution_dt_0.1_euler.png)

![Phase-space trajectory, other methods](figures/pendulum_solution_dt_0.1_other.png)

![Phase-space trajectory, Gauss(2) variants](figures/pendulum_solution_dt_0.1_gauss2.png)

## Long scenario (Δt = 1, t ≤ 10 000)

### Energy error

![Energy error, Euler methods](figures/pendulum_energy_error_dt_1.0_euler.png)

![Energy error, other methods](figures/pendulum_energy_error_dt_1.0_other.png)

![Energy error, Gauss(2) variants](figures/pendulum_energy_error_dt_1.0_gauss2.png)

The symplectic and midpoint/trapezoidal methods — together with the partitioned Gauss(2) variants —
remain bounded over the full horizon; the non-symplectic Euler and explicit-midpoint methods drift.
Both half precisions carry the full `t ≤ 10 000` horizon here, far past where a `T`-typed clock stops
advancing (`t ≈ 2048` in Float16, `t ≈ 256` in BFloat16), which on its own would make the implicit
methods fail outright — see [Time stepping in a local frame](@ref) and [Initial guess](@ref). This scenario also
holds the *only* remaining failure among the four Hamiltonian problems: `Implicit Euler` at `Float16`
throws a `NaN` in the Newton direction. A first-order dissipative method at `Δt = 1` is the least
promising combination in the study, so it is a fitting place for the last one to sit.

### Solution error

![Solution error, Euler methods](figures/pendulum_solution_error_dt_1.0_euler.png)

![Solution error, other methods](figures/pendulum_solution_error_dt_1.0_other.png)

![Solution error, Gauss(2) variants](figures/pendulum_solution_error_dt_1.0_gauss2.png)

Over the long horizon the drifting methods depart from the reference, whereas the bounded geometric
and Gauss(2) methods keep a controlled trajectory error at the precisions where they still run.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/pendulum_solution_dt_1.0_euler.png)

![Phase-space trajectory, other methods](figures/pendulum_solution_dt_1.0_other.png)

![Phase-space trajectory, Gauss(2) variants](figures/pendulum_solution_dt_1.0_gauss2.png)
