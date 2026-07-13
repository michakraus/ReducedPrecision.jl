# Lotka–Volterra 2D

The Lotka–Volterra predator–prey system is a **degenerate Lagrangian** system: its Lagrangian is
linear in the velocities, so it has no ordinary Hamiltonian phase-space form and is posed instead as
an implicit / degenerate-Lagrangian (IODE / LODE) problem. This example uses the
`LotkaVolterra2dSingular` module — the *singular* gauge of the Lagrangian, which is the one the
`CMDVI` degenerate variational integrator requires.

Because it is degenerate, the explicit, symplectic-Euler and diagonally-implicit Runge–Kutta methods
of the other examples **do not apply**. The comparison here is instead between four flavours of the
**implicit midpoint rule** that *do* apply to such systems and differ only in implementation detail:

- **Implicit Midpoint** — the plain 1-stage Gauss / midpoint Runge–Kutta rule;
- **VPRK Gauss(1)** — the variational partitioned Runge–Kutta form of the same rule;
- **PMVI Midpoint** — the position–momentum variational integrator;
- **CMDVI** — the continuous midpoint degenerate variational integrator.

There is no closed-form solution; the solution error is measured against a Float64 `Gauss(8)`
reference (in the coarse scenario computed at the fine step and subsampled to the output grid). All
four methods are type-pure at every precision.

## Short scenario (Δt = 0.01, t ≤ 10)

### Energy error

![Energy error, variational integrators](figures/lotka_volterra_2d_energy_error_dt_0.01_variational.png)

At Float32 and Float64 all four midpoint variants conserve energy well and trace visibly different
fine structure — CMDVI in particular follows a distinct curve. At **Float16** the `VPRK` and `PMVI`
solves hit the same half-precision `NaN`-direction breakdown seen on the other stiff problems, while
implicit midpoint and CMDVI still run.

### Solution error

![Solution error, variational integrators](figures/lotka_volterra_2d_solution_error_dt_0.01_variational.png)

Against the `Gauss(8)` reference the four variants keep a small trajectory error at Float32/Float64,
with CMDVI again separating from the others; at Float16 only the surviving methods contribute.

### Configuration-space trajectory

![Trajectory, variational integrators](figures/lotka_volterra_2d_solution_dt_0.01_variational.png)

## Coarse scenario (Δt = 0.1, t ≤ 100)

### Energy error

![Energy error, variational integrators](figures/lotka_volterra_2d_energy_error_dt_0.1_variational.png)

At the coarser step and longer horizon the differences between the variational integrators, and the
effect of precision on their long-time energy behaviour, are accentuated.

### Solution error

![Solution error, variational integrators](figures/lotka_volterra_2d_solution_error_dt_0.1_variational.png)

Over the longer horizon the trajectory error against the reference likewise separates the variants,
exposing how each variational integrator accumulates error at the coarse step.

### Configuration-space trajectory

![Trajectory, variational integrators](figures/lotka_volterra_2d_solution_dt_0.1_variational.png)
