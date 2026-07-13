# Harmonic Oscillator

The harmonic oscillator ``H(q,p) = p^2/(2m) + k q^2 / 2`` is the reference test case: it is linear,
it has a closed-form solution (so the solution error is measured against the **analytic**
solution), and it is where the type-purity pipeline was first verified. All twelve methods (the
eight Euler/Runge–Kutta methods plus the four partitioned Gauss(2) variants) run at all three
precisions and pass the precision-purity gate.

## Short scenario (Δt = 0.1, t ≤ 100)

### Energy error

![Energy error, Euler methods](figures/harmonic_oscillator_energy_error_dt_0.1_euler.png)

![Energy error, other methods](figures/harmonic_oscillator_energy_error_dt_0.1_other.png)

![Energy error, Gauss(2) variants](figures/harmonic_oscillator_energy_error_dt_0.1_gauss2.png)

The qualitative picture is textbook: **explicit Euler** grows without bound, **implicit Euler**
dissipates towards a constant relative error of order one, and the **symplectic Euler** methods keep
a *bounded, oscillating* energy error. Among the higher-order methods, implicit midpoint and
Crank–Nicolson nearly conserve energy, with their noise floor set by the precision
(≈ `1e-2` for Float16, `1e-6` for Float32, `1e-15` for Float64), while explicit midpoint drifts and
RK4 sits in between. The four partitioned Gauss(2) variants differ only in implementation detail —
symplectic-by-construction (`SPRK`) versus by-duplication (`PRK`), and whether the rounding-error
compensation coefficients ``â, b̂, ĉ`` are retained or zeroed; on this linear system the four are
almost indistinguishable, and the differences grow on the nonlinear problems.

### Solution error

![Solution error, Euler methods](figures/harmonic_oscillator_solution_error_dt_0.1_euler.png)

![Solution error, other methods](figures/harmonic_oscillator_solution_error_dt_0.1_other.png)

![Solution error, Gauss(2) variants](figures/harmonic_oscillator_solution_error_dt_0.1_gauss2.png)

Measured against the analytic solution, the trajectory error mirrors the energy behaviour: the
symplectic Euler methods keep a *bounded, oscillating* error (≈ `1e-2`), explicit Euler grows and
implicit Euler saturates near order one, while the higher-order and Gauss(2) methods stay well below
the Euler pair. Reduced precision mainly raises the noise floor rather than changing the ranking.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/harmonic_oscillator_solution_dt_0.1_euler.png)

![Phase-space trajectory, other methods](figures/harmonic_oscillator_solution_dt_0.1_other.png)

![Phase-space trajectory, Gauss(2) variants](figures/harmonic_oscillator_solution_dt_0.1_gauss2.png)

In phase space the behaviour is unmistakable: the symplectic methods stay on a closed (slightly
deformed) ellipse — the level set of a nearby modified Hamiltonian — the reference is the exact
circle, explicit Euler spirals **outward**, and implicit Euler spirals **inward**.

## Long scenario (Δt = 1, t ≤ 10 000)

### Energy error

![Energy error, Euler methods](figures/harmonic_oscillator_energy_error_dt_1.0_euler.png)

![Energy error, other methods](figures/harmonic_oscillator_energy_error_dt_1.0_other.png)

![Energy error, Gauss(2) variants](figures/harmonic_oscillator_energy_error_dt_1.0_gauss2.png)

At the coarse step and long horizon the contrast is dramatic: explicit Euler and explicit midpoint
diverge exponentially (reaching ≈ `1e300` in Float64, clipped at the plot's `1e5` ceiling), while
the symplectic methods and the implicit midpoint / Crank–Nicolson rules — and the partitioned
Gauss(2) variants — remain bounded over the *entire* ``10^4`` time units. In Float16 the implicit
methods fail once the time grid can no longer resolve `Δt` (see [Findings](@ref)).

### Solution error

![Solution error, Euler methods](figures/harmonic_oscillator_solution_error_dt_1.0_euler.png)

![Solution error, other methods](figures/harmonic_oscillator_solution_error_dt_1.0_other.png)

![Solution error, Gauss(2) variants](figures/harmonic_oscillator_solution_error_dt_1.0_gauss2.png)

Over the long horizon the divergent explicit methods leave the analytic reference entirely, whereas
the bounded symplectic, implicit-midpoint and Gauss(2) methods retain a usable trajectory error for
the full ``10^4`` time units at Float32/Float64.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/harmonic_oscillator_solution_dt_1.0_euler.png)

![Phase-space trajectory, other methods](figures/harmonic_oscillator_solution_dt_1.0_other.png)

![Phase-space trajectory, Gauss(2) variants](figures/harmonic_oscillator_solution_dt_1.0_gauss2.png)
