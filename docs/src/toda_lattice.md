# Toda Lattice

The Toda lattice is a periodic ``N``-site lattice with nearest-neighbour exponential interactions
``H(q,p) = \tfrac12 p \cdot p + \alpha \sum_n e^{q_n - q_{n+1}}`` — a classic completely integrable
system. Like the double pendulum it is `EulerLagrange`-generated (no `::Type{T}` constructor, inputs
built at precision `T` by hand), and its Hamiltonian additionally takes the lattice size `N`. A
modest lattice (**`N = 16`**) is used so the sweep stays tractable; the default ``N = 200`` soliton
example would make the high-order reference and the implicit solves prohibitively slow. The
trajectory plots use the phase space ``(q_1, p_1)`` of the first lattice site, and the bump initial
condition keeps ``q \in [0,1]`` so the exponentials are well-behaved.

## Short scenario (Δt = 0.1, t ≤ 100)

### Energy error

![Energy error, Euler methods](figures/toda_lattice_energy_error_dt=0.1_euler.png)

![Energy error, other methods](figures/toda_lattice_energy_error_dt=0.1_other.png)

All methods run at **all three precisions** — the bounded initial data makes the Toda lattice
markedly more Float16-friendly than the double pendulum. The now-familiar pattern holds: symplectic
and midpoint/trapezoidal methods keep the energy error bounded (≈ `1e-5`–`1e-6` at Float32/Float64),
while explicit midpoint and RK4 drift and the Euler methods grow/dissipate.

### Partitioned Gauss(2) midpoint variants

![Energy error, Gauss(2) midpoint variants](figures/toda_lattice_energy_error_dt=0.1_midpoint.png)

The four partitioned-Gauss(2) variants are compared on the Toda lattice; the differences between the
symplectic and duplicated tableaus and between keeping or zeroing ``â, b̂, ĉ`` are most visible at the
higher precisions where the energy-error floor is not precision-limited.

### Phase-space trajectory

![Phase-space trajectory, Euler methods](figures/toda_lattice_solution_dt=0.1_euler.png)

The first site traces a quasi-periodic orbit; the symplectic methods stay on it while explicit Euler
spirals out and implicit Euler spirals in, just as for the single oscillator.

## Long scenario (Δt = 1, t ≤ 10 000)

![Energy error, Euler methods](figures/toda_lattice_energy_error_dt=1.0_euler.png)

![Energy error, other methods](figures/toda_lattice_energy_error_dt=1.0_other.png)

The `Gauss(8)` reference converges even at `Δt = 1`, so the full plot set is produced. Implicit
midpoint and Crank–Nicolson keep the energy bounded (≈ `1e-5`) over the whole horizon while explicit
midpoint and RK4 drift; in Float16 the two implicit methods fail on the long-horizon time-grid
saturation.

![Phase-space trajectory, Euler methods](figures/toda_lattice_solution_dt=1.0_euler.png)
