# Kelvin-Helmholtz instability
Numerical simulation of the KH instability with SU2

## Geometry, Boundary and Initial conditions
The Euler equations are soved on a square of edge equal to 1. The computational grid is a uniform structured mesh with 1000x1000 points. Periodic boundary conditions are enforced on all sides of the square. The initial condition is set to trigger the instability as:
$$
(\rho0, u0, v0, p0) = \begin{cases}
(2, -0.5, 0.01 sin(2 \pi x), 2.5) \quad for |y| \leq 0.25 \\
(1, 0.5, 0.01 sin(2 \pi x), 2.5) \quad for |y| > 0.25
\end{cases}
$$
The working fluid is an ideal gas with $\gamma = 1.4$ and $R = 1$.

## Numerical setup
The convective fluxes are discretized using the Approximate Riemann Solver of Roe. The scheme is made second order in space by employing the MUSCL reconstruction. Gradients are computed via the weighted least squares approach, anche the slope of the resconstruction is limited via the Van Albada limiter. A dual time stepping approach is used to evolve the equations in time. The physical time stepping is performed via a BDF2 scheme, making the simulation second order in time. The internal loop of the dual time stepping is integrated with the implicit euler scheme. The internal loop is iterated until the RMS of the residual of the density equation falls below the value of $10^{-7}$. The physical time step used is $0.004$. 
## Simulation result
![](FIG/rhowcont.gif)
