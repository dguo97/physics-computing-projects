# Circuit & Linear-System Solver

Solves a general 4×4 system and a 9-unknown resistor circuit written as Kirchhoff current/voltage equations.

## Physics
**Governing laws.** Kirchhoff's current law ($\sum I=0$ at every junction), Kirchhoff's voltage law ($\sum V=0$ around every loop) and Ohm's law ($V=IR$) are all *linear* in the unknown currents and voltage, so the circuit collapses to

$$A\,\mathbf{x}=\mathbf{b}$$

**Why this method.** For a small dense linear system, a direct LU solve (`np.linalg.solve`) is exact up to rounding and needs no iteration or initial guess. As an independent check, the solution is substituted into loop equations that were *not* put in the matrix; they should vanish.

**What the result means.** A solution that satisfies the extra loop equations (residuals below $10^{-9}$) shows the equations are consistent and the circuit is fully determined.

## Run
```bash
python linear_solver.py
python linear_solver.py --which circuit    # or: system
```
Run any script with `--help` to list its options.

## Result
All loop-check residuals are ~0 (|r| < 1e-9).

## Tested by
`test_linear_systems_solve_their_equations` in [`tests/test_projects.py`](../../tests/test_projects.py).
