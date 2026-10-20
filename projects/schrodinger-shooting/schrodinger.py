"""Shooting-method solution of the 1D Schrodinger equation for a particle in a
box of width 2a whose floor is 0 on [0, a] and V0 on [a, 2a]
(infinite walls at x=0 and x=2a).

  - RK4 integrates psi'' = (V - E) psi / (hbar^2 / 2m)
  - bisection finds the E that gives psi(2a) = 0
  - Simpson's rule normalises psi
Energies in eV, lengths in the units implied by hc = 197e6 eV.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import simpson

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

A = 0.75e5
V0 = 5.0
MC2 = 511e3
HC = 197e6
CONST = HC**2 / (2 * MC2)   # hbar^2 / 2m
STEPS = 1000
X_END = 2 * A


def potential(x, v0=V0):
    """Potential energy [eV] at position x."""
    if 0 <= x <= A:
        return 0.0
    if A < x <= X_END:
        return v0
    return np.inf


def _rhs(state, v, energy):
    """Derivatives of (psi, psi') with the potential `v` held constant."""
    psi, dpsi = state
    return np.array([dpsi, (v - energy) * psi / CONST])


def solve(energy, steps=STEPS, v0=V0):
    """Integrate from x=0 with psi=0, psi'=1.  Returns (x, psi) on steps+1 points."""
    h = X_END / steps
    x = np.linspace(0.0, X_END, steps + 1)
    psi = np.empty(steps + 1)
    state = np.array([0.0, 1.0])
    psi[0] = state[0]
    for i in range(steps):
        # V is piecewise constant, so evaluate it once per step at the midpoint:
        # a step then never straddles the jump at x = a (it sits on a grid point
        # for an even step count) and RK4 keeps its 4th-order accuracy.
        v = potential(x[i] + 0.5 * h, v0)
        k1 = h * _rhs(state, v, energy)
        k2 = h * _rhs(state + 0.5 * k1, v, energy)
        k3 = h * _rhs(state + 0.5 * k2, v, energy)
        k4 = h * _rhs(state + k3, v, energy)
        state = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        psi[i + 1] = state[0]
    return x, psi


def find_energy(e_low=50.0, e_high=100.0, target=1e-10, steps=STEPS, v0=V0):
    """Bisect for the energy where psi(2a) = 0 inside [e_low, e_high]."""
    f_low = solve(e_low, steps, v0)[1][-1]
    if f_low * solve(e_high, steps, v0)[1][-1] > 0:
        raise ValueError("psi(2a) has the same sign at both ends: no eigenvalue bracketed")
    while abs(e_high - e_low) > target:
        mid = 0.5 * (e_low + e_high)
        f_mid = solve(mid, steps, v0)[1][-1]
        if f_low * f_mid > 0:
            e_low, f_low = mid, f_mid
        else:
            e_high = mid
    return 0.5 * (e_low + e_high)


def schroedinger(e_low=50.0, e_high=100.0, steps=STEPS, v0=V0):
    """Save normalised wave function and probability plots; return the energy."""
    energy = find_energy(e_low, e_high, steps=steps, v0=v0)
    x, psi = solve(energy, steps, v0)
    psi /= np.sqrt(simpson(psi**2, x=x))   # normalise
    prob = psi**2

    fig, ax = plt.subplots()
    ax.plot(x, psi)
    ax.set(title="Normalised wave function", xlabel="x", ylabel="psi(x)")
    fig.savefig(OUTPUT / "schrodinger_wavefunction.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(x, prob)
    ax.set(title="Normalised probability density", xlabel="x", ylabel="|psi(x)|^2")
    fig.savefig(OUTPUT / "schrodinger_probability.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    assert abs(simpson(prob, x=x) - 1) < 1e-9
    return energy


def main(argv=None):
    p = argparse.ArgumentParser(description="Shooting-method eigenenergy for a box with a step floor.")
    p.add_argument("--v0", type=float, default=V0, help="height of the step in eV (default 5)")
    p.add_argument("--e-low", type=float, default=50.0, help="lower end of the energy bracket in eV (default 50)")
    p.add_argument("--e-high", type=float, default=100.0, help="upper end of the energy bracket in eV (default 100)")
    p.add_argument("--steps", type=int, default=STEPS, help="RK4 steps (default 1000)")
    args = p.parse_args(argv)
    print(f"E = {schroedinger(args.e_low, args.e_high, args.steps, args.v0):.8f} eV")


if __name__ == "__main__":
    main()
