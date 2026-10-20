"""Find an allowed energy of a particle by solving the transcendental
equation  k + tan(k w) k' = 0  with Newton's method.

    k  = sqrt(2 m c^2 E) / hc,   k' = sqrt(2 m c^2 (E - V0)) / hc
Energies in eV, lengths in the units implied by hc = 197e6 eV (x units).
"""
import argparse

import numpy as np

MC2 = 511e3   # electron rest energy [eV]
HC = 197e6    # [eV * length unit]


def energy(v0, width, guess=50.0, tol=1e-12, max_iter=100, return_history=False):
    """Newton iteration for E starting at `guess` (needs guess > v0).

    Uses a central finite-difference derivative (scipy.misc.derivative, used in
    the original, was removed in SciPy 1.12) and a convergence test on the step
    size (the original `|f| < 1e-15` test could loop forever).
    With return_history=True also returns the list of |Newton step| per iteration.
    """
    def f(en):
        k = np.sqrt(2 * MC2 * en) / HC
        kp = np.sqrt(2 * MC2 * (en - v0)) / HC
        return k + np.tan(k * width) * kp

    def df(en, dx=1e-6):
        return (f(en + dx) - f(en - dx)) / (2 * dx)

    x = guess
    steps = []
    for _ in range(max_iter):
        step = f(x) / df(x)
        x -= step
        steps.append(abs(step))
        if abs(step) < tol * max(1.0, abs(x)):
            return (x, steps) if return_history else x
    raise RuntimeError("Newton iteration did not converge")


def main(argv=None):
    p = argparse.ArgumentParser(description="Newton's method for an allowed energy of a piecewise-constant potential.")
    p.add_argument("--v0", type=float, default=10.0, help="well depth V0 in eV (default 10)")
    p.add_argument("--width", type=float, default=10e5, help="well width (default 1e6)")
    p.add_argument("--guess", type=float, default=50.0, help="starting energy in eV, must exceed V0 (default 50)")
    args = p.parse_args(argv)
    print(f"E = {energy(args.v0, args.width, args.guess):.10f} eV")


if __name__ == "__main__":
    main()
