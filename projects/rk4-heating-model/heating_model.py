"""Solve dy/dx = a*exp(-2x) - b*(y^4 - 1), y(0)=0, with 4th-order Runge-Kutta.

Models a body heated by a decaying source while radiating to its surroundings.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)


def ff(x, y, apar, bpar):
    """Right-hand side of the ODE."""
    return apar * np.exp(-2 * x) - bpar * (y**4 - 1)


def solve_rk4(apar, bpar, xmax, steps=100_000):
    """RK4 on [0, xmax] with `steps` equal steps.

    Returns (xs, ys, y_end): the grid and solution at x_i = i*h (i < steps),
    plus the value at x = xmax itself.
    """
    h = xmax / steps
    xs = np.linspace(0.0, xmax, steps, endpoint=False)
    ys = np.empty(steps)
    y = 0.0
    for i, x in enumerate(xs):  # preallocated (np.append in a loop was O(n^2))
        ys[i] = y
        k1 = h * ff(x, y, apar, bpar)
        k2 = h * ff(x + 0.5 * h, y + 0.5 * k1, apar, bpar)
        k3 = h * ff(x + 0.5 * h, y + 0.5 * k2, apar, bpar)
        k4 = h * ff(x + h, y + k3, apar, bpar)
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return xs, ys, y


def temperature(apar, bpar, xmax, steps=100_000):
    """Integrate on [0, xmax]; save the plot; return (max y, x at max y)."""
    xs, ys, _ = solve_rk4(apar, bpar, xmax, steps)

    fig, ax = plt.subplots()
    ax.plot(xs, ys)
    ax.set(title="Temperature against time", xlabel="time", ylabel="temperature")
    fig.savefig(OUTPUT / "rk4_temperature.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    i = int(np.argmax(ys))
    return ys[i], xs[i]


def main(argv=None):
    p = argparse.ArgumentParser(description="RK4 heating/cooling model: dy/dx = a e^(-2x) - b (y^4 - 1).")
    p.add_argument("--a", type=float, default=5.0, help="heating strength (default 5)")
    p.add_argument("--b", type=float, default=0.5, help="radiative cooling strength (default 0.5)")
    p.add_argument("--xmax", type=float, default=5.0, help="end time (default 5)")
    p.add_argument("--steps", type=int, default=100_000, help="RK4 steps (default 100000)")
    args = p.parse_args(argv)
    ymax, xpeak = temperature(args.a, args.b, args.xmax, args.steps)
    print(f"max temperature {ymax:.5f} at t = {xpeak:.5f}")


if __name__ == "__main__":
    main()
