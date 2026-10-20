"""Accuracy check: drift of the conserved Jacobi constant versus RK4 step size.

The exact solution conserves the Jacobi constant, so its drift along the orbit
measures the integrator's error without needing a reference solution.
RK4 is 4th order, so the drift should shrink as h^4.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import three_body as tb

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

STEP_COUNTS = [40_000, 80_000, 160_000, 320_000]


def convergence_order(step_counts=STEP_COUNTS):
    """Return (h, Jacobi drift, fitted slope of log(drift) vs log(h))."""
    h = tb.T_END / np.array(step_counts, dtype=float)
    drift = np.array([np.ptp(tb.jacobi_constant(tb.integrate(n)[1])) for n in step_counts])
    slope = np.polyfit(np.log(h), np.log(drift), 1)[0]
    return h, drift, slope


def main():
    h, drift, slope = convergence_order()
    fig, ax = plt.subplots()
    ax.loglog(h, drift, "o-", label=f"Jacobi-constant drift (slope {slope:.2f})")
    ax.loglog(h, drift[0] * (h / h[0]) ** 4, "k--", label=r"$\propto h^4$")
    ax.set(xlabel="step size h", ylabel="max - min of Jacobi constant",
           title="Three-body RK4: energy-like drift vs step size")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("drift:", drift, f"\nfitted order: {slope:.2f} (theory: 4)")


if __name__ == "__main__":
    main()
