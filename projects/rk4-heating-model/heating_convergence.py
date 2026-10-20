"""Accuracy check: the global error of RK4 should fall as h^4.

Compares y(x_max) from our RK4 at a range of step sizes with a high-accuracy
reference from scipy's DOP853 and fits the slope on a log-log plot.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import heating_model as hm

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

A, B, XMAX = 5.0, 0.5, 5.0
STEP_COUNTS = [25, 50, 100, 200, 400, 800, 1600]


def convergence_order(step_counts=STEP_COUNTS):
    """Return (h, error, fitted slope of log(error) vs log(h))."""
    ref = solve_ivp(lambda x, y: hm.ff(x, y[0], A, B), (0, XMAX), [0.0],
                    method="DOP853", rtol=1e-13, atol=1e-15).y[0, -1]
    h = XMAX / np.array(step_counts, dtype=float)
    err = np.array([abs(hm.solve_rk4(A, B, XMAX, n)[2] - ref) for n in step_counts])
    slope = np.polyfit(np.log(h), np.log(err), 1)[0]
    return h, err, slope


def main():
    h, err, slope = convergence_order()
    fig, ax = plt.subplots()
    ax.loglog(h, err, "o-", label=f"RK4 error (fitted slope {slope:.2f})")
    ax.loglog(h, err[0] * (h / h[0]) ** 4, "k--", label=r"$\propto h^4$")
    ax.set(xlabel="step size h", ylabel="|y(x_max) - reference|",
           title="RK4 convergence: error falls as h^4")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"fitted convergence order: {slope:.2f} (theory: 4)")


if __name__ == "__main__":
    main()
