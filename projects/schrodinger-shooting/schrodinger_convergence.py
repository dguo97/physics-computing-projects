"""Accuracy check: how does the eigenenergy change as the RK4 step is refined?

No exact answer is needed: we look at |E(h) - E(h/2)| as h shrinks.  Because the
potential is piecewise constant and each RK4 step sees a single value of V, the
solver keeps its 4th-order accuracy (changes shrink ~16x per halving of h).
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import schrodinger as sc

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

STEP_COUNTS = [250, 500, 1000, 2000]


def energies(step_counts=STEP_COUNTS):
    """Eigenenergy found with each step count."""
    return np.array([sc.find_energy(steps=n) for n in step_counts])


def convergence_order(step_counts=STEP_COUNTS):
    """Return (h, |E(h) - E(h/2)|, fitted slope on a log-log plot)."""
    e = energies(step_counts)
    h = sc.X_END / np.array(step_counts, dtype=float)
    diff = np.abs(np.diff(e))
    return h[:-1], diff, np.polyfit(np.log(h[:-1]), np.log(diff), 1)[0]


def main():
    e = energies()
    h = sc.X_END / np.array(STEP_COUNTS, dtype=float)
    diff = np.abs(np.diff(e))              # |E(h) - E(h/2)|
    slope = np.polyfit(np.log(h[:-1]), np.log(diff), 1)[0]
    fig, ax = plt.subplots()
    ax.loglog(h[:-1], diff, "o-", label=f"|E(h) - E(h/2)| (slope {slope:.2f})")
    ax.loglog(h[:-1], diff[0] * (h[:-1] / h[0]) ** 4, "k--", label=r"$\propto h^4$")
    ax.set(xlabel="RK4 step size h", ylabel="change in eigenenergy [eV]",
           title="Schrodinger shooting: convergence with step size")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("energies:", e, f"\nobserved order: {slope:.2f}")


if __name__ == "__main__":
    main()
