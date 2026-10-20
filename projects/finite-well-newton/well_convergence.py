"""Accuracy check: Newton's method converges quadratically near the root.

Plots the size of each Newton step; once close to the root the number of correct
digits roughly doubles per iteration until it hits the floating-point floor.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import well_energy

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)


def newton_steps(v0=10, width=10e5, guess=50.0):
    """Return (final energy, |step| for each Newton iteration)."""
    return well_energy.energy(v0, width, guess, return_history=True)


def main():
    e, steps = newton_steps()
    fig, ax = plt.subplots()
    ax.semilogy(range(1, len(steps) + 1), steps, "o-")
    ax.set(xlabel="Newton iteration", ylabel="|step| in energy [eV]",
           title=f"Newton convergence (E = {e:.6f} eV)")
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("steps:", [f"{s:.1e}" for s in steps])


if __name__ == "__main__":
    main()
