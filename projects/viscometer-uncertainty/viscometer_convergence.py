"""Accuracy check: the precision of the *estimated mean* viscosity improves as 1/sqrt(N).

The spread of a single virtual experiment (about 0.77 Pa s) is fixed by the input
uncertainties; running more experiments only pins down the mean better.  We repeat
the N-experiment simulation many times and look at the scatter of the resulting means.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import viscometer

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

SAMPLE_COUNTS = [10**2, 10**3, 10**4, 10**5]


def scatter_of_means(sample_counts=SAMPLE_COUNTS, runs=60):
    """Return (scatter of the mean at each N, single-experiment spread)."""
    scatter = []
    for n in sample_counts:
        means = [viscometer.viscosity(n, seed=1000 * n + s)[0] for s in range(runs)]
        scatter.append(np.std(means))
    spread = viscometer.viscosity(10**6, seed=0)[1]
    return np.array(scatter), spread


def main():
    scatter, spread = scatter_of_means()
    n = np.array(SAMPLE_COUNTS, dtype=float)
    slope = np.polyfit(np.log(n), np.log(scatter), 1)[0]
    fig, ax = plt.subplots()
    ax.loglog(n, scatter, "o-", label=f"scatter of the mean (slope {slope:.2f})")
    ax.loglog(n, spread / np.sqrt(n), "k--", label=r"theory $\sigma/\sqrt{N}$")
    ax.set(xlabel="simulated experiments N", ylabel="std of the mean viscosity [Pa s]",
           title="Viscometer: precision of the mean vs N")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"fitted slope {slope:.2f} (theory -0.5), single-experiment spread {spread:.3f}")


if __name__ == "__main__":
    main()
