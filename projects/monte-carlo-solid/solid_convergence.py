"""Accuracy check: Monte-Carlo error falls as 1/sqrt(N).

The solid's volume is known exactly (pi^2 / 2), so we measure the RMS error of
the estimate over many independent runs for increasing sample counts N.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import solid_mass

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

EXACT_VOLUME = np.pi**2 / 2
SAMPLE_COUNTS = [10**3, 10**4, 10**5, 10**6]


def rms_errors(sample_counts=SAMPLE_COUNTS, runs=60):
    """RMS error of the volume estimate at each N, over `runs` seeded repeats."""
    out = []
    for n in sample_counts:
        est = np.array([solid_mass.calculate_mass(n, seed=1000 * n + s)[0] for s in range(runs)])
        out.append(np.sqrt(np.mean((est - EXACT_VOLUME) ** 2)))
    return np.array(out)


def main():
    err = rms_errors()
    n = np.array(SAMPLE_COUNTS, dtype=float)
    slope = np.polyfit(np.log(n), np.log(err), 1)[0]
    p = EXACT_VOLUME / (4 * np.pi)                      # acceptance probability
    theory = 4 * np.pi * np.sqrt(p * (1 - p) / n)       # binomial standard error
    fig, ax = plt.subplots()
    ax.loglog(n, err, "o-", label=f"RMS error of volume (slope {slope:.2f})")
    ax.loglog(n, theory, "k--", label=r"theory $4\pi\sqrt{p(1-p)/N}$")
    ax.set(xlabel="samples N", ylabel="RMS error", title="Monte-Carlo error falls as 1/sqrt(N)")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"fitted slope {slope:.2f} (theory -0.5)")


if __name__ == "__main__":
    main()
