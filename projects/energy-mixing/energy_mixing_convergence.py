"""Accuracy check: how fast does the energy distribution narrow?

Each collision removes (e_i - e_j)^2 / 2 from the sum of squares, so on average
the variance obeys  Var(n) = Var(0) * (1 - 1/(N-1))^n  ~  Var(0) exp(-n/(N-1)).
We average many simulated runs and compare with this prediction.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import energy_mixing as em

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)


def mean_variance_history(runs=30, n_collisions=em.N_COLLISIONS, every=250):
    """Average variance-vs-collisions curve over `runs` independent simulations."""
    histories = [em.variance_history(n_collisions, every, seed=s) for s in range(runs)]
    return histories[0][0], np.mean([h[1] for h in histories], axis=0)


def main():
    counts, var = mean_variance_history()
    theory = em.predicted_variance(counts, var[0])
    fig, ax = plt.subplots()
    ax.semilogy(counts, var, "o", ms=3, label="simulation (mean of 30 runs)")
    ax.semilogy(counts, theory, "k--", label=r"theory $V_0 (1-\frac{1}{N-1})^n$")
    ax.set(xlabel="collisions", ylabel="variance of energies",
           title="Energy mixing: exponential narrowing")
    ax.legend()
    fig.savefig(OUTPUT / "convergence.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"final variance: simulation {var[-1]:.1f}, theory {theory[-1]:.1f}")


if __name__ == "__main__":
    main()
