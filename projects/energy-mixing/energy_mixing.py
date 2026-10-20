"""Energy mixing among 2000 objects: 1000 at energy 0.1 and 1000 at 1000.
In each collision two random distinct objects share their combined energy equally.
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

N_OBJECTS = 2000
N_COLLISIONS = 10_000


def initial_energies():
    return np.concatenate([np.full(N_OBJECTS // 2, 0.1), np.full(N_OBJECTS // 2, 1000.0)])


def _collide(energy, rng):
    """One random pairwise collision, in place."""
    i, j = rng.choice(N_OBJECTS, size=2, replace=False)
    energy[i] = energy[j] = (energy[i] + energy[j]) / 2


def mixing(n_collisions=N_COLLISIONS, seed=None):
    """Run the simulation, save the energy histogram, return (mean, std) of the energies.

    The energies array is created inside the function, so repeated calls start
    fresh (the original used a global array that carried state between calls).
    """
    rng = np.random.default_rng(seed)
    energy = initial_energies()
    for _ in range(n_collisions):
        _collide(energy, rng)

    fig, ax = plt.subplots()
    ax.hist(energy, 100, color="blue")
    ax.set(xlabel="Energy", ylabel="Number of objects")
    fig.savefig(OUTPUT / "energy_mixing.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return energy.mean(), energy.std()


def variance_history(n_collisions=N_COLLISIONS, every=250, seed=None):
    """Return (collision counts, variance of the energies) sampled every `every` collisions."""
    rng = np.random.default_rng(seed)
    energy = initial_energies()
    counts, variances = [0], [energy.var()]
    for n in range(1, n_collisions + 1):
        _collide(energy, rng)
        if n % every == 0:
            counts.append(n)
            variances.append(energy.var())
    return np.array(counts), np.array(variances)


def predicted_variance(counts, initial_variance):
    """Expected variance after n random pair-averaging collisions: V0 (1 - 1/(N-1))^n."""
    return initial_variance * (1 - 1 / (N_OBJECTS - 1)) ** np.asarray(counts)


def main(argv=None):
    p = argparse.ArgumentParser(description="Random pairwise energy sharing among 2000 objects.")
    p.add_argument("--collisions", type=int, default=N_COLLISIONS, help="number of collisions (default 10000)")
    p.add_argument("--seed", type=int, default=None, help="random seed for reproducible runs")
    args = p.parse_args(argv)
    mean, sd = mixing(args.collisions, args.seed)
    print(f"mean energy = {mean:.4f} (conserved: {initial_energies().mean():.2f}), std = {sd:.2f}")


if __name__ == "__main__":
    main()
