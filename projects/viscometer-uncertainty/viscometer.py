"""Monte-Carlo error propagation for the falling-sphere viscometer.

    eta = 2 (rho_s - rho_l) g r^2 t / (9 l)

Each measured quantity is drawn uniformly within its stated uncertainty.
"""
import argparse
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)

N_REPEAT = 10_000
G = 9.81  # m/s^2


def f(rho_s, rho_l, r, t, l):
    """Dynamic viscosity [Pa s] from sphere/liquid densities, radius, fall time and distance."""
    return 2 * (rho_s - rho_l) * G * r**2 * t / (9 * l)


def viscosity(n=N_REPEAT, seed=None, liquid_volume_scale=1.0):
    """Simulate `n` experiments; save histogram; return (mean, standard deviation).

    The liquid volume (0.1 +/- 0.0005) is used in units of `liquid_volume_scale` m^3.
    The default 1.0 takes it as m^3, exactly as the original script did, which makes
    the liquid density only 1.4 kg/m^3.  Pass 1e-3 if the volume is really in litres.
    """
    rng = np.random.default_rng(seed)

    def measure(value, uncertainty):
        return value + rng.uniform(-1, 1, n) * uncertainty

    diameter = measure(0.02, 0.0005)          # m
    r = diameter / 2
    rho_s = measure(0.021, 0.0005) / (4 / 3 * np.pi * r**3)   # sphere mass / volume
    rho_l = measure(0.14, 0.0005) / (measure(0.1, 0.0005) * liquid_volume_scale)  # liquid mass / volume
    l = measure(0.1, 0.0005)                  # fall distance, m
    t = measure(3.2, 0.05)                    # fall time, s

    values = f(rho_s, rho_l, r, t, l)

    fig, ax = plt.subplots()
    ax.hist(values, 100)
    ax.set(title="Distribution of random error", xlabel="viscosity [Pa s]",
           ylabel="Number of experiments")
    fig.savefig(OUTPUT / "viscosity_distribution.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

    return values.mean(), values.std()


def main(argv=None):
    p = argparse.ArgumentParser(description="Monte-Carlo uncertainty of a falling-sphere viscosity measurement.")
    p.add_argument("--n", type=int, default=N_REPEAT, help="number of simulated experiments (default 10000)")
    p.add_argument("--seed", type=int, default=None, help="random seed for reproducible runs")
    p.add_argument("--liquid-volume-in-litres", action="store_true",
                   help="treat the 0.1 liquid volume as litres instead of m^3 (default: m^3)")
    args = p.parse_args(argv)
    t0 = time.time()
    scale = 1e-3 if args.liquid_volume_in_litres else 1.0
    mean, sd = viscosity(args.n, args.seed, scale)
    print(f"viscosity = {mean:.4f} ± {sd:.4f} Pa s  ({time.time() - t0:.2f} s)")


if __name__ == "__main__":
    main()
