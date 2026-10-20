"""Monte-Carlo estimate of the volume and mass of the solid
    x^2 + y^2 < cos^2(z),  -pi/2 <= z <= pi/2,
with density rho = 10 |sin z| exp(-0.2|x|) (1 - y^2).
"""
import argparse
import numpy as np

N_SAMPLES = 1_000_000


def rho(x, y, z):
    """Density at (x, y, z)."""
    return 10 * np.abs(np.sin(z)) * np.exp(-0.2 * np.abs(x)) * (1 - y**2)


def in_surface(x, y, z):
    """True where the point lies inside the solid."""
    return x**2 + y**2 < np.cos(z) ** 2


def calculate_mass(n=N_SAMPLES, seed=None):
    """Return (volume, mass) from `n` uniform samples in the bounding box
    x, y in [-1, 1], z in [-pi/2, pi/2]  (box volume 4*pi)."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1, 1, n)
    y = rng.uniform(-1, 1, n)
    z = rng.uniform(-np.pi / 2, np.pi / 2, n)

    box_volume = 2 * 2 * np.pi
    inside = in_surface(x, y, z)
    volume = box_volume * inside.mean()
    mass = box_volume * np.mean(np.where(inside, rho(x, y, z), 0.0))
    return volume, mass


def main(argv=None):
    p = argparse.ArgumentParser(description="Monte-Carlo volume and mass of a solid with varying density.")
    p.add_argument("--n", type=int, default=N_SAMPLES, help="number of samples (default 1000000)")
    p.add_argument("--seed", type=int, default=None, help="random seed for reproducible runs")
    args = p.parse_args(argv)
    v, m = calculate_mass(args.n, args.seed)
    print(f"volume ≈ {v:.4f} (analytic {np.pi**2 / 2:.4f})   mass ≈ {m:.4f}")


if __name__ == "__main__":
    main()
