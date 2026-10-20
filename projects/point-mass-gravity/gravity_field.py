"""Gravitational force on a test mass from a set of point masses.

Data columns: x, y, z, mass.  `obj` = [x, y, z, mass] of the test object.
"""
import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

G = 6.674e-11  # m^3 kg^-1 s^-2


def force(filename, obj):
    """Return the force vector [Fx, Fy, Fz] on `obj` due to the system in `filename`."""
    data = np.loadtxt(filename, delimiter=",", ndmin=2)  # ndmin: also accept a single-row file
    positions, masses = data[:, :3], data[:, 3]

    separation = positions - np.asarray(obj[:3], dtype=float)  # object -> each mass
    distance = np.linalg.norm(separation, axis=1)
    return (G * obj[3] * masses[:, None] * separation / distance[:, None] ** 3).sum(axis=0)


def main(argv=None):
    p = argparse.ArgumentParser(description="Gravitational force on a test mass from a set of point masses.")
    p.add_argument("file", nargs="?", type=Path, default=DATA / "point_masses.csv",
                   help="CSV with columns x, y, z, mass (default: bundled sample)")
    p.add_argument("--object", nargs=4, type=float, metavar=("X", "Y", "Z", "MASS"),
                   default=[1, 1, 1, 1], help="test object position and mass (default 1 1 1 1)")
    args = p.parse_args(argv)
    f = force(args.file, args.object)
    print(f"force = ({f[0]:.4e}, {f[1]:.4e}, {f[2]:.4e}) N,  |F| = {np.linalg.norm(f):.4e} N")


if __name__ == "__main__":
    main()
