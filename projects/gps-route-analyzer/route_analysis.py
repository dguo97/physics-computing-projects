"""Plot a GPS journey and compute the distance travelled.

Data columns: x [m], y [m], elevation [m], timestamp.
Outputs: path in the x-y plane, elevation vs distance, total distance.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)


def walk(filename=DATA / "journey_gps.csv"):
    """Save the x-y path and elevation profile; return total x-y distance [m]."""
    data = np.loadtxt(filename, delimiter=",")
    x, y, elevation = data[:, 0], data[:, 1], data[:, 2]

    fig, ax = plt.subplots()
    ax.plot(x / 1000, y / 1000)  # metres -> km
    ax.set(xlabel="x coordinate [km]", ylabel="y coordinate [km]",
           title="Path in x-y plane")
    fig.savefig(OUTPUT / "journey_path_xy.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

    step = np.hypot(np.diff(x), np.diff(y))
    distance = np.concatenate(([0.0], np.cumsum(step)))

    fig, ax = plt.subplots()
    ax.plot(distance, elevation)
    ax.set(xlabel="Distance [m]", ylabel="Elevation [m]",
           title="Elevation against distance")
    fig.savefig(OUTPUT / "journey_elevation.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

    return distance[-1]


def main(argv=None):
    p = argparse.ArgumentParser(description="Plot a GPS route and compute the distance travelled.")
    p.add_argument("file", nargs="?", type=Path, default=DATA / "journey_gps.csv",
                   help="CSV with columns x[m], y[m], elevation[m], ... (default: bundled sample)")
    args = p.parse_args(argv)
    print(f"Total distance: {walk(args.file):.1f} m")


if __name__ == "__main__":
    main()
