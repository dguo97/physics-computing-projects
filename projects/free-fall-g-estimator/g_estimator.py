"""Estimate g from drop data (height, fall time, time uncertainty).

Model: s = g t^2 / 2  =>  t = sqrt(2 s / g), fitted by weighted least squares.
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUTPUT = HERE / "output"
OUTPUT.mkdir(exist_ok=True)


def fall_time(height, gacc):
    """Time to fall `height` from rest under acceleration `gacc`."""
    return np.sqrt(2 * height / gacc)


def gravitational_acceleration(filename=DATA / "free_fall_times.csv"):
    """Return (g, uncertainty on g) and save the data + fit plot."""
    data = np.loadtxt(filename, delimiter=",")
    height, time, sigma = data[:, 0], data[:, 1], data[:, 2]

    popt, pcov = curve_fit(fall_time, height, time, p0=9.81,
                           sigma=sigma, absolute_sigma=True)
    g, g_err = popt[0], np.sqrt(pcov[0, 0])

    fig, ax = plt.subplots()
    ax.errorbar(height, time, yerr=sigma, fmt="o", ms=3, label="data")
    h = np.linspace(0, height.max() * 1.05, 500)
    ax.plot(h, fall_time(h, g), "r-", label=f"fit: g = {g:.2f} ± {g_err:.2f} m/s²")
    ax.set(xlabel="Height [m]", ylabel="Time [s]", title="Fall time against height")
    ax.legend()
    fig.savefig(OUTPUT / "gravity_fit.png", dpi=150, bbox_inches='tight')
    plt.close(fig)

    return g, g_err


def main(argv=None):
    p = argparse.ArgumentParser(description="Estimate g from free-fall data.")
    p.add_argument("file", nargs="?", type=Path, default=DATA / "free_fall_times.csv",
                   help="CSV with columns height[m], time[s], time_uncertainty[s] (default: bundled sample)")
    args = p.parse_args(argv)
    g, err = gravitational_acceleration(args.file)
    print(f"g = {g:.4f} ± {err:.4f} m/s^2")


if __name__ == "__main__":
    main()
