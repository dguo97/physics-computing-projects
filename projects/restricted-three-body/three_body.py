"""Restricted three-body problem: a light satellite in the field of two bodies
(Earth-Moon-like, mu = 0.012277471), solved with 4th-order Runge-Kutta in the
rotating frame.  Initial conditions give the Arenstorf periodic orbit.
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

MU = 0.012277471
T_END = 18.0
STEPS = 100_000
INITIAL_STATE = (0.994, 0.0, 0.0, -2.0015851)   # x, vx, y, vy
PERIOD = 17.0652165601579625588917206249        # known period of this orbit


def equations(state):
    """d/dt of [x, vx, y, vy] (the system is autonomous, so no t argument)."""
    x, vx, y, vy = state
    r = np.sqrt((x + MU) ** 2 + y**2)
    s = np.sqrt((x - 1 + MU) ** 2 + y**2)
    return np.array([
        vx,
        x + 2 * vy - (1 - MU) * (x + MU) / r**3 - MU * (x - 1 + MU) / s**3,
        vy,
        y - 2 * vx - (1 - MU) * y / r**3 - MU * y / s**3,
    ])


def jacobi_constant(sol):
    """Conserved quantity of the rotating-frame problem; `sol` has columns x, vx, y, vy."""
    x, vx, y, vy = np.asarray(sol).T
    return (x**2 + y**2 + 2 * (1 - MU) / np.hypot(x + MU, y)
            + 2 * MU / np.hypot(x - 1 + MU, y) - vx**2 - vy**2)


def integrate(steps=STEPS, t_end=T_END):
    """RK4 from the initial state; returns (t, solution[steps, 4]). No plotting."""
    h = t_end / steps
    t = np.linspace(0.0, t_end, steps, endpoint=False)
    sol = np.empty((steps, 4))
    state = np.array(INITIAL_STATE)

    for i in range(steps):  # preallocated (np.vstack in a loop was O(n^2))
        sol[i] = state
        k1 = h * equations(state)
        k2 = h * equations(state + 0.5 * k1)
        k3 = h * equations(state + 0.5 * k2)
        k4 = h * equations(state + k3)
        state = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return t, sol


def satellite(steps=STEPS, t_end=T_END):
    """Integrate, save three plots, return (t, solution[steps, 4])."""
    t, sol = integrate(steps, t_end)

    fig, ax = plt.subplots()
    ax.plot(sol[:, 0], sol[:, 2])
    ax.set(title="Trajectory of the light object in the x-y plane",
           xlabel="x coordinate", ylabel="y coordinate")
    fig.savefig(OUTPUT / "three_body_trajectory.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(t, sol[:, 0], label="x")
    ax.plot(t, sol[:, 2], label="y")
    ax.set(title="x and y coordinates against time", xlabel="time", ylabel="coordinate")
    ax.legend()
    fig.savefig(OUTPUT / "three_body_coordinates.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.plot(sol[:, 1], sol[:, 3])
    ax.set(title="Velocity in y against velocity in x",
           xlabel="velocity in x", ylabel="velocity in y")
    fig.savefig(OUTPUT / "three_body_velocity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    return t, sol


def main(argv=None):
    p = argparse.ArgumentParser(description="Restricted three-body problem (Arenstorf orbit) with RK4.")
    p.add_argument("--steps", type=int, default=STEPS, help="RK4 steps (default 100000)")
    p.add_argument("--t-end", type=float, default=T_END, help="end time (default 18; orbit period is 17.065)")
    args = p.parse_args(argv)
    t0 = time.time()
    _, sol = satellite(args.steps, args.t_end)
    drift = np.ptp(jacobi_constant(sol))
    print(f"final state {sol[-1]}  Jacobi-constant drift {drift:.2e}  ({time.time() - t0:.1f} s)")


if __name__ == "__main__":
    main()
