"""Velocity and position (with numerical error estimates) from a given
acceleration profile, using nested numerical integration.

a(t): smooth speed-up until t=10, coasting until t=40, then braking.
"""
import argparse
import numpy as np
from scipy import integrate

T_END = 50.0
BREAKPOINTS = [10, 40]  # where the acceleration formula changes


def acceleration(time):
    """Acceleration [m/s^2] at `time` [s]."""
    if time < 10:
        return 2.5 / (1 + np.exp((time - 8) / 0.8))
    if time <= 40:  # original used `10 < time < 40`, which skipped t == 10 and t == 40
        return 0.0
    return -2.5 * (1 - 1 / (1 + np.exp((time - 47) / 1.2)))


def velocity(time):
    """Return (velocity, integration error estimate) at `time`."""
    pts = [b for b in BREAKPOINTS if b < time] or None
    return integrate.quad(acceleration, 0, time, points=pts)


def velocity_value(time):
    return velocity(time)[0]


def velocity_error(time):
    return velocity(time)[1]


def motion(time=T_END):
    """Return ((x, x_err), (v, v_err)) at `time`.

    x_err = outer-quad error + integral of the inner velocity errors,
    i.e. a +/- bound (the original script returned the full width, 2x this).
    """
    pts = [b for b in BREAKPOINTS if b < time] or None
    x, x_quad_err = integrate.quad(velocity_value, 0, time, points=pts)
    x_prop_err, _ = integrate.quad(velocity_error, 0, time, points=pts)
    return (x, x_quad_err + x_prop_err), velocity(time)


def main(argv=None):
    p = argparse.ArgumentParser(description="Position and velocity from a piecewise acceleration profile.")
    p.add_argument("--time", type=float, default=T_END, help="end time in seconds (default 50)")
    args = p.parse_args(argv)
    (x, x_err), (v, v_err) = motion(args.time)
    print(f"x({args.time:g}) = {x:.6f} ± {x_err:.2e} m")
    print(f"v({args.time:g}) = {v:.6f} ± {v_err:.2e} m/s")


if __name__ == "__main__":
    main()
