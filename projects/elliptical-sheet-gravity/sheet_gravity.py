"""Mass of a uniform elliptical sheet and the gravitational force it exerts on a
10 kg point mass.

Sheet: semi-axes a=10 m (x), b=5 m (y), lying in the z=0 plane, density 1000 kg/m^2.
"""
import argparse
import numpy as np
from scipy import integrate

G = 6.674e-11
A, B = 10.0, 5.0       # semi-axes [m]
SIGMA = 1000.0         # surface density [kg/m^2]
M_TEST = 10.0          # test mass [kg]


def _x_limit(y):
    """Half-width of the ellipse at height y."""
    return A * np.sqrt(max(0.0, 1 - y**2 / B**2))


def _integrate_over_sheet(func):
    # dblquad integrates func(inner, outer): inner = x, outer = y
    return integrate.dblquad(func, -B, B, lambda y: -_x_limit(y), _x_limit)[0]


def force(position):
    """Return [sheet mass, Fx, Fy, Fz] for a point mass at `position` = (x, y, z)."""
    px, py, pz = position

    def component(index):
        def integrand(x, y):
            d = np.array([x - px, y - py, -pz])
            return G * M_TEST * SIGMA * d[index] / np.linalg.norm(d) ** 3
        return integrand

    mass = _integrate_over_sheet(lambda x, y: SIGMA)
    return [mass] + [_integrate_over_sheet(component(i)) for i in range(3)]


def main(argv=None):
    p = argparse.ArgumentParser(description="Gravity of a uniform elliptical sheet on a 10 kg point mass.")
    p.add_argument("--position", nargs=3, type=float, metavar=("X", "Y", "Z"), default=[1, 1, 1],
                   help="position of the point mass in metres (default 1 1 1)")
    args = p.parse_args(argv)
    m, fx, fy, fz = force(args.position)
    print(f"mass = {m:.1f} kg (analytic {SIGMA * np.pi * A * B:.1f})")
    print(f"force = ({fx:.4e}, {fy:.4e}, {fz:.4e}) N")


if __name__ == "__main__":
    main()
