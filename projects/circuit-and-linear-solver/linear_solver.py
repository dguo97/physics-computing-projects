"""Solve linear systems: a generic 4x4 example and a resistor-circuit problem."""
import argparse
import numpy as np


def equ_system():
    """Solve a 4x4 system A x = b."""
    a = np.array([[1, 2, 3, 4],
                  [2, 3, 4, 7],
                  [3, 4, 6, 9],
                  [4, 5, 7, 8]], dtype=float)
    b = np.array([1, 8, 3, 5], dtype=float)
    return np.linalg.solve(a, b)


def circuit():
    """Unknown currents and voltage in the given circuit (9 unknowns, 9 equations:
    Kirchhoff current/voltage laws)."""
    a = np.array([[0, -1, 0, 1, 0, 0, 0, 0, 0],
                  [0, 0, 0, 1, -1, 1, 0, 0, 0],
                  [0, 0, 1, 0, 0, 0, 1, -1, 0],
                  [0, 1, 0, 0, -1, 0, 0, 0, 1],
                  [0, 0, 0, 0, 0, -1, 0, 1, 0],
                  [0, 1, 0, 0, 5, 0, 0, 0, 0],
                  [-1, 0, 0, 0, 0, 8, 0, 0, 0],
                  [1, -1, 10, 0, 0, 0, 0, 0, 0],
                  [0, 0, -10, 0, 0, 0, 5, 0, 0]], dtype=float)
    b = np.array([0.1, 0, 0, 0, 0.1, -12, 0.1, 11.9, 5])
    return np.linalg.solve(a, b)


def circuit_loop_residuals(x):
    """Outer-loop equations that are NOT in the matrix; all should be ~0."""
    return [x[2] + x[6] - x[8],
            x[0] - x[1] + 5 * x[6] - 16.9,
            x[0] - x[1] - 5 * x[4] - 8 * x[5] - 11.9]


def main(argv=None):
    p = argparse.ArgumentParser(description="Solve the example linear system and/or the resistor circuit.")
    p.add_argument("--which", choices=["system", "circuit", "both"], default="both")
    args = p.parse_args(argv)
    if args.which in ("system", "both"):
        print("4x4 system:", equ_system())
    if args.which in ("circuit", "both"):
        sol = circuit()
        print("circuit   :", sol)
        print("loop-check residuals:", np.round(circuit_loop_residuals(sol), 12))


if __name__ == "__main__":
    main()
