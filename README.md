# Computational Physics Mini-Projects

Thirteen small, self-contained Python projects in numerical physics: data fitting,
numerical integration, ODE solvers, Monte-Carlo simulation, orbital mechanics and
quantum mechanics. Each project has its own folder with a README (governing equation,
why the method fits the problem, what the result means), a command-line interface, its
code, and where relevant its data, figures and an accuracy study.

These are a couple of projects I made during university, later cleaned up, fixed and tested.

**Author:** David Guo  ·  **License:** [MIT](LICENSE)

![tests](../../actions/workflows/tests.yml/badge.svg)

## Projects

| Project | What it does |
|---|---|
| [GPS Route Analyzer](projects/gps-route-analyzer) | Turns raw GPS samples (x, y, elevation) into a route map, an elevation profile and a total distance. |
| [Free-Fall g Estimator](projects/free-fall-g-estimator) | Estimates the gravitational acceleration from drop experiments with per-point timing uncertainties. |
| [Motion from Acceleration](projects/motion-from-acceleration) | Position and velocity, with numerical error bounds, from a piecewise acceleration profile (speed-up, coast, brake). |
| [Circuit & Linear-System Solver](projects/circuit-and-linear-solver) | Solves a general 4×4 system and a 9-unknown resistor circuit written as Kirchhoff current/voltage equations. |
| [RK4 Heating Model](projects/rk4-heating-model) | Integrates dy/dx = a·e^(−2x) − b(y⁴ − 1): a body heated by a decaying source while radiating to its surroundings. |
| [Point-Mass Gravity](projects/point-mass-gravity) | Net gravitational force on a test mass from a cloud of point masses read from a CSV. |
| [Elliptical Sheet Gravity](projects/elliptical-sheet-gravity) | Mass of a uniform elliptical plate and the force it exerts on a 10 kg mass above it. |
| [Restricted Three-Body Problem](projects/restricted-three-body) | A light satellite moving in the field of two massive bodies (μ = 0.01228) in the rotating frame — the Arenstorf periodic orbit. |
| [Monte-Carlo Solid](projects/monte-carlo-solid) | Volume and mass of the solid x² + y² < cos²z with a spatially varying density. |
| [Viscometer Uncertainty](projects/viscometer-uncertainty) | Monte-Carlo error propagation for the falling-sphere method of measuring a liquid's viscosity. |
| [Energy-Mixing Simulation](projects/energy-mixing) | 2000 objects — half cold (0.1), half hot (1000) — repeatedly pair up at random and share their energy equally. |
| [Finite-Well Energy (Newton)](projects/finite-well-newton) | Finds an allowed energy by solving the transcendental condition k + k′·tan(kw) = 0. |
| [Schrödinger Shooting Method](projects/schrodinger-shooting) | Eigenenergy and normalised wave function for a particle in a box with a step in its floor. |

Each project folder also has an `ABOUT.txt`, a plain-language explanation of what it does.

## Quick start
```bash
pip install -r requirements-dev.txt
pytest                # 34 checks (~1 min)
python run_all.py     # run every project and every accuracy study
python projects/viscometer-uncertainty/viscometer.py --n 50000 --seed 1
python projects/schrodinger-shooting/schrodinger.py --v0 8 --e-low 60 --e-high 120
```
### Running one project
Either give the script's path from the repository root, or move into the project's folder first:
```bash
python projects/rk4-heating-model/heating_model.py       # from the repository root
cd projects/rk4-heating-model && python heating_model.py  # or from inside the folder
```
Each project's own README lists its commands and options.

### Run commands for every project
Run these from the repository root (add `--help` to any of them for options).

| Project | Run | Accuracy study |
|---|---|---|
| GPS Route Analyzer | `python projects/gps-route-analyzer/route_analysis.py` | – |
| Free-Fall g Estimator | `python projects/free-fall-g-estimator/g_estimator.py` | – |
| Motion from Acceleration | `python projects/motion-from-acceleration/motion.py` | – |
| Circuit & Linear-System Solver | `python projects/circuit-and-linear-solver/linear_solver.py` | – |
| RK4 Heating Model | `python projects/rk4-heating-model/heating_model.py` | `python projects/rk4-heating-model/heating_convergence.py` |
| Point-Mass Gravity | `python projects/point-mass-gravity/gravity_field.py` | – |
| Elliptical Sheet Gravity | `python projects/elliptical-sheet-gravity/sheet_gravity.py` | – |
| Restricted Three-Body Problem | `python projects/restricted-three-body/three_body.py` | `python projects/restricted-three-body/three_body_convergence.py` |
| Monte-Carlo Solid | `python projects/monte-carlo-solid/solid_mass.py` | `python projects/monte-carlo-solid/solid_convergence.py` |
| Viscometer Uncertainty | `python projects/viscometer-uncertainty/viscometer.py` | `python projects/viscometer-uncertainty/viscometer_convergence.py` |
| Energy-Mixing Simulation | `python projects/energy-mixing/energy_mixing.py` | `python projects/energy-mixing/energy_mixing_convergence.py` |
| Finite-Well Energy (Newton) | `python projects/finite-well-newton/well_energy.py` | `python projects/finite-well-newton/well_convergence.py` |
| Schrödinger Shooting Method | `python projects/schrodinger-shooting/schrodinger.py` | `python projects/schrodinger-shooting/schrodinger_convergence.py` |

Every project script has a command-line interface; run it with `--help` to see the options.

## Run on Google Colab
No local setup is needed. Colab already has NumPy, SciPy and Matplotlib, so you only need to fetch the code.

1. Open [colab.research.google.com](https://colab.research.google.com) and create a **New notebook**.
2. Paste this into the first cell and run it 
   ```python
   !git clone https://github.com/dguo97/physics-computing-projects.git
   %cd physics-computing-projects
   !pip install -q -r requirements.txt
   ```
3. Run any project with `!python`, using the same options as on the command line:
   ```python
   !python projects/viscometer-uncertainty/viscometer.py --n 50000 --seed 1
   !python projects/schrodinger-shooting/schrodinger.py
   ```
4. Scripts save their figures as PNG files instead of opening a window, so display one in the notebook like this:
   ```python
   from IPython.display import Image
   Image("projects/schrodinger-shooting/output/schrodinger_wavefunction.png")
   ```
   Or show every figure of a project at once:
   ```python
   import glob
   from IPython.display import display, Image
   for f in sorted(glob.glob("projects/restricted-three-body/output/*.png")):
       display(Image(f))
   ```
5. To run everything, or the tests:
   ```python
   !python run_all.py                       # all projects and accuracy studies (about a minute)
   !pip install -q -r requirements-dev.txt
   !pytest -q                               # 34 checks
   ```

Files in Colab are deleted when the session ends. To keep the figures, download them from the file browser on the left, or mount Google Drive.

## Accuracy studies
Seven projects include a short study showing the method behaves as theory predicts:

| Project | What is checked | Observed |
|---|---|---|
| RK4 Heating Model | error vs step size | order 3.94 (theory 4) |
| Restricted Three-Body | Jacobi-constant drift vs step size | order 4.3 (theory 4) |
| Schrödinger Shooting | eigenenergy change vs step size | order 4.00 (theory 4) |
| Monte-Carlo Solid | RMS error vs samples | slope −0.48 (theory −0.5) |
| Viscometer Uncertainty | precision of the mean vs experiments | slope −0.51 (theory −0.5) |
| Energy-Mixing | variance decay vs prediction | 1674 vs 1678 |
| Finite-Well Newton | Newton step sizes | quadratic convergence |

## Tests
`tests/test_projects.py` checks each project against something independent: an analytic
result (sheet mass = σπab, solid volume ≈ π²/2, inverse-square law), a conservation law
(Jacobi constant, total energy), a periodicity condition (the three-body orbit closes after
one period), another solver (`scipy.integrate.solve_ivp`), the convergence rates above, and
that every command-line interface runs.

## Layout
```
projects/<name>/   code, README.md, data/ (if any), output/ (figures)
tests/             pytest suite
run_all.py         runs every project and accuracy study
```
