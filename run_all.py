"""Run every mini project and print a one-line result for each.
Figures are written to each project's own output/ folder."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for project in sorted((ROOT / "projects").iterdir()):
    if project.is_dir():
        sys.path.insert(0, str(project))

import route_analysis, g_estimator, motion, linear_solver, heating_model  # noqa: E402
import gravity_field, sheet_gravity, three_body, solid_mass, viscometer     # noqa: E402
import energy_mixing, well_energy, schrodinger                              # noqa: E402
import heating_convergence, three_body_convergence, schrodinger_convergence  # noqa: E402
import solid_convergence, viscometer_convergence, energy_mixing_convergence, well_convergence  # noqa: E402)

print("GPS route distance [m]   :", route_analysis.walk())
print("g, error [m/s^2]         :", g_estimator.gravitational_acceleration())
print("(x, err), (v, err)       :", motion.motion())
print("4x4 system               :", linear_solver.equ_system())
print("circuit                  :", linear_solver.circuit())
print("heating peak, time       :", heating_model.temperature(5, 0.5, 5))
print("force on point mass [N]  :", gravity_field.force(gravity_field.DATA / "point_masses.csv", [1, 1, 1, 1]))
print("elliptical sheet [m, F]  :", sheet_gravity.force([1, 1, 1]))
three_body.satellite(); print("three-body plots         : saved")
print("solid volume, mass       :", solid_mass.calculate_mass())
print("viscosity mean, sd       :", viscometer.viscosity())
print("energy mixing mean, sd   :", energy_mixing.mixing())
print("finite-well energy [eV]  :", well_energy.energy(10, 10e5))
print("Schrodinger energy [eV]  :", schrodinger.schroedinger())

print("\n-- accuracy studies (figures saved as convergence.png in each project) --")
for module in (heating_convergence, three_body_convergence, schrodinger_convergence,
               solid_convergence, viscometer_convergence, energy_mixing_convergence,
               well_convergence):
    module.main()
