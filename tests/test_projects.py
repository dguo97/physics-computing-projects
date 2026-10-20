"""Physics sanity checks: each test compares a project's output against an
analytic result, a conservation law, or an independent solver."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp, trapezoid

import route_analysis
import g_estimator
import motion
import linear_solver
import heating_model
import gravity_field
import sheet_gravity
import three_body
import solid_mass
import viscometer
import energy_mixing
import well_energy
import schrodinger
import importlib


MODULES = [route_analysis, g_estimator, heating_model, three_body,
           viscometer, energy_mixing, schrodinger]


@pytest.fixture(autouse=True)
def _no_figure_side_effects(tmp_path, monkeypatch):
    """Send any figures a test generates to a temp dir, not the project's output/."""
    for mod in MODULES:
        monkeypatch.setattr(mod, "OUTPUT", tmp_path)


# ---- data analysis ---------------------------------------------------------
def test_route_distance_matches_reference():
    assert route_analysis.walk() == pytest.approx(14348.68, abs=0.01)


def test_g_estimator_recovers_known_g(tmp_path):
    rng = np.random.default_rng(1)
    h = np.linspace(1, 20, 40)
    sigma = np.full_like(h, 0.01)
    t = np.sqrt(2 * h / 9.81) + rng.normal(0, 0.01, h.size)
    csv = tmp_path / "drops.csv"
    np.savetxt(csv, np.column_stack([h, t, sigma]), delimiter=",")
    g, err = g_estimator.gravitational_acceleration(csv)
    assert abs(g - 9.81) < 4 * err
    assert 0 < err < 0.2


# ---- numerical methods -----------------------------------------------------
def test_motion_is_consistent_with_direct_integration():
    (x, x_err), (v, v_err) = motion.motion(50)
    t = np.linspace(0, 50, 200_001)
    a = np.array([motion.acceleration(ti) for ti in t[::50]])
    ts = t[::50]
    vel = np.concatenate(([0], np.cumsum((a[1:] + a[:-1]) / 2 * np.diff(ts))))
    assert v == pytest.approx(vel[-1], rel=1e-3)
    assert x == pytest.approx(trapezoid(vel, ts), rel=1e-3)
    assert x_err < 1e-3 and v_err < 1e-3


def test_velocity_constant_while_coasting():
    assert motion.velocity(15)[0] == pytest.approx(motion.velocity(35)[0], abs=1e-9)


def test_linear_systems_solve_their_equations():
    x = linear_solver.equ_system()
    a = np.array([[1, 2, 3, 4], [2, 3, 4, 7], [3, 4, 6, 9], [4, 5, 7, 8]])
    assert a @ x == pytest.approx([1, 8, 3, 5])
    assert np.abs(linear_solver.circuit_loop_residuals(linear_solver.circuit())).max() < 1e-9


def test_rk4_matches_scipy_solver():
    ymax, xmax = heating_model.temperature(5, 0.5, 5)
    ref = solve_ivp(lambda x, y: heating_model.ff(x, y[0], 5, 0.5), (0, 5), [0.0],
                    rtol=1e-10, atol=1e-12, dense_output=True)
    xs = np.linspace(0, 5, 100_001)
    ys = ref.sol(xs)[0]
    assert ymax == pytest.approx(ys.max(), rel=1e-6)
    assert xmax == pytest.approx(xs[ys.argmax()], abs=1e-3)


# ---- gravity & orbits ------------------------------------------------------
def test_two_equal_masses_cancel_at_midpoint(tmp_path):
    csv = tmp_path / "masses.csv"
    np.savetxt(csv, [[-1, 0, 0, 5.0], [1, 0, 0, 5.0]], delimiter=",")
    assert gravity_field.force(csv, [0, 0, 0, 1.0]) == pytest.approx([0, 0, 0], abs=1e-25)


def test_single_mass_obeys_inverse_square_law(tmp_path):
    csv = tmp_path / "mass.csv"
    np.savetxt(csv, [[3.0, 4.0, 0.0, 1e6]], delimiter=",")
    f = gravity_field.force(csv, [0, 0, 0, 2.0])
    expected = gravity_field.G * 1e6 * 2.0 / 25.0
    assert np.linalg.norm(f) == pytest.approx(expected, rel=1e-12)
    assert f[0] / f[1] == pytest.approx(3 / 4)  # pulled toward the mass


def test_sheet_mass_equals_sigma_pi_a_b():
    mass = sheet_gravity.force([0, 0, 1])[0]
    assert mass == pytest.approx(sheet_gravity.SIGMA * np.pi * sheet_gravity.A * sheet_gravity.B, rel=1e-6)


def test_sheet_force_symmetry_and_far_field():
    _, fx, fy, fz = sheet_gravity.force([0, 0, 1])
    assert fx == pytest.approx(0, abs=1e-15) and fy == pytest.approx(0, abs=1e-15)
    assert fz < 0  # attracted down toward the sheet
    mass, _, _, fz_far = sheet_gravity.force([0, 0, 1000])
    point_mass = sheet_gravity.G * sheet_gravity.M_TEST * mass / 1000**2
    assert fz_far == pytest.approx(-point_mass, rel=1e-3)


@pytest.fixture(scope="module")
def orbit(tmp_path_factory):
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(three_body, "OUTPUT", tmp_path_factory.mktemp("figs"))
        return three_body.satellite()


def test_three_body_orbit_closes_after_one_period(orbit):
    t, sol = orbit
    period = 17.0652165601579625588917206249
    i = int(round(period / (t[1] - t[0])))
    assert sol[i] == pytest.approx(sol[0], abs=0.02)


def test_three_body_conserves_jacobi_constant(orbit):
    _, sol = orbit
    assert np.ptp(three_body.jacobi_constant(sol)) < 1e-5


# ---- Monte Carlo -----------------------------------------------------------
def test_solid_volume_near_pi_squared_over_two():
    volume, mass = solid_mass.calculate_mass(seed=42)
    assert volume == pytest.approx(np.pi**2 / 2, abs=0.03)
    assert mass > 0


def test_viscosity_mean_close_to_nominal_value():
    r = 0.01
    rho_s = 0.021 / (4 / 3 * np.pi * r**3)
    nominal = viscometer.f(rho_s, 0.14 / 0.1, r, 3.2, 0.1)
    mean, sd = viscometer.viscosity(seed=7)
    assert mean == pytest.approx(nominal, rel=0.01)
    assert 0 < sd / mean < 0.05


def test_energy_mixing_conserves_energy_and_narrows_distribution():
    mean, sd = energy_mixing.mixing(seed=3)
    start = energy_mixing.initial_energies()
    assert mean == pytest.approx(start.mean())
    assert sd < start.std()


def test_energy_mixing_runs_are_independent():
    assert energy_mixing.mixing(seed=5) == energy_mixing.mixing(seed=5)


# ---- quantum ---------------------------------------------------------------
def test_well_energy_is_a_root_of_the_equation():
    v0, width = 10, 10e5
    e = well_energy.energy(v0, width)
    k = np.sqrt(2 * well_energy.MC2 * e) / well_energy.HC
    kp = np.sqrt(2 * well_energy.MC2 * (e - v0)) / well_energy.HC
    assert abs(k + np.tan(k * width) * kp) < 1e-9


def test_schrodinger_energy_and_normalisation():
    e = schrodinger.find_energy()
    assert e == pytest.approx(69.19849364, abs=1e-7)
    x, psi = schrodinger.solve(e)
    assert abs(psi[-1]) < 1e-6 * np.abs(psi).max()  # boundary condition at 2a


# ---- convergence studies (accuracy claims in the READMEs) -------------------
def test_rk4_heating_is_fourth_order():
    conv = importlib.import_module("heating_convergence")
    _, _, slope = conv.convergence_order()
    assert 3.7 < slope < 4.3


def test_three_body_jacobi_drift_shrinks_like_h4():
    conv = importlib.import_module("three_body_convergence")
    _, drift, slope = conv.convergence_order([40_000, 80_000, 160_000])
    assert np.all(np.diff(drift) < 0)
    assert 3.8 < slope < 5.0


def test_schrodinger_energy_converges_with_fourth_order():
    conv = importlib.import_module("schrodinger_convergence")
    _, diff, slope = conv.convergence_order([250, 500, 1000])
    assert 3.5 < slope < 4.5


def test_monte_carlo_error_falls_like_inverse_sqrt_n():
    conv = importlib.import_module("solid_convergence")
    err = conv.rms_errors([10**3, 10**4, 10**5], runs=40)
    slope = np.polyfit(np.log([1e3, 1e4, 1e5]), np.log(err), 1)[0]
    assert -0.65 < slope < -0.35


def test_viscometer_mean_precision_improves_like_inverse_sqrt_n():
    conv = importlib.import_module("viscometer_convergence")
    scatter, spread = conv.scatter_of_means([10**2, 10**3, 10**4], runs=40)
    slope = np.polyfit(np.log([1e2, 1e3, 1e4]), np.log(scatter), 1)[0]
    assert -0.65 < slope < -0.35
    assert scatter[-1] == pytest.approx(spread / 100, rel=0.35)


def test_energy_variance_decays_as_predicted():
    conv = importlib.import_module("energy_mixing_convergence")
    counts, var = conv.mean_variance_history(runs=20, n_collisions=6000, every=500)
    theory = energy_mixing.predicted_variance(counts, var[0])
    assert var[-1] == pytest.approx(theory[-1], rel=0.1)


def test_newton_converges_quadratically():
    conv = importlib.import_module("well_convergence")
    _, steps = conv.newton_steps()
    assert steps[-1] < 1e-10 and len(steps) < 15
    tail = steps[-4:]                       # last steps roughly square each time
    assert tail[1] < 10 * tail[0] ** 2 and tail[2] < 10 * tail[1] ** 2


# ---- command-line interfaces -------------------------------------------------
@pytest.mark.parametrize("module, args, expected", [
    (viscometer, ["--n", "2000", "--seed", "1"], "viscosity"),
    (solid_mass, ["--n", "20000", "--seed", "1"], "volume"),
    (energy_mixing, ["--collisions", "500", "--seed", "1"], "mean energy"),
    (well_energy, ["--v0", "10"], "E ="),
    (motion, ["--time", "20"], "x(20)"),
    (sheet_gravity, ["--position", "0", "0", "2"], "mass"),
    (gravity_field, ["--object", "0", "0", "0", "2"], "force"),
    (linear_solver, ["--which", "system"], "4x4"),
    (heating_model, ["--steps", "2000"], "max temperature"),
])
def test_cli_runs(module, args, expected, capsys):
    module.main(args)
    assert expected in capsys.readouterr().out
