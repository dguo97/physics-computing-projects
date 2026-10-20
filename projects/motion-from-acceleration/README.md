# Motion from Acceleration

Position and velocity, with numerical error bounds, from a piecewise acceleration profile (speed-up, coast, brake).

A plain-language summary is in [ABOUT.txt](ABOUT.txt).

## Physics
**Governing relations.**

$$v(t)=\int_0^t a(t')\,dt',\qquad x(t)=\int_0^t v(t')\,dt'$$

The acceleration is a smooth logistic switch-on (until $t=10$ s), zero while coasting ($10$–$40$ s), then a logistic switch-on of braking centred at 47 s.

**Why this method.** The logistic profile has no convenient closed-form double integral, so we use adaptive quadrature (`scipy.integrate.quad`), which returns an error estimate with every integral. The jumps in the formula at $t=10$ and $t=40$ are passed to `quad` as break points so it does not waste evaluations (or lose accuracy) across them. The inner integral's error is propagated into the position error.

**What the result means.** After 50 s the object has covered about 894 m and is still moving at 12.1 m/s, because braking only begins around $t\approx47$ s. The numerical error is about $10^{-5}$ m, negligible next to any physical uncertainty in the model.

## Run
From the repository root, `cd projects/motion-from-acceleration` first, then:
```bash
python motion.py              # state at t = 50 s
python motion.py --time 30    # state at t = 30 s
```
Run any script with `--help` to list its options.

## Result
x(50 s) = **894.41 m**, v(50 s) = **12.114 m/s**.

## Tested by
`test_motion_is_consistent_with_direct_integration, test_velocity_constant_while_coasting` in [`tests/test_projects.py`](../../tests/test_projects.py).
