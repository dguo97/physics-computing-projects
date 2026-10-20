# Elliptical Sheet Gravity

Mass of a uniform elliptical plate and the force it exerts on a 10 kg mass above it.

## Physics
**Governing integral.** Summing the pull of every area element of the sheet (surface density $\sigma$, semi-axes $a=10$ m, $b=5$ m):

$$\mathbf{F}=G\,m\,\sigma\iint_{\text{ellipse}}\frac{\mathbf{r}'-\mathbf{r}}{|\mathbf{r}'-\mathbf{r}|^{3}}\,dA,\qquad M=\sigma\iint dA=\sigma\pi ab$$

**Why this method.** For a point off the symmetry axis the force involves elliptic integrals, so we integrate numerically with nested adaptive quadrature (`scipy.integrate.dblquad`), integrating over $x$ first for each $y$ with the ellipse as the limits. The exactly known mass $\sigma\pi ab$ is the check that the integration region is right.

**What the result means.** The $z$ component is negative: the sheet pulls the mass down towards itself. On the central axis the $x$ and $y$ components vanish by symmetry, and far from the sheet the force approaches that of a point mass $M$ (both tested).

## Run
```bash
python sheet_gravity.py
python sheet_gravity.py --position 0 0 2
```
Run any script with `--help` to list its options.

## Result
Mass **157 079.6 kg** = σπab exactly.

## Tested by
`test_sheet_mass_equals_sigma_pi_a_b, test_sheet_force_symmetry_and_far_field` in [`tests/test_projects.py`](../../tests/test_projects.py).
