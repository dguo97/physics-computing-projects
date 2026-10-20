# Point-Mass Gravity

Net gravitational force on a test mass from a cloud of point masses read from a CSV.

A plain-language summary is in [ABOUT.txt](ABOUT.txt).

## Physics
**Governing law.** Newtonian gravity is linear, so the forces from the individual masses simply add (superposition):

$$\mathbf{F}=\sum_i G\,\frac{m_i M\,(\mathbf{r}_i-\mathbf{r})}{|\mathbf{r}_i-\mathbf{r}|^{3}}$$

**Why this method.** No approximation is needed for a finite set of point masses, so we evaluate the sum exactly, vectorised over all masses at once.

**What the result means.** The force vector points towards the mass-weighted pull of the cloud. Its magnitude (about $10^{-11}$ N for a 1 kg test mass here) is tiny because $G$ is small and the masses are small. Two equal masses cancel at their midpoint, and a single mass reproduces the inverse-square law (both tested).

## Run
```bash
python gravity_field.py
python gravity_field.py masses.csv --object 0 0 5 2    # test mass 2 kg at (0, 0, 5)
```
Reads its input from `data/`. Run any script with `--help` to list its options.

## Result
Force on a 1 kg mass at (1, 1, 1): **(−4.4e-11, 2.5e-12, −7.6e-12) N**.

## Tested by
`test_two_equal_masses_cancel_at_midpoint, test_single_mass_obeys_inverse_square_law` in [`tests/test_projects.py`](../../tests/test_projects.py).
