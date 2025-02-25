
# Shortcut

 - Run `./fit.sh` to fit the excitations in the given range.
 - Run `./decouple_<format>.sh` to unpack and decouple the Fourier-transformed densities (given in `readmeFiles/`, either as `cube` or `BTcompact` format) to gain the true transition densities.

# Fitting

- Get the initial input file: `cp readmeFiles/eval.yaml.0 ./eval.yaml`
- Fit, e.g., doing four iterations by calling `./BTevalSpec.py --niter=4 fit`
  The latter will produce the following add-line objectives (one after the other)
  ![addLineObj0](readmeFiles/addLineObj0.png   "Add-line objective (for fit iteration 1, initial)"))]
  ![addLineObj1](readmeFiles/addLineObj1.png   "Add-line objective (for fit iteration 2)"))]
  ![addLineObj2](readmeFiles/addLineObj2.png   "Add-line objective (for fit iteration 3)"))]
  ![addLineObj3](readmeFiles/addLineObj3.png   "Add-line objective (for fit iteration 4)"))]
  ![addLineObj4](readmeFiles/addLineObj4.png   "Add-line objective (for fit iteration 5, with error suppression)"))]
- There is a small line left at about `0.13 Ry`, which is smaller than two shape-error lines at about `0.20` Ry.
  Fit this line by suppressing the shape-error lines using `./BTevalSpec.py --niter=1 --wref=0.01 --signif --skip fit`.
  The latter call supresses the inital fit (`--skip`) and computes significances (`--signif`).
- The results are given in `readmeFiles/eval.yaml.final.fit` and `readmeFiles/excit_1.dat.final.fit`.
  The latter shows that all excitations have overall high significances.

# Transition-density decoupling

Computing the Fourier transform of the density `FT[n](r,omega)` at the excitation energies yields transition denties that are perturbed by the spectral overlap from transition densities of adjacent excitations.
Moreover, they are not properly normalized.

- To get the unperturbed transition densities `n_j`, you need the Fourier-transformed density at (or close to) the excitation energies.
  The latter are given either as Gaussian cube or BTcompact files in `readmeFiles/densft_<type>.tar.gz`.
  Expand the latter using `tar -xzf readmeFiles/densft_<type>.tar.gz`.
- Add the `DENSFT` section to `eval.yaml` (cf. `readmeFiles/eval.yaml.densft_<type>`).
  The number of density-energy pairs must equal the number of excitations.
  The given energies must be those at which the Fourier-transformed densities are evaluated (which may be different from but should be close to the actual excitation energies).
  Then run `./BTevalSpec.py decouple`, which finally returns the transition densities and the output below.
- Information about the transition dipoles that are evaluated from the transition densities via `mu = -e*int r*n_j d^3r` is added to `eval.yaml` and compared to the fitted ones as a consistency check.

```
Transition density | abs. norm | real norm | imag norm | sin^2(dAng) | dAbs
S1                      0.0259     99.92 %      0.08 %      0.00        0.90 %
S2                      0.0317    100.00 %      0.00 %      0.00        0.00 %
S3                      0.0255    100.00 %      0.00 %      0.00       -0.13 %
S4                      0.0299    100.00 %      0.00 %      0.00        0.00 %
S5                      0.0315    100.00 %      0.00 %      0.00        0.01 %
S6                      0.0200     99.95 %      0.05 %      0.00        0.41 %
S7                      0.0303    100.00 %      0.00 %      0.00       -0.03 %
S8                      0.0384    100.00 %      0.00 %      0.00       -0.48 %
S9                      0.0307     99.75 %      0.25 %      0.00        2.00 %
S10                     0.0292     99.90 %      0.10 %      0.00       -0.06 %
```

- Col 1) Excitation name
- Col 2) int |n_j| d^3r
- Col 3) Real-valued fraction (must be close to 100%)
- Col 4) Imag-valued fraction (must be close to   0%)
- Col 5) Compares the angle between the transition dipole from the fit and the one evaluated from the transition density (must be close to 0)
- Col 6) Compares the modulus of the transition dipole from the fit and the one evaluated from the transition density (must be close to 0)

If any of the columns 4-6 shows unrealistic values (significantly different from 0), something went wrong.
