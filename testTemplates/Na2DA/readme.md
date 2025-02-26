# Description

The two dipole-moment files provided stem from a calculation of an imperfect H-aggregate consisting of two Na2 dimers, of which one is considered the donor and the other the acceptor.
the connection axis betweent the molecules is the `x` axis.
The donor (in the `x<0` halfspace) is aligned in `z`-direction.
The acceptor (in the `x>0` halfspace) is slightly rotated around the connectino axis and shows a slightly different bond legnth.
The two dipole moments have been taken by integrating the two halfspaces separatelly.

Since a single Na2 has one dominant excitation along its axis, we expect to see is the spectrum of an imperfect H-aggregate with two lines: One of which shows a symmetric coupling between the two Na2 excitations (at lower energy) and one an anti-symmetric coupling (at higher energy).
The symmetry can be derived from the signs that the excitations show in the different molecules' halfspaces.

The external laser excitation (electrid dipole field) is focused on the donor molecule so that (i) both excitations are excited by this pulse and (ii) are visible in the molecule-specific dipole moments.

In the current version of `BTevalSpec.py`, the computed oscillator strength is not the real one.

# Shortcut

Run `./run.sh`

# Fitting

- Get the initial input file: `cp readmeFiles/eval.yaml.0 ./eval.yaml`

## Iteration 1

One iteration

```bash
./BTevalSpec.py --niter=1 fit
```

produces the add-line objective

![addLineObj0](readmeFiles/addLineObj0.png   "Add-line objective (for fit iteration 1, initial)")

and two excitations:

```
cat readmeFiles/excit_1.dat.1
# name  | energy[Ry] |  strength  | phase[pi]  | ... signifFit|signifAng|signifExc|signifErr|signifRng|signifPha
S1            0.14912      0.01484      0.40642  ...   0.00      0.93      0.00      1.00      1.00      1.00
S2            0.15940      0.03733      0.74448        0.00      0.99      0.00      1.00      1.00      1.00
```

Inspecting the donor- and acceptor specific Pade spectra (z-component, real part), reveals that the two peaks show different signs on the two molecules:

![PadeDA](readmeFiles/padeByArea_z.1.png   "Pade (z,real) spectrum of donor and acceptor")

As expected from the H-aggregate, the symmetrically coupled excitation (equal signs) is at higher energy and the anti-symmetrically coupled one at lower energy.


After this fit, the comprehensive fit error is already as small as `0.02`, i.e., most of the spectrum is already well described by these two excitations.
One can still find more excitations, which may be interesting to show in the following.
Note, however, that these need not be physically correct.
The following just demonstrates the script functionality.

## Iterations 2 & 3

To separate real excitations from the line-shape errors from the already fitted, much larger exitations, you can use the error-suppression and call

```bash
./BTevalSpec.py --niter=1 --wref=1  --skip fit
./BTevalSpec.py --niter=1 --wref=10 --skip fit
```

The latter gives a comprehensive fit error of `0.015` and produces the add-line objectives

![addLineObj1_wref](readmeFiles/addLineObj1_wref.png   "Add-line objective (for fit iteration 2, wref=1)")
![addLineObj2_wref1ß](readmeFiles/addLineObj2_wref10.png   "Add-line objective (for fit iteration 3, wref=10)")

## Check significances & remove S6 & re-fit

Looking into `excit_1.dat` reveals that `S6` shows a very low significance (see `readmeFiles/excit_1.dat.3iter`, note that `signifFit` and `signifExc` are zero since the `--signif` flag was not given):

```text
# name  | energy[Ry] |  strength  | phase[pi]  | ... signifFit|signifAng|signifExc|signifErr|signifRng|signifPha
S1            0.14912      0.01485      0.40642  ...  0.00      0.93      0.00      1.00      1.00      1.00
S2            0.15939      0.03738      0.74436       0.00      0.99      0.00      1.00      1.00      1.00
S3            0.16421      0.00126      0.90283       0.00      0.97      0.00      0.99      1.00      1.00
S4            0.20566      0.00927      1.26650       0.00      0.45      0.00      0.93      1.00      1.00
S5            0.22648      0.00203      0.95154       0.00      0.56      0.00      0.38      1.00      1.00
S6            0.27046      0.05076      1.39838       0.00      0.09(!)   0.00    -59.58(!)   1.00      1.00
S7            0.29500      0.00364      1.20590       0.00      1.00      0.00      0.80      1.00      1.00
S8            0.30767      0.00476      1.62248       0.00      0.87      0.00      1.00      0.94      1.00
```

Therefore, it is a good idea to remove `S6` and fit again:

```bash
./BTevalSpec.py rm S6
./BTevalSpec.py --signif fit
```

This produces the add-line objective

![addLineObj2_wref1ß](readmeFiles/addLineObj4_wref10_afterRmS6Refit.png   "Final add-line objective (wref=10)")

The comprehensive fit error hardly changes after removing `S6`, however, the significances are now reasonably:

```text
# name  | energy[Ry] |  strength  | phase[pi]  | ... strengthError  | signifFit|signifAng|signifExc|signifErr|signifRng|signifPha
S1            0.14912      0.01485      0.40642  ...     0.0000212     1.00      0.93      1.00      1.00      1.00      1.00
S2            0.15939      0.03738      0.74436          0.0000184     1.00      0.99      1.00      1.00      1.00      1.00
S3            0.16421      0.00126      0.90282          0.0000262     1.00      0.97      0.99      0.99      1.00      1.00
S4            0.20566      0.00927      1.26650          0.0010854     1.00      0.45(c)   0.98      0.93      1.00      1.00
S5            0.22648      0.00203(b)   0.95156          0.0021873(b)  1.00      0.56(c)   0.93      0.38(a)   1.00      1.00
S6            0.29500      0.00365      1.20588          0.0012226     1.00      1.00      0.98      0.81      1.00      1.00
S7            0.30767      0.00477      1.62261          0.0000267     1.00      0.87      0.98      1.00      0.93      1.00
```

`S5` still shows a small `signifErr` (see "(a)" marker), which stems from the relatively large fit error from the line-height fit (compare `strengthErr` vs. `strength` for this excitation, "(b)" marker).
However, the excitation seems to be real concerning the other significances.
One could test this further by removing `S5` and re-fitting.

The smaller value of `sgnifAng` of `S4` and `S5` (see "(c)" marker) is not problematic if the other significances are fine.
Adding more excitations, however, would be pure speculation.