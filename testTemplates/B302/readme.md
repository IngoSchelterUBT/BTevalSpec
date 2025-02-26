# Description

This test tries to fit as many excitations as possible and test the capabilities of this script.
The dipole moment stems from a boost calculation on a bacteriochlorophyll a.

# Shortcut

Run `./run.sh`

# Fitting

- Get the initial input file: `cp readmeFiles/eval.yaml.0 ./eval.yaml`.
- Fit, e.g., calling

```bash
./BTevalSpec.py --niter=4 fit
./BTevalSpec.py --niter=4 --wref=0.1 fit
```

  The latter will produce the following add-line objectives (one after the other, note how the lines get smaller) and results with 56 excitations in a comprehensive fit error of `0.011`.
  ![addLineObj0](readmeFiles/addLineObj0.png   "Add-line objective (for fit iteration 1, initial)")
  ![addLineObj1](readmeFiles/addLineObj1.png   "Add-line objective (for fit iteration 2)")
  ![addLineObj2](readmeFiles/addLineObj2.png   "Add-line objective (for fit iteration 3)")
  ![addLineObj3](readmeFiles/addLineObj3.png   "Add-line objective (for fit iteration 4)")
  ![addLineObj4](readmeFiles/addLineObj4_wrefE-1.png   "Add-line objective (for fit iteration 5, error supression with wref=0.1)")
  ![addLineObj5](readmeFiles/addLineObj5_wrefE-1.png   "Add-line objective (for fit iteration 6, error supression with wref=0.1)")
  ![addLineObj6](readmeFiles/addLineObj6_wrefE-1.png   "Add-line objective (for fit iteration 7, error supression with wref=0.1)")
  ![addLineObj7](readmeFiles/addLineObj7_wrefE-1.png   "Add-line objective (for fit iteration 8, error supression with wref=0.1)")

- During the fit, the computed significances look good. However, the last iteration takes a while, which can indicate that excitations are insignificant.
- Computing the significances reveals overall good significances with the exception of three excitations `S2`, `S40`, and `S41` (cf. `readmeFiles/excit_1.dat.8.signif` and `readmeFiles/eval.yaml.8.signif`).
  If two adjacent lines have low significances, this often indicates that one of them is redundant. Remove `S2` and `S40`and re-fit.

```bash
./BTevalSpec.py rm S2,S40
./BTevalSpec.py --signif fit
```

This produces the files `readmeFiles/eval.yaml.final` and `readmeFiles/excit_1.dat.final` with 53 excitations.

In principle, one could now go on and try adding more excitations, e.g., with `wref=1.0`. Also, `S26` and `S27` have `signifFit` values of about `0.5`, which could indicate that one is redundant.

```bash
./BTevalSpec.py rm S27        #remove S27
./BTevalSpec.py --signif fit  #re-fit
./BTevalSpec.py --wref=1. fit #Simulate addint new excitations with high error supression
...
```