# Shortcut

Run `./run.sh`

# Fitting

- Get the initial input file: `cp readmeFiles/eval.yaml.0 ./eval.yaml` or generate one by calling `./BTevalSpec.py new` and add the dipole-moment file and change the fitrange to `0.1:0.28`.
- Fit, e.g., doing three iterations by calling

```bash
./BTevalSpec.py --niter=3 fit
```

  The latter will produce the following add-line objectives (one after the other, note how the lines get smaller)
  ![addLineObj0](readmeFiles/addLineObj0.png   "Add-line objective (for fit iteration 1, initial)")
  ![addLineObj1](readmeFiles/addLineObj1.png   "Add-line objective (for fit iteration 2)")
  ![addLineObj2](readmeFiles/addLineObj2.png   "Add-line objective (for fit iteration 3)")

- There is some wiggly error left, therefore one can try suppress the line-shape errors of already fitted lines using

```bash
./BTevalSpec.py --wref=1. --skip fit
```

The latter produces the following add-line objective, which does not show an proper additional line.
  ![addLineObj3](readmeFiles/addLineObj3_wref.png   "Add-line objective (for simulating fit iteration 4 with error-supression)")

Note that the excitation `S9` from `Na4_laser_transdens` test is missing here.
