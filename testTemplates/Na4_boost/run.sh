#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=3 --signif fit
./BTevalSpec.py --wref=1.0 --skip fit #Don't add a line, just see if error-suppression helps
