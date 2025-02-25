#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=4 fit
./BTevalSpec.py --niter=1 --wref=0.01 --signif --skip fit
