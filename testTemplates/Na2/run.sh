#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=1 fit
./BTevalSpec.py --wref=1.0 --skip --niter=1 --signif fit
