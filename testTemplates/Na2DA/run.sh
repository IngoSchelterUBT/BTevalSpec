#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=1 fit
./BTevalSpec.py --niter=1 --wref=1  --skip fit
./BTevalSpec.py --niter=1 --wref=10 --skip fit
./BTevalSpec.py rm S6
./BTevalSpec.py --signif fit

