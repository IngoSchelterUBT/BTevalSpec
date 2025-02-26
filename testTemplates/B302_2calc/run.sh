#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=4 fit
./BTevalSpec.py --niter=4 --wref=0.1 fit
./BTevalSpec.py --skip --signif fit
#./BTevalSpec.py rm S2,S40,S41
#./BTevalSpec.py fit
#./BTevalSpec.py --wref=1.0 fit #inspect add-line objective with higher error suppression
