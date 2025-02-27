#!/bin/bash

cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=4 fit
./BTevalSpec.py --niter=4 --wref=0.1 fit #So far, no indication of insignificant lines
./BTevalSpec.py --skip --signif fit      #Just compute significances
#S2, S40, and S41 have lower significances, but not completely unreasonably -> remove and fit again
./BTevalSpec.py rm S2,S40,S41
./BTevalSpec.py fit
./BTevalSpec.py --wref=1.0 fit #inspect add-line objective with higher error suppression
