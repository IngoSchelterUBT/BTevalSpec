#!/bin/bash

tar -xzf readmeFiles/transd.tar.gz
cp readmeFiles/eval.yaml.0 ./eval.yaml
./BTevalSpec.py --niter=5 --signif fit
./BTevalSpec.py --imag decouple
