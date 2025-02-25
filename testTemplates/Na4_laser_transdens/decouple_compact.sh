#!/bin/bash

#Unpack densft compact files
tar -xzf readmeFiles/densft_compact.tar.gz
#Prepend DENSFT section to eval.yaml
cat eval.yaml readmeFiles/eval.yaml.densft_compact >> eval.yaml.tmp
mv eval.yaml.tmp eval.yaml
#Call BTevalSpec
./BTevalSpec.py decouple
#Clean up
rm densft*.compact
