#!/bin/bash

#Unpack densft cube files
tar -xzf readmeFiles/densft_cube.tar.gz
#Prepend DENSFT section to eval.yaml
cat readmeFiles/eval.yaml.densft_cube eval.yaml >> eval.yaml.tmp
mv eval.yaml.tmp eval.yaml
#Call BTevalSpec
./BTevalSpec.py decouple
#Clean up
rm densft*.cube
