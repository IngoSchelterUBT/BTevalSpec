#!/usr/bin/env python
#------------------------------------------------------------------------------#
# author: Rian Richter
# Fit of spectrum
import numpy as np
import re
import sys
import os.path
import getopt
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#only for tests
from pprint import pprint


#Import Modules
import config
import dipole

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():
  #Inizialize the object containing all the configuration information in eval.yaml
  conf = config.Config()
  #Inizialize a list of objects containing all configurations of all dipole files
  dip = []
  for i, fileName in enumerate(conf.dipoleFiles):
    dip.append(dipole.Dipole(fileName))

  pprint(vars(dip[0]))

  #Do Fourier Transformation of the dipole moment file(s)


#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
