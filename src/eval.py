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


#Import own Modules
import config
import dipole
import fourierTransform as fourier
import padeApprox
import inout

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():
  #write template config
  inout.writeEmptyConfig()
  #Inizialize the object containing all the configuration information in eval.yaml
  conf = config.Config()
  #Inizialize a list of objects containing all configurations of all dipole files
  dip = []
  if conf.fourier or conf.pade:
    for i, fileName in enumerate(conf.dipoleFiles):
      dip.append(dipole.Dipole(fileName,i+1))
  elif not conf.pade and not conf.fourier and (conf.fit or conf.fit_guess):
    #Construct dummy dip list. In this case the dipole moment is not needed
    #to be read.
    for i, fileName in enumerate(conf.dipoleFiles):
      dip.append(0)

  #pprint(vars(dip[0]))

  #Do Fourier Transformation of the dipole moment file(s) of if only fit is true
  #than read the Fourier Transformation
  if conf.fourier or conf.fit:
    ft = []
    for i, fileName in enumerate(conf.dipoleFiles):
      ft.append(fourier.FT(conf,dip[i],i+1))

  #Do Pade Approximation of the dipole moment file(s)
  #To Do: insert boolean for making guess in fit
  if conf.pade or conf.fit_guess:
    pade = []
    for i, fileName in enumerate(conf.dipoleFiles):
      pade.append(padeApprox.Pade(conf,dip[i],i+1))


#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
