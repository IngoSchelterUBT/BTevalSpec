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
import specGuess
import specFit
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
      dip.append(dipole.Dipole(fileName,i))
    if conf.numDipoleFiles == 3:
      dip.append(dipole.Dipole(conf.dipoleFiles,3,calcFlag='trace'))
  elif not conf.pade and not conf.fourier and (conf.fit or conf.fit_guess):
    #Construct dummy dip list. In this case the dipole moment is not needed
    #to be read.
    for i in range(conf.numDipoleFiles+1):
      dip.append(0)


  #Do Fourier Transformation of the dipole moment file(s) of if only fit is true
  #than read the Fourier Transformation
  ft = []
  if conf.fourier or conf.fit:
    if conf.fourier: inout.cleanFT() #delete Osci and PW folder
    for i, fileName in enumerate(conf.dipoleFiles):
      ft.append(fourier.FT(conf,dip[i],i))
    if conf.numDipoleFiles == 3:
      ft.append(fourier.FT(conf,dip[3],3,calcFlag='trace'))

  #Do Pade Approximation of the dipole moment file(s)
  pade = []
  if conf.pade or conf.fit_guess:
    if conf.pade: inout.cleanPade() #delete PADE folder
    for i, fileName in enumerate(conf.dipoleFiles):
      pade.append(padeApprox.Pade(conf,dip[i],i))
    if conf.numDipoleFiles == 3:
      pade.append(padeApprox.Pade(conf,dip[3],3,calcFlag='trace'))
  elif not conf.pade and conf.fit and not conf.fit_guess:
    #construct dummy pade list
    for i in range(conf.numDipoleFiles+1):
      pade.append(0)
      

  #Do Guess for fit of the spectrum
  if conf.fit:
    guess = []
    for i, fileName in enumerate(conf.dipoleFiles):
      guess.append(specGuess.Guess(conf,ft[i],pade[i],i))
    if conf.numDipoleFiles == 3:
      guess.append(specGuess.Guess(conf,ft[3],pade[3],3,calcFlag='trace'))

  #Do Fit of the spectrum
  if conf.fit:
    fit = []
    for i, fileName in enumerate(conf.dipoleFiles):
      fit.append(specFit.Fit(conf,ft[i],guess[i],i))
    if conf.numDipoleFiles == 3:
      fit.append(specFit.Fit(conf,ft[3],guess[3],3,calcFlag='trace'))


  input("Press [enter] to end and close all plots!")
#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
