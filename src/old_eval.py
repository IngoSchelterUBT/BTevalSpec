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
#matplotlib.use("TkAgg")
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
import handleTrace

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
    dip = [None]*len(conf.dipoleFiles)
    for i, fileName in enumerate(conf.dipoleFiles):
      dip[i] = dipole.Dipole(fileName,i)
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
    ft = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: ft.append(None)
    if conf.fourier: inout.cleanFT() #delete Osci and PW folder
    for i in range(len(ft)):
      if i == 3:
        calcFlag = 'trace'
      else:
        calcFlag = 'no'
      ft[i] = fourier.FT(conf,dip[i],i,calcFlag)

  #Do Pade Approximation of the dipole moment file(s)
  pade = []
  if conf.pade or conf.fit_guess:
    pade = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: pade.append(None)
    if conf.pade: inout.cleanPade() #delete PADE folder
    for i in range(len(pade)):
      if i == 3:
        calcFlag = 'trace'
      else:
        calcFlag = 'no'
      pade[i] = padeApprox.Pade(conf,dip[i],i,calcFlag)
  elif not conf.pade and conf.fit and not conf.fit_guess:
    #construct dummy pade list
    for i in range(conf.numDipoleFiles+1):
      pade.append(0)
      

  #Do Guess for fit of the spectrum
  if conf.fit:
    guess = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: guess.append(None)
    for i in range(len(guess)):
      if i == 3:
        calcFlag = 'trace'
      else:
        calcFlag = 'no'
      guess[i] = specGuess.Guess(conf,ft[i],pade[i],i,calcFlag)

  #Do Fit of the spectrum
  if conf.fit:
    fit = [None]*len(conf.dipoleFiles)
    for i, fileName in enumerate(conf.dipoleFiles):
      fit[i] = specFit.Fit(conf,ft[i],guess[i],i)
    if conf.numDipoleFiles == 3:
      handleTrace.guessTrace(guess[3], fit) 
      fit.append(specFit.Fit(conf,ft[3],guess[3],3,calcFlag='trace'))

  #write excitation lines
  inout.writeExcitations(conf,fit)


  input("Press [enter] to end and close all plots!")

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
