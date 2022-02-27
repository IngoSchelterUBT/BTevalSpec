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
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import concurrent.futures

#only for tests
from pprint import pprint

#Import own Modules
import config
import dipole
import fourierTransform as fourier
import padeApprox
import specGuess
import specFit
import spectrum
import inout
import handleTrace
import errorHandler as err

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():
  #write template config
  if not os.path.isfile('eval.yaml'):
    inout.writeEmptyConfig()
    err.err(1,'There was no eval.yaml file, now a template is created!')
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
  else:
    #Construct dummy dip list. In this case the dipole moment is not needed
    #to be read.
    for i in range(conf.numDipoleFiles+1):
      dip.append(0)


  #Do Fourier Transformation of the dipole moment file(s) of if only fit is true
  #than read the Fourier Transformation
  ft = []
  if conf.fourier or conf.fit or conf.plot_result:
    ft = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: ft.append(None)
    if conf.fourier: inout.cleanFT() #delete Osci and PW folder
    future = ft #create a future object list
    with concurrent.futures.ThreadPoolExecutor() as executer:
      for i in range(len(ft)):
        if i == 3:
          calcFlag = 'trace'
        else:
          calcFlag = 'no'
        future[i] = executer.submit(fourier.FT, conf, dip[i], i, calcFlag)
        ft[i] = future[i].result()
        #ft[i] = fourier.FT(conf,dip[i],i,calcFlag)

  #Do Pade Approximation of the dipole moment file(s)
  pade = []
  if conf.pade or conf.fit_guess:
    pade = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: pade.append(None)
    if conf.pade: inout.cleanPade() #delete PADE folder
    future = pade #create a future object list
    with concurrent.futures.ThreadPoolExecutor() as executer:
        for i in range(len(pade)):
          if i == 3:
            calcFlag = 'trace'
          else:
            calcFlag = 'no'
          future[i] = executer.submit(padeApprox.Pade, conf, dip[i], i, calcFlag)
          pade[i] = future[i].result()
          #pade[i] = padeApprox.Pade(conf,dip[i],i,calcFlag)
  else:
    #construct dummy pade list
    for i in range(conf.numDipoleFiles+1):
      pade.append(0)


  #Do Guess for fit of the spectrum
  if conf.fit or conf.plot_result:
    guess = [None]*len(conf.dipoleFiles)
    if conf.numDipoleFiles == 3: guess.append(None)
    future = guess #create a future object list
    with concurrent.futures.ThreadPoolExecutor() as executer:
        for i in range(len(guess)):
          if i == 3:
            calcFlag = 'trace'
          else:
            calcFlag = 'no'
          future[i] = executer.submit(specGuess.Guess, conf, ft[i], pade[i], i, calcFlag)
          guess[i] = future[i].result()
          #guess[i] = specGuess.Guess(conf,ft[i],pade[i],i,calcFlag)


  #Do Fit of the spectrum
  if conf.fit or conf.plot_result:
    fit = [None]*len(conf.dipoleFiles)
    future = fit
    with concurrent.futures.ThreadPoolExecutor() as executer:
        for i, fileName in enumerate(conf.dipoleFiles):
          future[i] = executer.submit(specFit.Fit, conf, ft[i], guess[i], i)
          fit[i] = future[i].result()
    if conf.numDipoleFiles == 3:
      if conf.fit: handleTrace.guessTrace(conf, guess[3], fit)
      fit.append(specFit.Fit(conf,ft[3],guess[3],3,calcFlag='trace'))
    
    #plotting the results has to be unparalleled (problem with starting
    #matplotlib gui)
    for i, f in enumerate(fit):
      if i == 3:
        calcFlag = 'trace'
      else:
        calcFlag = 'no'
      f.plotFit(calcFlag)

    #create object spectrum
    spec = spectrum.Spectrum(conf,fit)
    #write excitation lines
    if conf.fit: inout.writeExcitations(spec)
    input("Press [enter] to end and close all plots!")

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
