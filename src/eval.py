#!/usr/bin/env python3
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
# Read dipole files
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():

    #--------------------------------------------------------------------------#
    # Read configuration from eval.yaml
    #--------------------------------------------------------------------------#
    if not os.path.isfile('eval.yaml'):
        inout.writeEmptyConfig()
        err.err(1,'There was no eval.yaml file, now a template is created!')
    conf = config.Config()

    #--------------------------------------------------------------------------#
    # Read dipole moments into dip[ncalc][narea] list:
    #--------------------------------------------------------------------------#
    #Initialize a list of objects containing all configurations of all dipole files
    dip = []
    #dipg= []
    for icalc in range(len(conf.files)):                #Go through calculations
        dip.append([])
        if conf.fourier or conf.pade:
            for iarea in range(len(conf.files[icalc])): #Go through areas
                dip[icalc].append(dipole.Dipole(conf.files[icalc][iarea]))
        else: #No nead to read in this case
            dip[icalc].append([None]*len(conf.files[icalc]))
#        dipg.append(dipole.sum(dip[icalc]))
#    if conf.numDipoleCalc==3.and.conf.numDipoleArea==1: #Assume boost in x/y/z directions
#        trace   = dipole.trace([dip[0][0],dip[1][0],dip[2][0])

    #--------------------------------------------------------------------------#
    # Fourier transform
    # If only fit is true, then read the Fourier Transformation instead
    #--------------------------------------------------------------------------#
    if conf.fourier or conf.fit or conf.plot_result:
        if conf.fourier: inout.cleanFT() #delete Osci and PW folder
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(conf.files)):            #Go through calculations
                for iarea in range(len(conf.files[icalc])): #Go through areas
                    if conf.fourier:
                        executer.submit(dip[icalc][iarea].ft(minpw=conf.minpw,window=conf.window,smooth=conf.smooth,rmDC=True)) #, conf, dip[i], i, calcFlag)
#                executer.submit(dipg[icalc].ft)
#            executer.submit(trace.ft)

    #--------------------------------------------------------------------------#
    # Write spectra
    #--------------------------------------------------------------------------#
    for icalc in range(len(conf.files)):            #Go through calculations
        for iarea in range(len(conf.files[icalc])): #Go through areas
            dip[icalc][iarea].writeSpectra()

#    #Do Pade Approximation of the dipole moment file(s)
#    pade = []
#    if conf.pade or conf.fit_guess:
#        pade = [None]*len(conf.files)
#        if conf.numDipoleFiles == 3: pade.append(None)
#        if conf.pade: inout.cleanPade() #delete PADE folder
#        future = pade #create a future object list
#        with concurrent.futures.ThreadPoolExecutor() as executer:
#                for i in range(len(pade)):
#                    if i == 3:
#                        calcFlag = 'trace'
#                    else:
#                        calcFlag = 'no'
#                    future[i] = executer.submit(padeApprox.Pade, conf, dip[i], i, calcFlag)
#                    pade[i] = future[i].result()
#                    #pade[i] = padeApprox.Pade(conf,dip[i],i,calcFlag)
#    else:
#        #construct dummy pade list
#        for i in range(conf.numDipoleFiles+1):
#            pade.append(0)
#
#
#    #Do Guess for fit of the spectrum
#    if conf.fit or conf.plot_result:
#        guess = [None]*len(conf.files)
#        if conf.numDipoleFiles == 3: guess.append(None)
#        future = guess #create a future object list
#        with concurrent.futures.ThreadPoolExecutor() as executer:
#                for i in range(len(guess)):
#                    if i == 3:
#                        calcFlag = 'trace'
#                    else:
#                        calcFlag = 'no'
#                    future[i] = executer.submit(specGuess.Guess, conf, ft[i], pade[i], i, calcFlag)
#                    guess[i] = future[i].result()
#                    #guess[i] = specGuess.Guess(conf,ft[i],pade[i],i,calcFlag)
#
#
#    #Do Fit of the spectrum
#    if conf.fit or conf.plot_result:
#        fit = [None]*len(conf.files)
#        future = fit
#        with concurrent.futures.ThreadPoolExecutor() as executer:
#                for i, fileName in enumerate(conf.files):
#                    future[i] = executer.submit(specFit.Fit, conf, ft[i], guess[i], i)
#                    fit[i] = future[i].result()
#        if conf.numDipoleFiles == 3:
#            if conf.fit: handleTrace.guessTrace(conf, guess[3], fit)
#            fit.append(specFit.Fit(conf,ft[3],guess[3],3,calcFlag='trace'))
#        
#        #plotting the results has to be unparalleled (problem with starting
#        #matplotlib gui)
#        for i, f in enumerate(fit):
#            if i == 3:
#                calcFlag = 'trace'
#            else:
#                calcFlag = 'no'
#            f.plotFit(calcFlag)
#
#        #create object spectrum
#        spec = spectrum.Spectrum(conf,fit)
#        #write excitation lines
#        if conf.fit: inout.writeExcitations(spec)
#        input("Press [enter] to end and close all plots!")

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
