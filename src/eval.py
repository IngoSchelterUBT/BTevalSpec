#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# author: Rian Richter
# Fit of spectrum
import numpy as np
import re
import sys
import os.path
import getopt
#import matplotlib
#matplotlib.use("Qt5Agg")
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import concurrent.futures

#only for tests
from pprint import pprint

#Import own Modules
import config
import dipole
#import specGuess
#import specFit
#import spectrum
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
        for iarea in range(len(conf.files[icalc])): #Go through areas
            dip[icalc].append(dipole.Dipole(conf.files[icalc][iarea]))
#        dipg.append(dipole.sum(dip[icalc]))
#    if conf.numDipoleCalc==3.and.conf.numDipoleArea==1: #Assume boost in x/y/z directions
#        trace   = dipole.trace([dip[0][0],dip[1][0],dip[2][0])

    #--------------------------------------------------------------------------#
    # Fourier transform
    #--------------------------------------------------------------------------#
    if conf.fourier:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].ft(minpw=conf.minpw,window=conf.window,smooth=conf.smooth,rmDC=True))
#                executer.submit(dipg[icalc].ft(minpw=conf.minpw,window=conf.window,smooth=conf.smooth,rmDC=True))
#            executer.submit(trace.ft(minpw=conf.minpw,window=conf.window,smooth=conf.smooth,rmDC=True))
    elif conf.fit or conf.plot_results:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["ft","pw"])
        
    #--------------------------------------------------------------------------#
    # Pade approxiimation
    #--------------------------------------------------------------------------#
    pade = []
    if conf.pade:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].pade(conf.pade_wmax,conf.pade_dw,thin=conf.pade_thin,smooth=conf.pade_smooth))
#                executer.submit(dipg[icalc].pade)
#            executer.submit(trace.pade)
    elif conf.fit_guess:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["pade"])

    #--------------------------------------------------------------------------#
    # Write spectra
    #--------------------------------------------------------------------------#
    for icalc in range(len(dip)):
        for iarea in range(len(dip[icalc])):
            if conf.fourier: dip[icalc][iarea].writeSpectra(what=["ft","pw"])
            if conf.pade   : dip[icalc][iarea].writeSpectra(what=["pade"]   )
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
