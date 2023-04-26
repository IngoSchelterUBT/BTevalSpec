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
import handleTrace
import errorHandler as err

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():

    #--------------------------------------------------------------------------#
    # In Future:
    # Process command-line arguments
    #--------------------------------------------------------------------------#
    ifile = "eval.yaml"

    #--------------------------------------------------------------------------#
    # Read configuration from eval.yaml
    #--------------------------------------------------------------------------#
    if not os.path.isfile(ifile):
        config.writeTemplate(ifile)
        err.err(1,'There was no '+ifile+' file, now a template is created!')
    conf = config.Config(ifile)

    #--------------------------------------------------------------------------#
    # Read dipole moments into dip[ncalc][narea] list:
    #--------------------------------------------------------------------------#
    dip = []
    for icalc in range(len(conf.dipfiles)):
        dip.append([])
        for iarea in range(len(conf.dipfiles[icalc])):
            dip[icalc].append(dipole.Dipole(conf.dipfiles[icalc][iarea]))

    #--------------------------------------------------------------------------#
    # Fourier transform
    #--------------------------------------------------------------------------#
    if conf.opt["FT"]["calc"]:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].ft(minpw=conf.opt["FT"]["minpw"],window=conf.opt["FT"]["window"],smooth=conf.opt["FT"]["smooth"],rmDC=True))
    elif conf.opt["Fit"]["calc"] or conf.opt["Fit"]["plot_results"]:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["ft","pw"])
        
    #--------------------------------------------------------------------------#
    # Pade approximation
    #--------------------------------------------------------------------------#
    pade = []
    if conf.opt["Pade"]["calc"]:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].pade(conf.opt["Pade"]["wmax"],conf.opt["Pade"]["dw"],thin=conf.opt["Pade"]["thin"],smooth=conf.opt["Pade"]["smooth"]))
    elif conf.opt["Fit"]["guess"]:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["pade"])

    #--------------------------------------------------------------------------#
    # Write spectra
    #--------------------------------------------------------------------------#
    for icalc in range(len(dip)):
        for iarea in range(len(dip[icalc])):
            if conf.opt["FT"  ]["calc"]: dip[icalc][iarea].writeSpectra(what=["ft","pw"])
            if conf.opt["Pade"]["calc"]: dip[icalc][iarea].writeSpectra(what=["pade"]   )

#    #--------------------------------------------------------------------------#
#    # Create initial guess for the spectrum fit
#    #--------------------------------------------------------------------------#
#    guess  = []
#    future = []
#    if conf.fit or conf.plot_result:
#        with concurrent.futures.ThreadPoolExecutor() as executer:
#            for icalc in range(len(dip)):
#                guess.append([])
#                future.append([])
#                for iarea in range(len(dip[icalc])):
#                    future[i].append(executer.submit(specGuess.Guess,dip[icalc][iarea].ft, dip[icalc][iarea].pade,...))
#                    guess[i].append(future[i].result())
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
