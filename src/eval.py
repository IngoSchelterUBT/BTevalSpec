#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# author: Rian Richter
# Fit of spectrum
import numpy as np
#import re
import sys
import os.path
#import getopt
#import matplotlib
#matplotlib.use("Qt5Agg")
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import concurrent.futures

#only for tests
#from pprint import pprint

#Import own Modules
import errorHandler as err
import config
import extern
import dipole
import excitations
import fit
#import specGuess
#import specFit
#import spectrum
#import handleTrace

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
    conf  = config.Config(ifile)
    excit = excitations.Excitations(ncalc=len(conf.dipfiles),narea=len(conf.dipfiles[0]),ncomp=3,exlist=conf.excit)

    #--------------------------------------------------------------------------#
    # Read dipole moments into dip[ncalc][narea] list:
    # In case of an external excitation with t_end>0:
    # Move the time frame so that it starts at t_end
    #--------------------------------------------------------------------------#
    dip = []
    for icalc in range(len(conf.dipfiles)):
        dip.append([])
        for iarea in range(len(conf.dipfiles[icalc])):
            dip[icalc].append(dipole.Dipole(conf.dipfiles[icalc][iarea],["x","y","z"]))

    #--------------------------------------------------------------------------#
    # Fourier transform
    #--------------------------------------------------------------------------#
    if conf.opt["FT"]["calc"]:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getFt(minpw=conf.opt["FT"]["minpw"],window=conf.opt["FT"]["window"],smooth=conf.opt["FT"]["smooth"],rmDC=True))
                    #dip[icalc][iarea].getFt(minpw=conf.opt["FT"]["minpw"],window=conf.opt["FT"]["window"],smooth=conf.opt["FT"]["smooth"],rmDC=True)
                    dip[icalc][iarea].writeSpectra(what=["ft","pw"])
        conf.opt["FT"]["calc"] = False #Next time: read by default
    else: #elif conf.opt["Fit"]["calc"] or conf.opt["Fit"]["plot_results"]:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["ft","pw"])
        
    #--------------------------------------------------------------------------#
    # Pade approximation
    #--------------------------------------------------------------------------#
    if conf.opt["Pade"]["calc"]:
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getPade(conf.opt["Pade"]["wmax"],conf.opt["Pade"]["dw"],thin=conf.opt["Pade"]["thin"],smooth=conf.opt["Pade"]["smooth"]))
                    #dip[icalc][iarea].getPade(conf.opt["Pade"]["wmax"],conf.opt["Pade"]["dw"],thin=conf.opt["Pade"]["thin"],smooth=conf.opt["Pade"]["smooth"])
                    dip[icalc][iarea].writeSpectra(what=["pade"])
        conf.opt["Pade"]["calc"] = False #Next time: read by default
    else: #elif conf.opt["Fit"]["guess"]:
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["pade"])

    #--------------------------------------------------------------------------#
    # Read and Fourier transform external excitation
    #--------------------------------------------------------------------------#
    ext = extern.Extern(conf.ext.get("profile",""))
    #if isinstance(ext,extern.Extern):
    ext.write()

    #--------------------------------------------------------------------------#
    # Create initial guess for the spectrum fit
    #--------------------------------------------------------------------------#
    dfit  = fit.Fit(dip,ext,excit,conf.opt["Fit"]["range"])
    if conf.opt["Fit"]["guess"]:
        excit = dfit.newGuess()
        conf.opt["Fit"]["guess"] = False #Next time: No new initial guess

#    if conf.opt["Fit"]["calc"] and conf.opt["Fit"]["guess"]: # or conf.opt["Fit"]["plot_result"]:
#        dfit.guess(init=excit,thres=conf.opt["Fit"]["guess_thres"])

#    #--------------------------------------------------------------------------#
#    # Fit
#    #--------------------------------------------------------------------------#
#    if conf.opt["Fit"]["calc"]:
#        fit.fit(init=excit,crit=conf.opt["Fit"]["relerr_crit"],maxIter=conf.opt["Fit"]["max_iter"])

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

    #--------------------------------------------------------------------------#
    # Update configuration file
    #--------------------------------------------------------------------------#
    conf.excit = [excit.exlist[i].todict() for i in range(len(excit.exlist))]
    conf.write(ifile)

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
