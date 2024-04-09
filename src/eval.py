#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# author: Rian Richter & Ingo Schelter
#------------------------------------------------------------------------------#
#
# TODO:
# - Implement command-line-argument handling
# - Go through each excitation and warn if there are close excitations with same dipole direction. This could be due to misshaped lines (e.g., non-linear effects due to too strong excitation)
# - Warn if excitations are insignificant (e.g, small oscillator strength, strong correlations between excitations)
# - Spectrum plot with amplitudes instead of osci strengths if the latter are not significant (e.g., if excitation restricted to an area)
# - Print #number of excitations and error after fit into file
# - If manually adding an excitation:
#  -> Note this automatically
#  -> Set skipfirst to false
#  -> Make a proper guess for the dipole moments and phase
# - Log commands and iterations automatically (save intermediate eval.yaml and dipole*fit.dat etc.)
#  -> e.g., automatically write and update a skript with all commands that are executed
#
#
# See errorHandler.py for usage information
#------------------------------------------------------------------------------#
import numpy as np
import sys
import os.path
#import matplotlib
#matplotlib.use("Qt5Agg")
#import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import concurrent.futures

#Import own Modules
import errorHandler as err
import config
import extern
import dipole
import excitations
import fit

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():

    #--------------------------------------------------------------------------#
    # In Future:
    # Process command-line arguments
    #--------------------------------------------------------------------------#
    ifile = "eval.yaml"
    dbg   = 1

    #--------------------------------------------------------------------------#
    # Read configuration from eval.yaml
    #--------------------------------------------------------------------------#
    if dbg>0: print("Read configuration",end="")
    if not os.path.isfile(ifile):
        config.writeTemplate(ifile)
        err.err(1,'There was no '+ifile+' file, now a template is created!')
    conf  = config.Config(ifile)
    if dbg>0: print(" - done")


    #--------------------------------------------------------------------------#
    # Read dipole moments into dip[ncalc][narea] list:
    # In case of an external excitation with t_end>0:
    # Move the time frame so that it starts at t_end
    #--------------------------------------------------------------------------#
    if dbg>0: print("Read dipole moment files",end="")
    dip = []
    for icalc in range(len(conf.dipfiles)):
        dip.append([])
        for iarea in range(len(conf.dipfiles[icalc])):
            dip[icalc].append(dipole.Dipole(conf.dipfiles[icalc][iarea],["x","y","z"]))
    if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Fourier transform
    #--------------------------------------------------------------------------#
    if conf.opt["FT"]["calc"]:
        if dbg>0: print("Calculate Fourier transform",end="")
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getFt(minpw=conf.opt["FT"]["minpw"],window=conf.opt["FT"]["window"],smooth=conf.opt["FT"]["smooth"],rmDC=True))
                    dip[icalc][iarea].writeSpectra(what=["ft","pw"])
        conf.opt["FT"]["calc"] = False #Next time: read by default
        if dbg>0: print(" - done")
    else: #elif conf.opt["Fit"]["calc"] or conf.opt["Fit"]["plot_results"]:
        if dbg>0: print("Read Fourier transform",end="")
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["ft","pw"])
        if dbg>0: print(" - done")
        
    #--------------------------------------------------------------------------#
    # Pade approximation
    #--------------------------------------------------------------------------#
    if conf.opt["Pade"]["calc"]:
        if dbg>0: print("Calculate Pade approximation",end="")
        if not conf.opt["Pade"].get("smooth",0.) > 0.: #Choose automatic smoothing
            conf.opt["Pade"]["smooth"] = float(np.log(100)/dip[0][0].tprop) # e^-s*T=0.01 <=> s=ln(100)/T
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getPade(conf.opt["Pade"]["wmax"],conf.opt["Pade"]["dw"],thin=conf.opt["Pade"]["thin"],smooth=conf.opt["Pade"]["smooth"]))
                    dip[icalc][iarea].writeSpectra(what=["pade"])
        conf.opt["Pade"]["calc"] = False #Next time: read by default
        if dbg>0: print(" - done")
    else: #elif conf.opt["Fit"]["guess"]:
        if dbg>0: print("Read Pade approximation",end="")
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["pade"])
        if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Read and Fourier transform external excitation
    #--------------------------------------------------------------------------#
    if dbg>0: print("Initialize extern potential",end="")
    ext = extern.Extern(conf.ext.get("profile",""),conf.ext.get("invertPhase",False),dip[0][0].efield,dip[0][0].text,[dip[icalc][0].epol for icalc in range(len(dip))])
    #if isinstance(ext,extern.Extern):
    ext.write()
    if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Initialize excitations structure
    #--------------------------------------------------------------------------#
    if dbg>0: print("Initialize excitations",end="")
    excit = excitations.Excitations(ncalc=len(conf.dipfiles),narea=len(conf.dipfiles[0]),ncomp=3,ext=ext,exlist=conf.excit)
    if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Create initial guess for the spectrum fit
    #--------------------------------------------------------------------------#
    dfit = fit.Fit(dip,ext,excit,conf.opt["Fit"]["range"])
    if conf.opt["Fit"]["guess"]:
        if dbg>0: print("Initial guess",end="")
        excit = dfit.newGuess(hf=conf.opt["Fit"]["guess_thres"],dbg=dbg)
        conf.opt["Fit"]["guess"] = False #Next time: No new initial guess
        if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Fit
    #--------------------------------------------------------------------------#
    if conf.opt["Fit"]["calc"]:
        if dbg>0: print("Fitting:")
        excit, fiterr = dfit.fit(dbg=dbg,tol=conf.opt["Fit"]["relerr_crit"],maxex=conf.opt["Fit"]["max_excit"],skipfirst=conf.opt["Fit"].get("skipfirst",False),allSignif=conf.opt["Fit"].get("significances",False),nsigma=conf.opt["Fit"].get("nsigma",2.),firstsingle=conf.opt["Fit"].get("firstsingle",False),resetErange=conf.opt["Fit"].get("reset_erange",False),fitphase=conf.opt["Fit"].get("fitphase",True))
        conf.opt["Fit"]["skipfirst"]    = True
        conf.opt["Fit"]["firstsingle"]  = False
        conf.opt["Fit"]["reset_erange"] = False
        conf.opt["Fit"]["fiterr"]       = float(fiterr)
        if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Update configuration file
    #--------------------------------------------------------------------------#
    if dbg>0: print("Update new configuration file")
    conf.excit = [excit.exlist[i].todict() for i in range(len(excit.exlist))]
    conf.write(ifile) 
    if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Write Fit files
    #--------------------------------------------------------------------------#
    if dbg>0: print("Write fit files")
    if conf.opt["Fit"]["calc"]: dfit.writeFit(dbg=dbg)
    excit.excitFiles(dip[0][0].tprop)
    if dbg>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Plot spectrum
    #--------------------------------------------------------------------------#
    if dbg>0: print("Plot spectrum")
    excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,fname="spectrum.png")
    for icalc in range(excit.ncalc):
        excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,jcalc=icalc,fname=f"spectrumEped_{icalc+1}.png")
    if conf.opt["Fit"].get("plot_result",False):
#        for iarea in range(excit.narea):
#            for icomp in range(excit.ncomp):
#                excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,jarea=iarea,jcomp=icomp)
        excit.plotPanels(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop)
    if dbg>0: print(" - done")

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
