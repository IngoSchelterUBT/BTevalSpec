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
#
# Call
#
#  ./eval.py [--debug=<dbg>] [<cmd-opt>] cmd arg
# 
# Create a new eval.yaml file
#  ./eval.py [--ext=<extern profile>] new <dipole files>
#    Generates a new eval.yaml (or other) file based on the dipole-moment files
#    Automatically groups the dipole-moment files into calculations and areas.
# Fourier transform
#  ./eval.py [--minpw=<pw>] [--smooth=<smooth>] [--window=<window>] ft
#    Fourier transforms the dipole moment files, writes the transformations, and updates eval.yaml. 
# Pade approx
#  ./eval.py [--wmax=<wmax>] [--dw=<dw>] [--smooth=<smooth>] [--thin=<thin>] pade
#    Pade approximates the dipole moment files, writes the transformations, and updates eval.yaml.
# Make a new guess based on the pade approximation
#  ./eval.py [<pade opt>] [--thres=<thres>] guess <lb,rb>
#    Creates a new guess, sets the plot range, and updates eval.yaml.
#    Does a prior Pade approximation with the given <pade opt> options if not already done.
# Remove excitation without fit
#  ./eval.py rm <exlist>
#    Remove <exlist> from the list of excitations, where exlist is a comma-separated list of excitation labels, e.g., S1,S2,S3
# Add excitation without fit
#  ./eval.py [--label=<label>] [--nsig=<nsig>] add [<energy>]
#    Add new excitation with given label (optional) at given energy (optional).
#    Do no fit but guess phase and dipoles moments.
#    Add excitations that are larger than the mean value of maxima of ft-fit plus <nsig> times the standard deviation of maxima heights
# Reset energy range to the standard pi/T interval around the current excitations' energy values
#  ./eval.py reset
# Fit
#  ./eval.py [<ft opt>] [<guess opt>] [--skip] [--range=<lb,rb>] [--imag] fit [<nadd>] [<listOfExcitations>]
#    Fit current excitations (if not --skip), add and fit <nadd> excitations one after the other, and update eval.yaml.
#    --imag: only use imaginary part for fitting; automatically true for boost excitation.
#    <nadd> defaults to +0
#    <exs>  Comma-separated list of excitations to fit
#    Does a prior Fourier transform if not already done.
#    Does a prior Guess if not already done; requires a range option.
# Error
#  ./eval.py error
#    Compute error value
# Significance
#  ./eval.py significance [<listOfExcitations>]
#    Compute significances of all lines or the excitations given in the comma-separated list of <listOfExcitations>
# Plot
#  ./eval.py [--exclude=<listOfExcitations>] plot <listOfMeasures>
#    Plots <listOfMeasure> in {pade, ft, fit, err, spectrum} with excitations
#    --exclude=<listOfExcitations>: Excitations that are excluded from the (fitted) data
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

    #--------------------------------------------------------------------------#
    # Read configuration from eval.yaml
    #--------------------------------------------------------------------------#
    if not os.path.isfile(ifile):
        config.writeTemplate(ifile)
        err.err(1,'There was no '+ifile+' file, now a template is created!')
    conf  = config.Config(ifile)
    excit = excitations.Excitations(narea=len(conf.dipfiles[0]),ncomp=3,exlist=conf.excit)

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
        if not conf.opt["Pade"].get("smooth",0.) > 0.: #Choose automatic smoothing
            conf.opt["Pade"]["smooth"] = float(np.log(100)/dip[0][0].tprop) # e^-s*T=0.01 <=> s=ln(100)/T
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getPade(conf.opt["Pade"]["wmax"],conf.opt["Pade"]["dw"],thin=conf.opt["Pade"]["thin"],smooth=conf.opt["Pade"]["smooth"]))
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
    dfit = fit.Fit(dip,ext,excit,conf.opt["Fit"]["range"])
    if conf.opt["Fit"]["guess"]:
        dfit.newGuess(hf=conf.opt["Fit"]["guess_thres"],dbg=2)
        conf.opt["Fit"]["guess"] = False #Next time: No new initial guess

    #--------------------------------------------------------------------------#
    # Fit
    #--------------------------------------------------------------------------#
    if conf.opt["Fit"]["calc"]:
        excit, fiterr = dfit.fit(dbg=2,tol=conf.opt["Fit"]["relerr_crit"],maxex=conf.opt["Fit"]["max_excit"],skipfirst=conf.opt["Fit"].get("skipfirst",False),signif=conf.opt["Fit"].get("significances",False),nsigma=conf.opt["Fit"].get("nsigma",2.))
        conf.opt["Fit"]["skipfirst"] = True
        conf.opt["Fit"]["fiterr"]    = float(fiterr)

    #--------------------------------------------------------------------------#
    # Update configuration file
    #--------------------------------------------------------------------------#
    conf.excit = [excit.exlist[i].todict() for i in range(len(excit.exlist))]
    conf.write(ifile) 

    #--------------------------------------------------------------------------#
    # Write Fit files
    #--------------------------------------------------------------------------#
    if conf.opt["Fit"]["calc"]: dfit.writeFit()

    #--------------------------------------------------------------------------#
    # Plot spectrum
    #--------------------------------------------------------------------------#
    excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,fname="spectrum.png")
    if conf.opt["Fit"].get("plot_result",False):
#        for iarea in range(excit.narea):
#            for icomp in range(excit.ncomp):
#                excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,jarea=iarea,jcomp=icomp)
        excit.plotPanels(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop)

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
