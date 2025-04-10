#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# Authors: Ingo Schelter* & Rian Richter
#   *Correspondence to ingo.schelter@uni-bayreuth.de
#------------------------------------------------------------------------------#
# See errorHandler.py for usage information or call
# ./BTevalSpec.py -h
#------------------------------------------------------------------------------#
import numpy as np
import sys
import os.path
from scipy.optimize import curve_fit
import concurrent.futures
import getopt
import ruamel.yaml

#Import own & third-party modules
sys.path.append("../../src/")
import errorHandler as err
import config
import extern
import dipole
import excitations
import fit
import transdens as td
import cubetools as ct
try:
  import BTgrid
  import BTclust
  import BTcompact
  BTDFT=True
except:
  err.warn("No BTDFT modules => BTcompact format not available")
  BTDFT=False

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main(argv):

    #--------------------------------------------------------------------------#
    # Process command-line arguments
    #--------------------------------------------------------------------------#
    done    = False

    #--------------------------------------------------------------------------#
    # Get all command-line arguments and options
    #--------------------------------------------------------------------------#
    #try:
    opts, args = getopt.getopt(argv,"hv:f:",["help","verbose=","file=","minpw=","smooth=","window=","rmDC","wmax=","dw=","thin=","thres=","nsig=","nadd=","niter=","energy=","guess=","nofix","skip","single","range=","wref=","imag","exclude=","invert","signif","fitphase","reset","crit=","jcalc="])
    #except getopt.GetoptError:
    #    err.err(1,"In processing command line (e.g. missing argument or unknown option)!")

    #--------------------------------------------------------------------------#
    # Get command argument
    #--------------------------------------------------------------------------#
    try:
        cmd = args[0]
    except:
        err.usage()
        exit(0)

    #--------------------------------------------------------------------------#
    # Process general command line options
    #--------------------------------------------------------------------------#
    ifile       = "eval.yaml"
    verbose     = 1
    got_verbose = False
    got_ifile   = False
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            err.usage(cmd=cmd)
            sys.exit()
        elif opt in ("-v", "--verbose"):
            if got_verbose: err.err(1,"Multiple verbose arguments!",cmd=cmd)
            verbose = int(arg)
            got_verbose = True
        elif opt in ("-f", "--file"):
            if got_ifile: err.err(1,"Multiple file arguments!",cmd=cmd)
            ifile   = arg
            got_ifile = True

    #--------------------------------------------------------------------------#
    # Handle command "new"
    if cmd == "new":
        #----------------------------------------------------------------------#
        # Create configuration file 
        config.writeTemplate(ifile)
        print(f"")
        print(f"Please add dipole-moment and extern-profile files to configuration file {ifile}")
        print(f"")
        done=True

    #----------------------------------------------------------------------#
    # Read configuration from input file
    if verbose>0: print(f"Read configuration file {ifile}",end="")
    if not os.path.isfile(ifile):
        err.err.err(1,"No configuration file '{ifile}'. Select an existing file using '-f <file>' option or create a new configuration file using the 'new' command",cmd=cmd)
    conf = config.Config(ifile)
    if verbose>0: print(" - done")

    if not done:
        #----------------------------------------------------------------------#
        # Read dipole moments into dip[ncalc][narea] list:
        # In case of an external excitation with t_end>0:
        # Move the time frame so that it starts at t_end
        #----------------------------------------------------------------------#
        if verbose>0: print("Read dipole moment files",end="")
        dip = []
        for icalc in range(len(conf.dipfiles)):
            dip.append([])
            for iarea in range(len(conf.dipfiles[icalc])):
                dip[icalc].append(dipole.Dipole(conf.dipfiles[icalc][iarea],["x","y","z"]))
        if verbose>0: print(" - done")

    if cmd == "ft" or conf.opt["FT"].get("calc",False) and not done:
        #--------------------------------------------------------------------------#
        # Fourier transform
        minpw  = conf.opt["FT"].get("minpw" ,17   )
        window = conf.opt["FT"].get("window",0.   )
        smooth = conf.opt["FT"].get("smooth",0.   )
        rmDC   = conf.opt["FT"].get("rmDC"  ,False)
        got_minpw  = False
        got_window = False
        got_smooth = False
        for opt, arg in opts:
            if opt in ("--minpw"):
                if got_minpw: err.err(1,"Multiple minpw arguments!",cmd=cmd)
                minpw = int(arg)
                got_minpw = True
            elif opt in ("--window"):
                if got_window: err.err(1,"Multiple window arguments!",cmd=cmd)
                window = float(arg)
                got_window = True
            elif opt in ("--smooth"):
                if got_smooth: err.err(1,"Multiple smooth arguments!",cmd=cmd)
                smooth = float(arg)
                got_smooth = True
            elif opt in ("--rmDC"):
                rmDC=True

        if verbose>0: print("Calculate Fourier transform",end="")
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getFt(minpw=minpw,window=window,smooth=smooth,rmDC=rmDC))
                    dip[icalc][iarea].writeSpectra(what=["ft","pw"])
        if verbose>0: print(" - done")
        conf.opt["FT"]["calc"]   = False #Next time: read by default
        conf.opt["FT"]["minpw"]  = minpw
        conf.opt["FT"]["window"] = window
        conf.opt["FT"]["smooth"] = smooth
        conf.opt["FT"]["rmDC"]   = rmDC
        if cmd=="ft": done = True
    elif not done and not cmd=="pade":
        if verbose>0: print("Read Fourier transform",end="")
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["ft","pw"])
        if verbose>0: print(" - done")

    if cmd == "pade" or conf.opt["Pade"].get("calc",False) and not done: 
        #----------------------------------------------------------------------#
        # Pade approximation
        wmax   = conf.opt["Pade"].get("wmax"  ,0.5   )
        dw     = conf.opt["Pade"].get("dw"    ,1.0e-5)
        smooth = conf.opt["Pade"].get("smooth",0.    )
        thin   = conf.opt["Pade"].get("thin"  ,0     )
        got_wmax   = False
        got_dw     = False
        got_smooth = False
        got_thin   = False
        for opt, arg in opts:
            if opt in ("--wmax"):
                if got_wmax: err.err(1,"Multiple wmax arguments!",cmd=cmd)
                wmax = float(arg)
                got_wmax = True
            elif opt in ("--dw"):
                if got_dw: err.err(1,"Multiple dw arguments!",cmd=cmd)
                dw = float(arg)
                got_dw = True
            elif opt in ("--smooth"):
                if got_smooth: err.err(1,"Multiple smooth arguments!",cmd=cmd)
                smooth = float(arg)
                got_smooth = True
            elif opt in ("--thin"):
                if got_thin: err.err(1,"Multiple thin arguments!",cmd=cmd)
                thin = float(arg)
                got_thin = True
        if verbose>0: print("Calculate Pade approximation",end="")
        if not smooth > 0.: #Choose automatic smoothing
            smooth = float(np.log(100)/dip[0][0].tprop) # e^-s*T=0.01 <=> s=ln(100)/T
        with concurrent.futures.ThreadPoolExecutor() as executer:
            for icalc in range(len(dip)):
                for iarea in range(len(dip[icalc])):
                    executer.submit(dip[icalc][iarea].getPade(wmax,dw,thin=thin,smooth=smooth))
                    dip[icalc][iarea].writeSpectra(what=["pade"])
        conf.opt["Pade"]["calc"]   = False #Next time: read by default
        conf.opt["Pade"]["wmax"]   = wmax 
        conf.opt["Pade"]["dw"]     = dw
        conf.opt["Pade"]["smooth"] = smooth
        conf.opt["Pade"]["thin"]   = thin
        if verbose>0: print(" - done")
        if cmd=="pade": done = True
    elif not done:
        if verbose>0: print("Read Pade approximation",end="")
        for icalc in range(len(dip)):
            for iarea in range(len(dip[icalc])):
                dip[icalc][iarea].readSpectra(what=["pade"])
        if verbose>0: print(" - done")

    if not done:
        #----------------------------------------------------------------------#
        # Read and Fourier transform external excitation
        if verbose>0: print("Initialize extern potential",end="")
        extprof = conf.ext.get("profile","")
        if dip[0][0].ext!="boost" and extprof=="": err.err(1,"Extern profile file required for non-boost calculation",cmd=cmd)
        ext = extern.Extern(extprof,conf.ext.get("invertPhase",False),dip[0][0].efield,dip[0][0].text,[dip[icalc][0].epol for icalc in range(len(dip))])
        ext.write()
        if verbose>0: print(" - done")

        #----------------------------------------------------------------------#
        # Create excit structure
        if verbose>0: print("Initialize excitations",end="")
        excit = excitations.Excitations(ncalc=len(conf.dipfiles),narea=len(conf.dipfiles[0]),ncomp=3,ext=ext,exlist=conf.excit)
        if verbose>0: print(" - done")

    #--------------------------------------------------------------------------#
    # Handle commands that change the configuration file or plot
    if cmd == "rm" and not done:
        try:
            excit.remove(rmname=args[1].split(","),dbg=verbose)
        except:
            err.err(1,"No excitations given as second argument or wrong format",cmd=cmd)
        done=True
    if cmd == "fix" and not done:
        invert  = False
        for opt, arg in opts:
            if opt in ("--invert"): invert=True
        excit.fix(whichname=args[1].split(","),dbg=verbose,inverse=invert,permanent=True)
        done=True
    if cmd == "release" and not done:
        invert  = False
        for opt, arg in opts:
            if opt in ("--invert"): invert=True
        excit.release(whichname=args[1].split(","),dbg=verbose,inverse=invert,permanent=True)
        done=True

    if not done:
        #----------------------------------------------------------------------#
        # Create fit structure
        fitrange  = conf.opt["Fit"].get("range",[0.0,0.4])
        wref      = 0.
        imagonly  = conf.opt["Fit"].get("imagonly",dip[0][0].ext=="boost")
        got_range = False
        got_wref  = False
        for opt, arg in opts:
            if opt in ("--range"):
                if got_range: err.err(1,"Multiple range arguments!",cmd=cmd)
                fitrange = [float(x) for x in arg.split(",")]
                got_range = True
                if len(fitrange)!=2: err.err(1,"Fit range must contain 2 float values, e.g., '--range=0.1,0.4'",cmd=cmd)
            elif opt in ("--wref"):
                if got_wref: err.err(1,"Multiple wref arguments!",cmd=cmd)
                wref = float(arg)
                got_wref = True
            elif opt in ("--imag"):
                imagonly = True
        dfit = fit.Fit(dip,ext,excit,fitrange,wref,imagonly=imagonly)
        conf.opt["Fit"]["range"]    = fitrange
        conf.opt["Fit"]["imagonly"] = imagonly

    if cmd == "plot" and not done:
        #----------------------------------------------------------------------#
        # Plot objective function
        dfit.plotObjective()
        #----------------------------------------------------------------------#
        # Plot spectrum
        excit.plot(gamma=np.pi/dip[0][0].tprop,fname="spectrum.png")
        for icalc in range(excit.ncalc):
            excit.plot(gamma=np.pi/dip[0][0].tprop,jcalc=icalc,fname=f"spectrum_ampl_{icalc+1}.png")
#        for iarea in range(excit.narea):
#            for icomp in range(excit.ncomp):
#                excit.plot(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop,jarea=iarea,jcomp=icomp)
        excit.plotPanels(conf.opt["Fit"]["range"],dw=0.00001,gamma=np.pi/dip[0][0].tprop)
        #----------------------------------------------------------------------#
        # Write Latex table
        excit.latexTable()
        excit.gnuTable()
        done=True

    #--------------------------------------------------------------------------#
    # Handle commands that actually do something
    #if (cmd == "guess" or conf.opt["Fit"].get("guess")!="no" or any("--guess" in opt for opt in opts)) and not done: 
    if (cmd == "guess" or any("--guess" in opt for opt in opts)) and not done: 
        #--------------------------------------------------------------------------#
        # Create initial guess for the spectrum fit
        guesstype     = "ft" #conf.opt["Fit"].get("guess","pade")
        #if guesstype=="no": guesstype="pade" #Choose pade as default if "no" is given by the eval.yaml
        thres         = 0.05
        nsig          = 2.
        got_thres     = False
        got_guesstype = False
        got_nsig      = False
        for opt, arg in opts:
            if opt in ("--guess"):
                if got_guesstype: err.err(1,"Multiple guess arguments!",cmd=cmd)
                guesstype = arg
                got_guesstype = True
            elif opt in ("--thres"):
                if got_thres: err.err(1,"Multiple thres arguments!",cmd=cmd)
                thres = float(arg)
                got_thres = True
            elif opt in ("--nsig"):
                if got_nsig: err.err(1,"Multiple nsig arguments!",cmd=cmd)
                nsig = float(arg)
                got_nsig = True
        if verbose>0: print("Initial guess",end="")
        excit = dfit.newGuess(hf=thres,guesstype=guesstype,nsigma=nsig,dbg=verbose)
        if verbose>0: print(" - done")
        #conf.opt["Fit"]["guess"] = "no" #Next time: No new initial guess
        if cmd=="guess": done = True

        #----------------------------------------------------------------------#
        # Write Fit files
        if verbose>0: print("Write fit(-guess) files")
        excit.excitFiles(dip[0][0].tprop)
        dfit.writeFit(dbg=verbose)
        if verbose>0: print(" - done")

    if cmd == "add" and not done:
        #----------------------------------------------------------------------#
        # Handle command options
        nofix   = False
        nex     = 0
        nsig    = 2.
        energy  = []
        got_nex     = False
        got_nsig    = False
        got_energy  = False
        for opt, arg in opts:
            if opt in ("--energy"):
                if got_nex or got_nsig or got_energy: err.err(1,"Multiple nex/nsig/energy arguments (exclude each other)!",cmd=cmd)
                energy = [float(x) for x in arg.split(",")]
                got_energy = True
            elif opt in ("--nsig"):
                if got_nex or got_nsig or got_energy: err.err(1,"Multiple nex/nsig/energy arguments (exclude each other)!",cmd=cmd)
                nsig = float(arg)
                got_nsig = True
            elif opt in ("--nadd"):
                if got_nex or got_nsig or got_energy: err.err(1,"Multiple nex/nsig/energy arguments (exclude each other)!",cmd=cmd)
                nex  = int(arg)
                got_nex  = True
            elif opt in ("--nofix"):
                nofix = True
        #----------------------------------------------------------------------#
        # Fix existing excitations
        if not nofix:
            if verbose>0: print("  - Fix excitations")
            excit.fix()
        #----------------------------------------------------------------------#
        # Add new excitations
        if verbose>0: print("  - Add new excitations:")
        if got_energy:
            excit, nadd = dfit.addEx(excit,dbg=verbose,addEnergies=energy)
        elif got_nex:
            excit, nadd = dfit.addEx(excit,dbg=verbose,nadd=nex)
        else:
            excit, nadd = dfit.addEx(excit,dbg=verbose,nsigma=nsig)
        if verbose>0: print("Added "+str(nadd))
        done=True

    if cmd == "fit" and not done:
        #----------------------------------------------------------------------#
        # Fit
        niter     = 1
        crit      = 0.
        skipfirst = False
        signif    = False
        single    = False
        nsig      = 2.
        nadd      = 0
        reset     = False
        fitphase  = conf.opt["Fit"].get("fitphase",False)
        got_niter = False
        got_crit  = False
        got_range = False
        got_nsig  = False
        got_nadd  = False
        for opt, arg in opts:
            if opt in ("--niter"):
                if got_niter: err.err(1,"Multiple niter arguments!",cmd=cmd)
                niter = int(arg)
                if niter>0: nadd = 1000 #almost unlimited number of excitations can be added
                got_niter = True
            if opt in ("--crit"):
                if got_nadd or got_crit: err.err(1,"Multiple nadd/crit arguments (exclusive)!",cmd=cmd)
                crit = float(arg)
                nadd = 1000 #almost unlimited number of excitations can be added
                if crit>1. or crit<0.:  err.err(1,"Criterion must be within [0;1]",cmd=cmd)
                got_crit = True
            elif opt in ("--nadd"):
                if got_nadd or got_crit: err.err(1,"Multiple nadd/crit arguments (exclusive)!",cmd=cmd)
                nadd = int(arg)
                got_nadd = True
            elif opt in ("--nsig"):
                if got_nsig: err.err(1,"Multiple nsig arguments!",cmd=cmd)
                nsig = float(arg)
                got_nsig = True
            elif opt in ("--skipfirst"):
                skipfirst = True
            elif opt in ("--reset"):
                reset = True
            elif opt in ("--fitphase"):
                fitphase = True
            elif opt in ("--single"):
                single = True
            elif opt in ("--signif"):
                signif = True

        if verbose>0: print("Fitting:")
        excit, fiterr = dfit.fit(dbg=verbose,tol=crit,addex=nadd,niter=niter,skipfirst=skipfirst,allSignif=signif,nsigma=nsig,firstsingle=single,resetErange=reset,fitphase=fitphase)
        conf.opt["Fit"]["fiterr"]   = float(fiterr)
        conf.opt["Fit"]["range"]    = fitrange
        conf.opt["Fit"]["fitphase"] = fitphase
        if verbose>0: print(" - done")

        #----------------------------------------------------------------------#
        # Write Fit files
        if verbose>0: print("Write fit files")
        excit.excitFiles(dip[0][0].tprop)
        dfit.writeFit(dbg=verbose)
        if verbose>0: print(" - done")

        done = True

    if cmd == "decouple":
        #Check if excitations exist
        nex = len(excit.exlist)
        if nex==0: err.err(1,"No excitation for decoupling",cmd=cmd)
        #Check if transition densities and energies exist and number=number of excitations
        densname = conf.densft.get("densft",[])
        densen   = conf.densft.get("densen",[])
        if any(nex!=nn for nn in [len(densname),len(densen)]): err.err(1,"Number of densities and energies must match number of excitations",cmd=cmd)
        #Check if calculation index exists
        jcalc = conf.densft.get("jcalc",0)
        got_jcalc = False
        for opt, arg in opts:
            if opt in ("--jcalc"):
                if got_jcalc: err.err(1,"Multiple jcalc arguments!",cmd=cmd)
                jcalc = int(arg)
                if jcalc>excit.ncalc: err.err(1,"jcalc is too large (or negative)",jcalc,cmd=cmd)
                got_jcalc = True
        if excit.ncalc > 1 and not got_jcalc: err.err(1,"You need to specify the calculation from which you got the given FT densities via the --jcalc= option.",cmd=cmd)
        #Read transition densities
        densft   = []
        densmeta = []
        if isinstance(densname[0],ruamel.yaml.comments.CommentedSeq):
            if os.path.splitext(densname[0][0])[1] == ".cube":
                ftype = "cube"
            else:
                err.err(1,"If the FT-density entries are lists, cube files are expected [real.cube,imag.cube]",cmd=cmd)
        else:
            if os.path.splitext(densname[0])[1] == ".compact":
                ftype="compact"
            elif os.path.splitext(densname[0])[1] == ".cube" and dfit.imagonly:
                ftype="cube"
            else:
                err.err(1,"Expect either one complex-valued compact file or a list of real/imag cube file for each excitation.",cmd=cmd)
        for fname in densname:
            if ftype=="cube":
                #In case of cube files, densname contains a list with real- and imag part cube files
                if dfit.imagonly:
                    if densname[0] is list:
                        data, meta = ct.read_cube(fname[1]) #Reading the imaginary part is sufficient in this case
                    else:
                        data, meta = ct.read_cube(fname) #Having only an imag part cube is fine in this case
                    #data *= 1.j     #This does not work since it is tried to save the result in the original data, which is a real-valued array
                    data = 1.j*data #This creates a new complex-valued arrey
                else:
                    data, meta = ct.read_imcube(fname[0],fname[1])

            elif ftype=="compact":
                #In case of compact files, densname contains complex-valued compact files
                try: 
                  comp = BTcompact.BTcompact(fname)
                  data, meta = comp.toCube(comp.readVal())
                except:
                  err.err(1,"BTcompact format not supported by this version",cmd=cmd)
            densft  .append(data)
            densmeta.append(meta)
        #Generate Hw
        energy  = np.zeros(nex,dtype=float)
        for iex, ex in enumerate(excit.exlist):
            energy[iex] = ex.energy
        Hw = ext.getVal(energy)
        #Call transdecoupling
        transDens = td.decouple(densft,densen,excit,dip[jcalc][0].tprop,dip[jcalc][0].efield,dip[jcalc][0].epol,Hw,jcalc=jcalc,dbg=verbose,imagonly=dfit.imagonly)
        transDip  = -np.sqrt(2)*td.getDipole(np.real(transDens),meta["org"],[meta["xvec"][0],meta["yvec"][1],meta["zvec"][2]]) #sqrt(2) = elementary charge in Ry units
        excit.setTransdensDipoles(transDip)
        if verbose>0:
            print(f"Transition density | abs. norm | real norm | imag norm | sin^2(dAng) | dAbs")
            for iex, ex in enumerate(excit.exlist):
                #Transition densities should be real-valued -> Check
                anorm = np.linalg.norm(        transDens[iex] )
                rnorm = np.linalg.norm(np.real(transDens[iex]))
                inorm = np.linalg.norm(np.imag(transDens[iex]))
                #Check if transition densities reproduce fitted transition dipoles 
                delta_ang = 1.-(np.dot(transDip[iex],ex.dipole)/(np.linalg.norm(transDip[iex])*np.linalg.norm(ex.dipole)))**2 #sin(x)^2 = 1-cos(x)^2 in [0:1]
                delta_abs = (np.linalg.norm(transDip[iex])-np.linalg.norm(ex.dipole))/np.linalg.norm(ex.dipole)
                print(f"{ex.name:19}  {anorm:9.4f}   {(rnorm/anorm)**2*100:7.2f} %   {(inorm/anorm)**2*100:7.2f} %   {delta_ang:7.2f}     {delta_abs*100:7.2f} %")
        #Write transition densities
        for iex, ex in enumerate(excit.exlist):
            #Only write real-part (since the transition density should be real-valued)
            fname = f"{ex.name}_{ex.energy:.5f}_TransDens.cube"
            ct.write_cube(np.real(transDens[iex]),meta,fname,comment1=f"{ex.name} transition density at {ex.energy} Ry / {ex.energy*13.606} eV",comment2="")
            if verbose>1:
                #Write imag part for debugging
                fname = ex.name+"_TransDens_imag.cube"
                ct.write_cube(np.imag(transDens[iex]),meta,fname,comment1=f"{ex.name} transition density at {ex.energy} Ry / {ex.energy*13.606} eV",comment2="")
        #Add jcalc/ftype to conf
        conf.densft["jcalc"] = jcalc
        #Done
        done = True

    if not done: err.err(1,"Unknown command "+cmd,cmd=cmd)

    #--------------------------------------------------------------------------#
    # Update configuration file
    #--------------------------------------------------------------------------#
    if verbose>0: print("Update configuration file")
    try:
        conf.excit = [excit.exlist[i].todict() for i in range(len(excit.exlist))]
    except:
        conf.excit = []
    conf.write(ifile) 
    if verbose>0: print(" - done")

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main(sys.argv[1:])
