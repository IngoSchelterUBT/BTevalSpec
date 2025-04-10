#==============================================================================#
# Fitting Class
#==============================================================================#
#intrinsic
import os
import sys
import numpy as np
from lmfit import Parameters, minimize, fit_report
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#own
import errorHandler as err
import dipole
import excitations
import extern
from mathtools import fspectrum, fspectrumDiff #, butter_lowpass_filter

class Fit:
    def __init__(self,dip,ext,excit,fitrange,wref,imagonly=False):
        self.fitrange = fitrange
        self.dip      = dip
        self.ext      = ext
        self.exttype  = self.dip[0][0].ext
        self.tprop    = self.dip[0][0].tprop
        self.fwhm     = np.pi/self.tprop
        self.ncalc    = len(self.dip)
        self.narea    = len(self.dip[0])
        self.ncomp    = len(self.dip[0][0].ft)
        self.excit    = excit
        self.freq     = np.array([dip[0][0].freq    [i] for i in range(len(dip[0][0].freq    )) if dip[0][0].freq    [i]>=self.fitrange[0] and dip[0][0].freq    [i]<=self.fitrange[1]])
        self.freqLoc  = np.copy(self.freq)
        self.dw       = (self.freq[-1]-self.freq[0])/(len(self.freq)-1)
        self.freqPade = np.array([dip[0][0].freqPade[i] for i in range(len(dip[0][0].freqPade)) if dip[0][0].freqPade[i]>=self.fitrange[0] and dip[0][0].freqPade[i]<=self.fitrange[1]])
        self.Nf       = len(self.freq)
        self.NfPade   = len(self.freqPade)
        self.imagonly = imagonly
        self.rc       = np.array([1]) if (self.dip[0][0].ext=="boost" or self.imagonly) else np.array([0,1]) #Use only imag(rc==1) part if boost excitation or real(rc==0) and imag(rc==1) part else
        self.nrc      = len(self.rc)
        self.ft       = []
        self.pw       = []
        self.pade     = []
        self.ftrc     = [] #Used for fit; contains only imag part for boost
        for icalc in range(self.ncalc):
            self.ft  .append([])
            self.pade.append([])
            self.pw  .append([])
            self.ftrc.append([])
            for iarea in range(self.narea):
                d = self.dip[icalc][iarea]
                self.ft  [icalc].append([])
                self.pade[icalc].append([])
                self.pw  [icalc].append([])
                self.ftrc[icalc].append([])
                for n in range(len(d.ft)):
                    self.ft  [icalc][iarea].append([])
                    self.pade[icalc][iarea].append([])
                    self.pw  [icalc][iarea].append(np.array([d.pw  [n][i] for i in range(len(d.pw  [n]   )) if d.freq    [i]>=self.fitrange[0] and d.freq    [i]<=self.fitrange[1]]))
                    self.ftrc[icalc][iarea].append([])
                    for irc in range(2):
                        self.ft  [icalc][iarea][n].append(np.array([d.ft  [n][        irc ][i] for i in range(len(d.ft  [n][        irc ])) if d.freq    [i]>=self.fitrange[0] and d.freq    [i]<=self.fitrange[1]]))
                        self.pade[icalc][iarea][n].append(np.array([d.pade[n][        irc ][i] for i in range(len(d.pade[n][        irc ])) if d.freqPade[i]>=self.fitrange[0] and d.freqPade[i]<=self.fitrange[1]]))
                    for irc in range(self.nrc):
                        self.ftrc[icalc][iarea][n].append(np.array([d.ft  [n][self.rc[irc]][i] for i in range(len(d.ft  [n][self.rc[irc]])) if d.freq    [i]>=self.fitrange[0] and d.freq    [i]<=self.fitrange[1]]))
        self.ft       = np.array(self.ft  )
        self.pw       = np.array(self.pw  )
        self.ftrc     = np.array(self.ftrc) #Only imag part if boost; used for fit
        self.ftrcLoc  = np.copy(self.ftrc)
        self.ftfit    = self.getFitFunc(self.excit)
        self.pwfit    = self.getPw(self.ftfit)
        self.ftfitrc  = self.getFitFunc(self.excit,self.rc) #only imag part if boost
        self.ftrcnorm = np.linalg.norm(self.ftrc.flatten()) #Use ftrc as reference
        self.ftrcnormLoc = np.linalg.norm(self.ftrcLoc.flatten()) #Use ftrc as reference
        self.fiterr   = self.getError(self.ftrc,self.ftfitrc,datnorm=self.ftrcnorm)
        self.wref     = wref
        self.scal     = self.getScaling(self.excit,dbg=0,wref=self.wref)
        self.breakMinimization = False
        self.runningError = 1000000.

    #--------------------------------------------------------------------------#
    # Get the fit function for a given set of frequencies and excitations
    #--------------------------------------------------------------------------#
    def getFitFunc(self,excit,rc=[0,1],freq=None,subtractfrom=None):
        # Setup new array
        if not isinstance(freq,np.ndarray): freq = self.freq
        nrc     = len(rc)
        nf      = len(freq)
        # Get relevant excitation measures
        T       = self.tprop
        Ef      = self.dip[0][0].efield
        nex     = len(excit.exlist)
        energy  = np.zeros(nex,dtype=float)
        phase   = np.zeros(nex,dtype=float)
        tmod    = np.zeros(nex,dtype=float)
        for iex, ex in enumerate(excit.exlist):
            energy[iex]        = ex.energy
            phase [iex]        = ex.phase
            tmod  [iex]        = ex.tmod
        try:
            Hw = self.ext.getVal(energy)
        except:
            raise
        ampl = np.zeros((self.ncalc,self.narea,self.ncomp,nex),dtype=float)
        for icalc in range(self.ncalc):
            Ep = self.dip[icalc][0].epol #Polarization can differ between calculations
            for iarea in range(self.narea):
                for n in range(self.ncomp):
                    for iex, ex in enumerate(excit.exlist):
                        #ampl[icalc][iarea][n][iex] = -Ef*np.abs(np.dot(Ep,ex.dipole))*np.abs(Hw[iex])*ex.dipoles[iarea][n]
                        ampl[icalc][iarea][n][iex] = -Ef*       np.dot(Ep,ex.dipole) *np.abs(Hw[iex])*ex.dipoles[iarea][n] #*1/hbar, which is one in Ry a.u.
        if isinstance(subtractfrom,np.ndarray):
            return fspectrumDiff(self.ncalc,self.narea,self.ncomp,np.array(rc),T,freq,energy,phase,tmod,ampl,subtractfrom)
        else:
            return fspectrum(self.ncalc,self.narea,self.ncomp,np.array(rc),T,freq,energy,phase,tmod,ampl)

    #--------------------------------------------------------------------------#
    # Compute Power spectrum from FT
    #--------------------------------------------------------------------------#
    def getPw(self,ft):
        pw = []
        for icalc in range(self.ncalc):
            pw.append([])
            for iarea in range(self.narea):
                pw[icalc].append([])
                for icomp in range(self.ncomp):
                    pw[icalc][iarea].append(np.zeros(self.Nf,dtype=float))
                    for i in range(self.Nf):
                        pw[icalc][iarea][icomp][i] = 2.*(ft[icalc][iarea][icomp][0][i]**2 + ft[icalc][iarea][icomp][1][i]**2)
        return np.array(pw)

    #--------------------------------------------------------------------------#
    # Errors dat vs. fit functions
    #--------------------------------------------------------------------------#
    def getError(self,dat,fit,datnorm=None):
        #fdat = np.array(dat).flatten()
        #ffit = np.array(fit).flatten()
        fdat = np.array(dat).ravel()
        ffit = np.array(fit).ravel()
        if isinstance(datnorm,float):
            return np.linalg.norm(np.subtract(fdat,ffit))/datnorm
        else:
            return np.linalg.norm(np.subtract(fdat,ffit))/np.linalg.norm(fdat)

    #--------------------------------------------------------------------------#
    # Make an initial guess based on the sum of the Pade power spectra and the FT spectra
    #--------------------------------------------------------------------------#
    def newGuess(self,guesstype="ft",hf=0.05,dbg=0,nsigma=2.):
        # Init empty excitation object and arrays
        self.excit = excitations.Excitations(ncalc=self.ncalc,narea=self.narea,ncomp=self.ncomp,ext=self.ext)
        # Pade-based guess
        if guesstype in ["pade"]:
            fdat       = self.pade
            f          = np.zeros(len(self.freqPade),dtype=float)
            # Sum up squared Pade spectra
            for icalc in range(self.ncalc):
                for iarea in range(self.narea):
                    for n in range(self.ncomp):
                        for irc in range(2):
                            for i in range(self.NfPade):
                                f[i] += fdat[icalc][iarea][n][irc][i]**2
            # Search most prominent peaks
            piT = np.pi/self.tprop
            peaks = self.findPeaks(self.freqPade,np.sqrt(f),dbg=dbg,hf=hf,dw=piT,minwidth=piT)
            # For each peak
            for ipeak in range(len(peaks)):
                #Get energy and guess phase and dipoles based on the FT
                energy  = self.freqPade[peaks[ipeak]]
                phase, dipoles = self.guessExcit(energy)
                #Add excitation to list
                self.excit.add(energy=energy,phase=phase,dipoles=dipoles,erange=piT)
        elif guesstype in ["ft",True,"true"]:
            self.excit, nadd = self.addEx(self.excit,dbg=dbg,nsigma=nsigma)
            if nadd==0:
                self.excit, nadd = self.addEx(self.excit,dbg=dbg,singleMax=True,nsigma=nsigma)
        elif guesstype!="no":
            err.err(1,"Unknown guess type "+str(guesstype))
        # Build fit-spectrum
        self.update()
        return self.excit

    #--------------------------------------------------------------------------#
    # Guess phase and dipole moments of an excitation at an approx. energy
    #--------------------------------------------------------------------------#
    def guessExcit(self,energy):
        # Get the excitation strength/phase at that energy
        Hw    = self.ext.getVal([energy])[0]
        phase = np.angle(Hw)%(2.*np.pi)
        # Get the closest index in the FT array
        idx = int(np.rint((energy-self.freq[0])/self.dw))
        # Find the calculation for which the excitation's peak is largest
        ampl        = [[[(self.ft[icalc][iarea][icomp][0][idx]-self.ftfit[icalc][iarea][icomp][0][idx])+1.0j*(self.ft[icalc][iarea][icomp][1][idx]-self.ftfit[icalc][iarea][icomp][1][idx]) for icomp in range(self.ncomp)] for iarea in range(self.narea)] for icalc in range(self.ncalc)]
        jcalc       = np.argmax([np.linalg.norm([sum([ampl[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)]) for icalc in range(self.ncalc)]) #Use this calculation for the dipole-moment guess
        # Approx. dipole directions
        T           = self.tprop
        Ef          = self.dip[    0][0].efield
        Ep          = self.dip[jcalc][0].epol
        heightArea  = [[np.imag(np.exp(-1.0j*phase)*     ampl[jcalc][iarea][icomp]) for iarea in range(self.narea)]   for icomp in range(self.ncomp)] #Rotate the height by e^-i*phi and take its imag part to get the line heights (with sign)
        height      = [ np.imag(np.exp(-1.0j*phase)*sum([ampl[jcalc][iarea][icomp]  for iarea in range(self.narea)])) for icomp in range(self.ncomp)]
        dipdir      = height/np.linalg.norm(height) # Direction of the dipole moment
        eped        = np.dot(dipdir,Ep) # Scalar product e_pol.e_mu
        #Prevent overestimating the dipole moment if the latter is almost orthogonal to the ext. field polarization:
        # Alternative 1: Hard filter
        #if eped < 0.1: eped=1.     
        # Alternative 2: use sqrt(eped) instead of eped. This leads to an systematic underestimation of the dipole moment, especially the closer the angle between excitation and dipole is to 90°
        #eped        = np.sqrt(abs(eped)) 
        # -> Use this to scale dipoles, not dipabs below

        # Approx. dipoles
        #dipabs  = np.sqrt(     np.linalg.norm( height                                                                                   )/(T*Ef*np.abs(eped)*np.abs(Hw))) #Absolute value of the total dipole moment
        #Use the following instead since the total transition dipole may vanish for excitations without dipole symmetry
        dipabs  = np.sqrt(sum([np.linalg.norm([heightArea[icomp][iarea] for icomp in range(self.ncomp)]) for iarea in range(self.narea)])/(T*Ef*np.abs(eped)*np.abs(Hw))) #Absolute value of the total dipole moment
        #Todo: If |a_j| << sum_l |a_jl| then use second option, else, use first option
        dipoles = []
        for iarea in range(self.narea):
            dipoles.append([])
            for icomp in range(self.ncomp):
                #dipoles[iarea].append(heightArea[icomp][iarea]/(T*Ef*eped/np.sqrt(np.abs(eped))*np.abs(Hw)*dipabs)) # hightarea[icomp] = T*Ef*eped*abs(Hw)*abs(mu)*mu[icomp] -> solve for mu[icomp]
                #Line above: sqrt(eped) is used to underestimate the dipoles more the less it is aligned with the ext. excitation polarization
                #Here the version without 1/sqrt(abs(eped))
                dipoles[iarea].append(heightArea[icomp][iarea]/(T*Ef*eped*np.abs(Hw)*dipabs)) # hightarea[icomp] = T*Ef*eped*abs(Hw)*abs(mu)*mu[icomp] -> solve for mu[icomp]

        return phase, dipoles
            
    #--------------------------------------------------------------------------#
    # Search for the highest peak in the given function
    #--------------------------------------------------------------------------#
    def findPeaks(self,freq,f,dbg=0,hf=0.,dw=None,minwidth=None,latex=False):
        pos, prop = find_peaks(f,height=hf*np.amax(f),width=minwidth/self.dw,distance=dw/self.dw)
      
        if dbg>1:
            print("Pade peak properties:")
            for i in range(len(pos)):
                print("  Peak ", i)
                print("    energy: ",pos[i])
                for key, val in prop.items():
                    print("    ",key,val)
        if dbg>0:
            if latex:
                plt.rcParams["text.usetex"] = True
                plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
                plt.ylabel(r"$\|\boldsymbol{\Delta}_{\mathrm{pade}}\|^\prime$ [a.u.]")
            else:
                plt.ylabel("Add-line objective [a.u.]")
            plt.xlabel("Energy [Ry]")
            plt.plot(freq,f,color="#36454F")
            plt.plot(freq[pos],f[pos],"x",markersize=10,color="#DAA520")
            plt.savefig("padeObjective.png",dpi=300)
            plt.show()
            if latex:
                plt.rcParams["text.latex.preamble"] = r""
                plt.rcParams["text.usetex"] = False

        return pos

    #--------------------------------------------------------------------------#
    # Fit objective
    # nfree is not used here but required since fitInter needs it
    #--------------------------------------------------------------------------#
    def fitObj(self,params,excit=None,noPhase=False,ext=None,dbg=0,breakmod=None,nfree=0):
        if not isinstance(excit,excitations.Excitations): excit = self.excit
        if not isinstance(ext  ,     extern.Extern     ): ext   = self.ext
        excit.updateFromParam(params,noPhase,ext,errors=False) #fit errors are not known, yet
        if self.breakMinimization==True:
            obj  = np.zeros(len(ffit),dtype=float)
        else:
            try:
                obj = self.getFitFunc(excit,self.rc,freq=self.freqLoc,subtractfrom=self.ftrcLoc)#.ravel() (array need not be flat as it seems)
            except:
                raise
        return obj

    #--------------------------------------------------------------------------#
    # Fit inter - is called in every minimizer iteration
    #--------------------------------------------------------------------------#
    def fitInter(self,params,iter,resid,excit=None,noPhase=False,ext=None,dbg=0,breakmod=5,nfree=0):
        if self.breakMinimization: return
        currentError = np.sum(np.abs(resid))/self.ftrcnormLoc
        if iter==-1:
            self.bestParams   = params
            self.bestError    = currentError
            self.runningError = currentError
        else:
            if currentError<self.bestError:
                self.bestError  = currentError
                self.bestParams = params
        if dbg>2: print("DEBUG: Objective Norm: ", iter, currentError)
        #if iter%(breakmod*len(params.valuesdict()))==0 and iter>0:
        #    if abs(currentError/self.runningError)>0.99:
        #        print("WARNING: Minimization stuck; abort and fall back to best parameter set",file=sys.stderr)
        #        self.breakMinimization = True
        #    if iter>-1: self.runningError = currentError

    #--------------------------------------------------------------------------#
    # Fit Atomic (and update excit)
    #--------------------------------------------------------------------------#
    def fitAtomic(self,excit,dbg=0,breakmod=5,fitrange=None,noPhase=False):
        noPhase = noPhase or self.exttype=="boost"
        params = excit.toParams(self.ext,noPhase=noPhase,noTmod=True) #fix phase if boost and fix time always
        #params = excit.toParams(self.ext) #free phases and time modifier
        nfree  = excit.countFree()
        if isinstance(fitrange,list):
            self.freqLoc = np.array([self.freq[i] for i in range(len(self.freq)) if self.freq[i]>=fitrange[0] and self.freq[i]<=self.fitrange[1]])
            nf = len(self.freqLoc)
            ftrcLoc = []
            for icalc in range(self.ncalc):
                ftrcLoc.append([])
                for iarea in range(self.narea):
                    ftrcLoc[icalc].append([])
                    for icomp in range(self.ncomp):
                        ftrcLoc[icalc][iarea].append([])
                        for irc in range(self.nrc):
                            ftrcLoc[icalc][iarea][icomp].append([self.ftrc[icalc][iarea][icomp][irc][i] for i in range(len(self.freq)) if self.freq[i]>=fitrange[0] and self.freq[i]<=self.fitrange[1]])
            self.ftrcLoc = np.array(ftrcLoc)
        else:
            self.freqLoc = np.copy(self.freq)
            self.ftrcLoc = np.copy(self.ftrc)
        #self.ftrcnormLoc = np.linalg.norm(self.ftrcLoc.flatten()) #Use ftrc as reference
        self.ftrcnormLoc = np.linalg.norm(self.ftrcLoc.ravel()) #Use ftrc as reference

        try:
            #fitres = minimize(self.fitObj,params,iter_cb=self.fitInter,kws={"excit":excit,"noPhase":noPhase,"ext":self.ext,"dbg":dbg,"breakmod":breakmod,"nfree":nfree})
            fitres = minimize(self.fitObj,params,kws={"excit":excit,"noPhase":noPhase,"ext":self.ext,"dbg":dbg,"breakmod":breakmod,"nfree":nfree})
        except:
            raise #self.excit is not updated in this case

        if self.breakMinimization: #If minimization stuck and was force-aborted, reset the parameters to the best ones reached
            params = self.bestParams
            self.breakMinimization = False
        if dbg>1: print(fit_report(fitres))
        excit.updateFromParam(fitres.params,noPhase,self.ext)

    #--------------------------------------------------------------------------#
    # Return error-scaling function
    #--------------------------------------------------------------------------#
    def getScaling(self,excit,dbg=0,wref=1.):
        T      = self.tprop
        T2     = T*T
        scal   = np.full((self.ncalc,self.narea,self.ncomp,self.Nf),1.) #Fill with 1
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    for i in range(self.Nf):
                        w = self.freq[i]
                        for iex, ex in enumerate(excit.exlist):
                            dw   = w-ex.energy
                            dw2  = dw*dw
                            eped = abs(np.dot(ex.dipole/np.linalg.norm(ex.dipole),self.dip[icalc][0].epol/np.linalg.norm(self.dip[icalc][0].epol)))
                            #eped = 1.
                            f    = ex.strength
                            nmu  = np.abs(ex.dipoles[iarea][icomp]/np.linalg.norm(ex.dipole))
                            scal[icalc][iarea][icomp][i] += wref*T*f*eped*nmu/np.sqrt(T2*dw2+1)
        return scal
    
    #--------------------------------------------------------------------------#
    # Compute objective function for adding new exictations
    # (obj0 is a version without error-compensating scaling)
    #--------------------------------------------------------------------------#
    def addExObj(self,fdat,ffit,scal):
        obj        = np.zeros(self.Nf)
        obj0       = np.zeros(self.Nf)
        objByArea  = np.zeros((self.narea,self.Nf))
        obj0ByArea = np.zeros((self.narea,self.Nf))
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    for irc in range(self.nrc):
                        obj                 += (np.abs(np.subtract(fdat[icalc][iarea][icomp][irc],ffit[icalc][iarea][icomp][irc]))/scal[icalc][iarea][icomp])**2
                        obj0                += (np.abs(np.subtract(fdat[icalc][iarea][icomp][irc],ffit[icalc][iarea][icomp][irc])))**2
                        objByArea [iarea,:] += (np.abs(np.subtract(fdat[icalc][iarea][icomp][irc],ffit[icalc][iarea][icomp][irc]))/scal[icalc][iarea][icomp])**2
                        obj0ByArea[iarea,:] += (np.abs(np.subtract(fdat[icalc][iarea][icomp][irc],ffit[icalc][iarea][icomp][irc])))**2
        return np.sqrt(obj), np.sqrt(obj0), np.sqrt(objByArea), np.sqrt(obj0ByArea)

    #--------------------------------------------------------------------------#
    # Add new excitation
    #--------------------------------------------------------------------------#
    def addEx(self,excit,dbg=0,addEnergies=[],nadd=0,nsigma=2.,maxadd=-1,latex=False):
        piT    = np.pi/self.tprop
        excit0 = excit.copy()
        fdat   = self.ftrc   
        ffit   = self.ftfitrc
        scal   = self.scal
        if maxadd>=0: nadd = min(nadd,maxadd)

        obj, obj0, objByArea, obj0ByArea = self.addExObj(fdat,ffit,scal)

        maxen  = []
        maxhei = []
        if len(addEnergies)>0:
            pos       = [min(range(len(self.freq)), key = lambda i: abs(self.freq[i]-en)) for en in addEnergies]
            maxen     = [self.freq[i] for i in pos]
            maxhei    = [obj      [i] for i in pos]
        else:
            pos, prop = find_peaks(obj)

            energies = []
            heights  = []
            for i in pos:
                energies.append(self.freq[i])
                heights .append(     obj [i])

            if nadd>0:
                idx = sorted(range(len(heights)), key = lambda i: heights[i])[-nadd:]
                maxen  = [energies[i] for i in idx]
                maxhei = [heights [i] for i in idx]
            else:
                minHeight  = np.amin(heights)
                meanHeight = np.mean(heights) #Mean value
                stdHeight  = np.std (heights) #Standard deviation
                baseHeight = meanHeight + nsigma*stdHeight
                maxen  = [energies[i] for i in range(len(energies)) if heights[i]>baseHeight]
                maxhei = [heights [i] for i in range(len(heights )) if heights[i]>baseHeight]
                if maxadd>=0 and len(maxen)>maxadd:
                    maxen  = [x for _,x in sorted(zip(maxhei,maxen ))][-maxadd:]
                    maxhei = [x for _,x in sorted(zip(maxhei,maxhei))][-maxadd:]

        for en in maxen:
            phase, dipoles = self.guessExcit(en)
            excit0.add(energy=en,phase=phase,dipoles=dipoles,erange=piT)

        if dbg>0:
            # Write data files
            head = 'Energy (Ry) | unscaled addExObj | scaled addExObj '
            np.savetxt("addExObj.dat",np.column_stack((self.freq,obj0,obj)),header=head)
            for iarea in range(self.narea):
                np.savetxt(f"addExObj_area{str(iarea+1).zfill(2)}.dat",np.column_stack((self.freq,obj0ByArea[iarea],objByArea[iarea])),header=head)

            # Plot
            # Objective by area
            for iarea in range(self.narea):
                fig, ax = plt.subplots()
                if latex:
                    plt.rcParams["text.usetex"] = True
                    plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
                    plt.ylabel(r"$\|\boldsymbol{\Delta}_s\|^\prime$ [a.u.]")
                else:
                    plt.ylabel(r"Add-line objective [a.u.]")
                ax.set_title(f"Add-line objective: Area {iarea:d}")
                ax.plot(self.freq,obj0ByArea[iarea,:],label="without scaling",color="#C0C0C3")
                ax.plot(self.freq,objByArea [iarea,:],label="with    scaling",color="#36454F")
                plt.xlabel("Energy [Ry]")
                plt.legend(loc="upper right")
                fig.savefig(f"addLineObjective_area{str(iarea+1).zfill(2)}.png",dpi=300)
            # Total objective
            fig, ax = plt.subplots()
            if latex:
                plt.rcParams["text.usetex"] = True
                plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
                plt.ylabel(r"$\|\boldsymbol{\Delta}_s\|^\prime$ [a.u.]")
            else:
                plt.ylabel(r"Add-line objective [a.u.]")
            try:
                ax.axhline(y=meanHeight                 ,color="#B22222",linestyle="-",label=r"\overline{h}")
                ax.axhline(y=meanHeight+nsigma*stdHeight,color="#B22222",linestyle=":",label=r"\overline{h}+2 \Delta h")
            except:
                pass
            ax.plot(maxen,maxhei,"x",markersize=10,color="#DAA520")
            ax.plot(self.freq,obj0,label="without scaling",color="#C0C0C3")
            ax.plot(self.freq,obj ,label="with    scaling",color="#36454F")
            plt.xlabel("Energy [Ry]")
            plt.legend(loc="upper right")
            fig.savefig("addLineObjective.png",dpi=300)
            # Show all
            plt.show()

#!            plt.plot(sorted(heights,reverse=True),"x",markersize=10)
#!            plt.axhline(y=meanHeight                      ,color="r",linestyle="-")
#!            plt.axhline(y=meanHeight+nsigma*stdHeight,color="r",linestyle=":")
#!            plt.savefig("heights.png")
#!            plt.show()
# Set different ylimits
#            ax = plt.gca()
#            ylim = ax.get_ylim()
#            ax.set_ylim([ylim[0]*0.05, ylim[1]*0.05])

            if latex:
                plt.rcParams["text.latex.preamble"] = r""
                plt.rcParams["text.usetex"] = False

        if dbg>0:
            excit0.print(long=dbg>1)
        return excit0, len(maxen)

    #--------------------------------------------------------------------------#
    # Fit wrapper
    #--------------------------------------------------------------------------#
    def fit(self,dbg=0,addex=0,tol=0.0,niter=0,skipfirst=False,allSignif=False,nsigma=2.,firstsingle=False,resetErange=False,fitphase=True):
        #------------------------------------------------------------------#
        # Init aux
        maxex = len(self.excit.exlist) + addex

        #------------------------------------------------------------------#
        # Report initial state
        if dbg>0:
            print("  - Initial state")
            self.reportFit(dbg=dbg)

        #----------------------------------------------------------------------#
        # Reset excitation-specific energy fit range
        if resetErange:
            if dbg>0: print("  - Reset energy range")
            self.excit.resetErange(self.fwhm)

        #----------------------------------------------------------------------#
        # Fit guess or input excitations
        if (self.exttype=="boost"): self.excit.zeroPhase()
        self.update(dbg=dbg)
        if not skipfirst:
            if firstsingle: #Fit single excitations
                if dbg>0: print("  - Fit single excitations")
                try:
                    #i) Sort excitations by effective size in the spectrum (strength*eped)
                    iexs = np.argsort([np.sqrt(sum([self.excit.exlist[iex].strengthEped[icalc]**2 for icalc in range(self.ncalc)])) for iex in range(len(self.excit.exlist))])
                    #ii) Fit single excitations one after the other
                    for iex in np.flip(iexs):              #Go through excitations from large to small
                        if dbg>0: print("    excit",iex)
                        self.excit.fix()                   #Fix all
                        self.excit.release(whichidx=[iex])    #Release single excitation iex
                        self.excit.exlist[iex].dipoles = np.full(self.excit.exlist[iex].dipoles.shape,0.001) #Set dipole moment of excitation to small value
                        self.fitAtomic(self.excit,dbg=dbg,noPhase=not fitphase) #Fix single excitation
                        self.excit.release()               #Release all excitations
                except:
                    raise #self.excit is not updated in this case
            else: #Fit all excitations at once
                if dbg>0: print("  - Fit excitations collectively")
                self.fitAtomic(self.excit,dbg=dbg,noPhase=not fitphase)   # Fit existing excitations
                self.update(dbg=dbg)      # Update some self.components
            self.reportFit(dbg=dbg)
        nex = len(self.excit.exlist)

        #----------------------------------------------------------------------#
        # Add and fit new excitations
        jiter = 0
        while len(self.excit.exlist)<maxex and jiter<niter:
            jiter += 1

            #------------------------------------------------------------------#
            # Add excitations
            if dbg>0: print("  - Fix excitations")
            self.excit.fix()                         # Temporarily fix all existing excitations
            if dbg>0: print("  - Add new excitations:")
            self.excit, nadd = self.addEx(self.excit,dbg=dbg,nsigma=nsigma,maxadd=maxex-len(self.excit.exlist)) # Add new excitations (also return this excitation)
            if dbg>0: print("Added "+str(nadd))
            if nadd==0:
                if dbg>0: print("  - Add single largest line")
                self.excit, nadd = self.addEx(self.excit,dbg=dbg,nadd=1) # Add single largest peak as excitation
                if dbg>0: print("Added "+str(nadd))

            #------------------------------------------------------------------#
            # Fit all new ones (all at once)
            try:
                if dbg>0: print("  - Fit new excitations",end="")
                self.fitAtomic(self.excit,dbg=dbg,noPhase=not fitphase) # Fit new excitations alone
                self.update(dbg=dbg)               # Update some self.components
                self.excit.release()               # Release all temporarily fixed excitations
                if dbg>0: print(" - done")
            except:
                raise #self.excit is not updated in this case

            #------------------------------------------------------------------#
            # Fit single excitations that still have zero error (those were not fitted properly, maybe due to their low weight)
            try:
                if dbg>0: print("  - Fit single new excitations that have zero error",end="")
                #i) Sort excitations by effective size in the spectrum (strength*eped)
                iexs = np.argsort([np.sqrt(sum([self.excit.exlist[iex].strengthEped[icalc]**2 for icalc in range(self.ncalc)])) for iex in range(len(self.excit.exlist))])
                #ii) Fit single excitations one after the other
                for iex in np.flip(iexs):              #Go through excitations from large to small
                    if self.excit.exlist[iex].strengthErr > 0.: continue
                    if dbg>0: print(iex,end="")
                    self.excit.fix()                   #Fix all
                    self.excit.release(whichidx=[iex])    #Release single excitation iex
                    self.excit.exlist[iex].dipoles = np.full(self.excit.exlist[iex].dipoles.shape,0.001) #Set dipole moment of excitation to small value
                    self.fitAtomic(self.excit,dbg=dbg,noPhase=not fitphase) #Fix single excitation
                    self.excit.release()               #Release all excitations
                if dbg>0: print(" - done")
            except:
                raise #self.excit is not updated in this case

#            #------------------------------------------------------------------#
#            # Attention: This may lead to an infinity loop since excitations are added in the next iteration again
#            # Exclude excitations Error>strength
#            rmidx = []
#            for iex, ex in enumerate(self.excit.exlist):
#                if abs(ex.strengthErr/ex.strength)>1.: rmidx.append(iex)
#            self.excit.remove(rmidx=rmidx)

            #------------------------------------------------------------------#
            # Fit all excitations
            try:
                if dbg>0: print("  - Fit all excitations",end="")
                self.fitAtomic(self.excit,dbg=dbg,noPhase=not fitphase)   # Fit all non-permanently fixed excitations
                self.update(dbg=dbg)      # Update some self.components
                if dbg>0: print(" - done")
            except:
                raise #self.excit is not updated in this case

            #------------------------------------------------------------------#
            # Report and check exit condition
            self.reportFit(dbg=dbg)
            if dbg>2: self.plotFitDebug(it)
            if self.fiterr<tol:
                if dbg>0: print("  - Reached tolerance -> exit")
                break
            
        #----------------------------------------------------------------------#
        # Simulate adding new excitations
        if dbg>0: print("  - Simulate adding new excitations:")
        extmp, nadd = self.addEx(self.excit,dbg=dbg,nsigma=nsigma)
        if dbg>0: print("Added "+str(nadd))
  
        #----------------------------------------------------------------------#
        # Compute significances (and update self.excit)
        self.update(dbg=dbg)      # Call update again if nothing has been done here
        if dbg>0: print("  - Calculate significances",end="")
        self.setSignificances(allSignif=allSignif,noPhase=not fitphase,dbg=dbg)
        if dbg>0: print(" - done")
        return self.excit, self.fiterr

    #--------------------------------------------------------------------------#
    # Update ftfit and fiterr based on self.excit and report
    #--------------------------------------------------------------------------#
    def update(self,dbg=0):
        self.ftfit   = self.getFitFunc(self.excit)
        self.pwfit   = self.getPw     (self.ftfit)
        self.ftfitrc = self.getFitFunc(self.excit,self.rc)
        self.fiterr  = self.getError  (self.ftrc,self.ftfitrc,datnorm=self.ftrcnorm)
        self.scal    = self.getScaling(self.excit,dbg=dbg,wref=self.wref)

    #--------------------------------------------------------------------------#
    # Report fit
    #--------------------------------------------------------------------------#
    def reportFit(self,dbg=0):
        print("--------------------------------")
        print(f"Nex {len(self.excit.exlist):3d} | Error {self.fiterr:7.4f}")
        if dbg>0:
            self.excit.print(long=dbg>1)
        print("--------------------------------")

    #--------------------------------------------------------------------------#
    # Compute excitation significances
    #--------------------------------------------------------------------------#
    def setSignificances(self,allSignif=False,noPhase=False,dbg=0):
        sN = self.getError(self.ftrc,self.ftfitrc,datnorm=self.ftrcnorm)
        for iex, ex in enumerate(self.excit.exlist):
            if dbg>0: print("    "+ex.name)
            self.setSignificance(self.excit,iex,sN=sN,rc=self.rc,allSignif=allSignif,noPhase=noPhase)

    #--------------------------------------------------------------------------#
    # Compute single significance
    #--------------------------------------------------------------------------#
    def setSignificance(self,exsN,iex0,sN=0.,rc=[0,1],allSignif=False,noPhase=False):
        piT = np.pi/self.tprop

        if allSignif:
            exs = exsN.copy()                                                 #Get a copy of the excitations object
            exs.restrict(erange=piT,neigh=True)                               #Restrict range of energy parameter before removing the excitation to ensure that a removed large excitation is not filled by a smaller one (leading to an erroneous low significance)
            exs.fixErange([exs.exlist[iex0].energy-2.*piT,exs.exlist[iex0].energy+2.*piT],inverse=True) #Fix excitations that are not closed than 2pi/T to the removed excitation
            fitrange = [exs.exlist[iex0].energy-4.*piT,exs.exlist[iex0].energy+4.*piT]                  #Restrict fit range to the removed excitation's energy +/- 4pi/T
            exs.remove(rmidx=[iex0])                                          #Remove the excitation in question
            
            if not sN>0.:
                ffitN = self.getFitFunc(exsN,rc=rc) #Complete fit function
            ffit0 = self.getFitFunc(exs ,rc=rc)     #Fit function with ex0 removed
            self.fitAtomic(exs,fitrange=fitrange,noPhase=noPhase)   #Fit exs new (no premature break out)
            ffit  = self.getFitFunc(exs ,rc=rc)     #Fit function with ex0 removed and re-fit

            exs = exsN.copy() #Reset exs structure
            exs.exlist[iex0].dipoles = np.full(exs.exlist[iex0].dipoles.shape,0.001) #Set dipole moment of excitation to small value
            exs.fix() #Fix all excitations
            exs.release(whichidx=[iex0]) #Release excitation in question
            self.fitAtomic(exs,fitrange=fitrange,noPhase=noPhase)   #Fit exs new (no premature break out)
            #ffit1 = self.getFitFunc(exs ,rc=rc)     #Fit function with ex0 alone fitted and all others fixed

            if not sN>0.:
                sN    = self.getError(self.ftrc,ffitN,datnorm=self.ftrcnorm) #Fit error of excitations exsN
            s0        = self.getError(self.ftrc,ffit0,datnorm=self.ftrcnorm) #Fit error of excitations exsN without ex0 (without re-fit)
            s         = self.getError(self.ftrc,ffit ,datnorm=self.ftrcnorm) #Fit error of excitations exsN without ex0 (with    re-fit)
            #s1        = self.getError(self.ftrc,ffit1,datnorm=self.ftrcnorm) #Fit error of excitations exsN with    ex0 alone refitted

        if allSignif:
            signifFit = (s-sN)/(s0-sN)
            #signifCon = np.sqrt((s0-sN)/(exsN.exlist[iex0].strength/sum([ex.strength for ex in exsN.exlist])))
            diperr    = np.linalg.norm(np.subtract(exs.exlist[iex0].dipoles,exsN.exlist[iex0].dipoles))/np.linalg.norm(exsN.exlist[iex0].dipoles)
            engerr    = abs(exs.exlist[iex0].energy-exsN.exlist[iex0].energy)/abs(exsN.exlist[iex0].energy)
            phaerr    = abs(exs.exlist[iex0].phase -exsN.exlist[iex0].phase )/max(abs(exsN.exlist[iex0].phase ),0.001) #0.001: Prevent "divide by zero"
            signifExc = 1.-np.linalg.norm([diperr,engerr,phaerr])/np.sqrt(3)
        else:
            signifFit =0.
            signifExc =0.
        dip       = exsN.exlist[iex0].dipole
        pol       = [self.dip[icalc][0].epol for icalc in range(self.ncalc)]
        signifAng = np.sqrt(np.max([abs(np.dot(dip/np.linalg.norm(dip),pol[icalc]/np.linalg.norm(pol[icalc]))) for icalc in range(self.ncalc)]))
        strerr    = abs(exsN.exlist[iex0].strengthErr/        exsN.exlist[iex0].strength        )
        #engerr    = abs(exsN.exlist[iex0].energyErr  /        exsN.exlist[iex0].energy          )
        engerr    = abs(exsN.exlist[iex0].energyErr  /        2.*piT                            ) #Use the expected energy range as reference (2pi/T)
        phaerr    = abs(exsN.exlist[iex0].phaseErr   /max(abs(exsN.exlist[iex0].phase   ),0.001))
        signifErr = 1.-np.linalg.norm([strerr,engerr,phaerr])/np.sqrt(3)
        w         = exsN.exlist[iex0].energy
        wc        = 0.5*(exsN.exlist[iex0].erange[1]+exsN.exlist[iex0].erange[0])
        dw        = 0.5*(exsN.exlist[iex0].erange[1]-exsN.exlist[iex0].erange[0])
        signifRng = 1.-(np.abs(w-wc)/np.abs(dw))**4
        Hw        = self.ext.getVal([exsN.exlist[iex0].energy])[0]
        phaana    = np.angle(Hw)%(2.*np.pi)
        phafit    = exsN.exlist[iex0].phase
        phadiff   = (phaana-phafit)%(2*np.pi) #in ]-2pi:2pi]
        if phadiff <= -np.pi: phadiff += 2*np.pi #Move into ]-pi:pi] range
        if phadiff >   np.pi: phadiff -= 2*np.pi
        signifPha = 1.-np.abs(phadiff)/np.pi

        exsN.exlist[iex0].setSignificance(signifFit,signifErr,signifAng,signifExc,signifRng,signifPha) #Compute and set significance

    #--------------------------------------------------------------------------#
    # Plot fit
    #--------------------------------------------------------------------------#
    def plotFitDebug(self,it=0):
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    for irc in range(2):
                        plt.title(f"it{it} calc{icalc} area{iarea} comp{icomp} rc{self.rc[irc]}")
                        plt.plot(self.freq,self.ft   [icalc][iarea][icomp][irc],label=f"Dat")
                        plt.plot(self.freq,self.ftfit[icalc][iarea][icomp][irc],label=f"Fit")
                        plt.show()

    #--------------------------------------------------------------------------#
    # Plot add-ine objective
    #--------------------------------------------------------------------------#
    def plotObjective(self,latex=False):
        scal        = np.ones(np.shape(self.scal))
        fdat        = self.ftrc
        ffit        = np.zeros(np.shape(self.ftfitrc))
        _, dat,_ ,_ = self.addExObj(fdat,ffit,scal)
        fdat        = np.zeros(np.shape(self.ftrc))
        ffit        = self.ftfitrc
        _, fit,_ ,_ = self.addExObj(fdat,ffit,scal)
        fdat        = self.ftrc   
        ffit        = self.ftfitrc
        _, err,_ ,_ = self.addExObj(fdat,ffit,scal)
        if latex:
            plt.rcParams["text.usetex"] = True
            plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
            plt.ylabel(r"$\|\boldsymbol{\Delta}_s\|^\prime$ [a.u.]")
        else:
            plt.ylabel(r"Add-line objective [a.u.]")
        plt.plot(self.freq,dat,label="Dat",color="#C0C0C3")
        plt.plot(self.freq,fit,label="Fit",color="#36454F")
        plt.plot(self.freq,err,label="Err",color="#DAA520")
# Set different ylimits
#        ax = plt.gca()
#        ylim = ax.get_ylim()
#        ax.set_ylim([ylim[0]*0.05, ylim[1]*0.05])
        plt.xlabel("Energy [Ry]")
        plt.legend(loc="upper right")
        plt.savefig("fitState.png",dpi=300)
        plt.show()
        head = 'Energy (Ry) | dat | fit | error '
        np.savetxt("fitState.dat",np.column_stack((self.freq,dat,fit,err)),header=head)
        if latex:
            plt.rcParams["text.latex.preamble"] = r""
            plt.rcParams["text.usetex"] = False

    #--------------------------------------------------------------------------#
    # Write fit and error functions
    #--------------------------------------------------------------------------#
    def writeFit(self,dbg=0):
        head = "Energy (Ry) | real(Fitf) | imag(Fitf)"
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                dip   = self.dip[icalc][iarea]
                dname = dip.dipname
                for icomp in range(self.ncomp):
                    fname = os.path.splitext(dname)[0]+'_fit_'+dip.descript[icomp]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq,self.ftfit[icalc][iarea][icomp][0],self.ftfit[icalc][iarea][icomp][1])),header=head)
        head = "Energy (Ry) | real(Errorf) | imag(Errorf)"
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                dip   = self.dip[icalc][iarea]
                dname = dip.dipname
                for icomp in range(self.ncomp):
                    fname = os.path.splitext(dname)[0]+'_err_'+dip.descript[icomp]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq,self.ft[icalc][iarea][icomp][0]-self.ftfit[icalc][iarea][icomp][0],self.ft[icalc][iarea][icomp][1]-self.ftfit[icalc][iarea][icomp][1])),header=head)
        if dbg>2:
            head = "Energy (Ry) | real(scalErrorf) | imag(scalErrorf)"
            for icalc in range(self.ncalc):
                for iarea in range(self.narea):
                    dip   = self.dip[icalc][iarea]
                    dname = dip.dipname
                    for icomp in range(self.ncomp):
                        fname = os.path.splitext(dname)[0]+'_scaerr_'+dip.descript[icomp]+'.dat'
                        np.savetxt(fname,np.column_stack((self.freq,(self.ft[icalc][iarea][icomp][0]-self.ftfit[icalc][iarea][icomp][0])/self.scal[icalc][iarea][icomp],(self.ft[icalc][iarea][icomp][1]-self.ftfit[icalc][iarea][icomp][1])/self.scal[icalc][iarea][icomp])),header=head)
            head = "Energy (Ry) | scal"
            for icalc in range(self.ncalc):
                for iarea in range(self.narea):
                    dip   = self.dip[icalc][iarea]
                    dname = dip.dipname
                    for icomp in range(self.ncomp):
                        fname = os.path.splitext(dname)[0]+'_sca_'+dip.descript[icomp]+'.dat'
                        np.savetxt(fname,np.column_stack((self.freq,self.scal[icalc][iarea][icomp])),header=head)
