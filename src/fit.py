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
from mathtools import fspectrum #, butter_lowpass_filter

class Fit:
    def __init__(self,dip,ext,excit,fitrange):
        self.minw     = fitrange[0]
        self.maxw     = fitrange[1]
        self.dip      = dip
        self.ext      = ext
        self.exttype  = self.dip[0][0].ext
        self.tprop    = self.dip[0][0].tprop
        self.fwhm     = np.pi/self.tprop
        self.ncalc    = len(self.dip)
        self.narea    = len(self.dip[0])
        self.ncomp    = len(self.dip[0][0].ft)
        self.excit    = excit
        self.freq     = np.array([dip[0][0].freq    [i] for i in range(len(dip[0][0].freq    )) if dip[0][0].freq    [i]>=self.minw and dip[0][0].freq    [i]<=self.maxw])
        self.dw       = (self.freq[-1]-self.freq[0])/(len(self.freq)-1)
        self.freqPade = np.array([dip[0][0].freqPade[i] for i in range(len(dip[0][0].freqPade)) if dip[0][0].freqPade[i]>=self.minw and dip[0][0].freqPade[i]<=self.maxw])
        self.Nf       = len(self.freq)
        self.NfPade   = len(self.freqPade)
        self.rc       = np.array([1]) if self.dip[0][0].ext=="boost" else np.array([0,1]) #Use only imag(rc==1) part if boost excitation or real(rc==0) and imag(rc==1) part else
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
                    self.pw  [icalc][iarea].append(np.array([d.pw  [n][i] for i in range(len(d.pw  [n]   )) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw]))
                    self.ftrc[icalc][iarea].append([])
                    for irc in range(2):
                        self.ft  [icalc][iarea][n].append(np.array([d.ft  [n][        irc ][i] for i in range(len(d.ft  [n][        irc ])) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw]))
                        self.pade[icalc][iarea][n].append(np.array([d.pade[n][        irc ][i] for i in range(len(d.pade[n][        irc ])) if d.freqPade[i]>=self.minw and d.freqPade[i]<=self.maxw]))
                    for irc in range(self.nrc):
                        self.ftrc[icalc][iarea][n].append(np.array([d.ft  [n][self.rc[irc]][i] for i in range(len(d.ft  [n][self.rc[irc]])) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw]))
        self.ft       = np.array(self.ft  )
        self.pw       = np.array(self.pw  )
        self.ftrc     = np.array(self.ftrc) #Only imag part if boost; used for fit
        self.ftfit    = self.getFitFunc(self.excit)
        self.pwfit    = self.getPw(self.ftfit)
        self.ftfitrc  = self.getFitFunc(self.excit,self.rc) #only imag part if boost
        self.ftrcnorm = np.linalg.norm(self.ftrc.flatten()) #Use ftrc as reference
        self.fiterr   = self.getError(self.ftrc,self.ftfitrc,datnorm=self.ftrcnorm)
        self.breakMinimization = False
        self.runningError = 1000000.

    #--------------------------------------------------------------------------#
    # Get the fit function for a given set of frequencies and excitations
    #--------------------------------------------------------------------------#
    def getFitFunc(self,excit,rc=[0,1]):
        # Setup new array
        nrc = len(rc)
        f  = np.zeros((self.ncalc,self.narea,self.ncomp,nrc,self.Nf),dtype=float)
        # Get relevant excitation measures
        T  = self.tprop
        Ef = self.dip[0][0].efield
        nex = len(excit.exlist)
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
                        ampl[icalc][iarea][n][iex] = Ef*np.abs(np.dot(Ep,ex.dipole))*np.abs(Hw[iex])*ex.dipoles[iarea][n]
        return fspectrum(self.ncalc,self.narea,self.ncomp,np.array(rc),self.Nf,T,self.freq,energy,phase,tmod,ampl)

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
        fdat = np.array(dat).flatten()
        ffit = np.array(fit).flatten()
        if isinstance(datnorm,float):
            return np.linalg.norm(np.subtract(fdat,ffit))/datnorm
        else:
            return np.linalg.norm(np.subtract(fdat,ffit))/np.linalg.norm(fdat)

    #--------------------------------------------------------------------------#
    # Make an initial guess based on the sum of the Pade power spectra and the FT spectra
    #--------------------------------------------------------------------------#
    def newGuess(self,hf=0.05):
        # Init empty excitation object and arrays
        self.excit = excitations.Excitations(narea=self.narea,ncomp=self.ncomp)
        fdat  = self.pade
        f     = np.zeros(len(self.freqPade),dtype=float)
        # Sum up squared Pade spectra
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for n in range(self.ncomp):
                    for irc in range(2):
                        for i in range(self.NfPade):
                            f[i] += fdat[icalc][iarea][n][irc][i]**2
        # Search most prominent peaks
        peaks = self.findPeaks(self.freqPade,np.sqrt(f),dbg=2,hf=hf)
        # For each peak
        for ipeak in range(len(peaks)):
            #Get energy and guess phase and dipoles based on the FT
            energy  = self.freqPade[peaks[ipeak]]
            phase, dipoles = self.guessExcit(energy)
            #Add excitation to list
            self.excit.add(energy=energy,phase=phase,dipoles=dipoles)

    #--------------------------------------------------------------------------#
    # Guess phase and dipole moments of an excitation at an approx. energy
    #--------------------------------------------------------------------------#
    def guessExcit(self,energy):
        # Get the excitation strength/phase at that energy
        Hw    = self.ext.getVal([energy])[0]
        # Set some abbreviations (only the polarization may differ between different calculations)
        t0    = self.dip[0][0].t0
        # Approx. phase: phi = phi_ext - t0*energy (the latter because dipole moment time frame is shifted by t0=excitation period)
        #phase = 0. if self.exttype=="boost" else (np.angle(Hw) - energy*t0)%(2.*np.pi)
        phase = (np.angle(Hw) - energy*t0)%(2.*np.pi)
        # Get the closest index in the FT array
        idx = int(np.rint((energy-self.freq[0])/self.dw))
        # Find the calculation for which the excitation's peak is largest
        ampl        = [[[self.ft[icalc][iarea][icomp][0][idx]+1.0j*self.ft[icalc][iarea][icomp][1][idx] for icomp in range(self.ncomp)] for iarea in range(self.narea)] for icalc in range(self.ncalc)]
        jcalc       = np.argmax([np.linalg.norm([sum([ampl[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)]) for icalc in range(self.ncalc)]) #Use this calculation for the dipole-moment guess
        # Approx. dipole directions
        T           = self.tprop
        Ef          = self.dip[    0][0].efield
        Ep          = self.dip[jcalc][0].epol
        heightArea  = [[np.imag(np.exp(-1.0j*phase)*     ampl[jcalc][iarea][icomp]) for iarea in range(self.narea)]   for icomp in range(self.ncomp)] #Rotate the height by e^-i*phi and take its imag part to get the line heights (with sign)
        height      = [ np.imag(np.exp(-1.0j*phase)*sum([ampl[jcalc][iarea][icomp]  for iarea in range(self.narea)])) for icomp in range(self.ncomp)]
        dipdir      = height/np.linalg.norm(height) # Direction of the dipole moment
        eped        = np.dot(dipdir,Ep) # Scalar product e_pol.e_mu
        # Approx. dipoles
        dipabs  = np.sqrt(np.linalg.norm(height)/(T*Ef*np.abs(eped)*np.abs(Hw))) #Absolute value of the total dipole moment
        dipoles = []
        for iarea in range(self.narea):
            dipoles.append([])
            for icomp in range(self.ncomp):
                dipoles[iarea].append(heightArea[icomp][iarea]/(T*Ef*eped*np.abs(Hw)*dipabs)) # hightarea[icomp] = T*Ef*eped*abs(Hw)*abs(mu)*mu[icomp] -> solve for mu[icomp]

        return phase, dipoles
            
    #--------------------------------------------------------------------------#
    # Search for the highest peak in the given function
    #--------------------------------------------------------------------------#
    def findPeaks(self,freq,f,dbg=0,hf=0.05):
        pos, prop = find_peaks(f,height=hf*np.amax(f))#,width=(0.,100))
      
        if dbg>1:
            plt.plot(freq,f)
            plt.plot(freq[pos], f[pos], "x")
            plt.show()

        return pos

    #--------------------------------------------------------------------------#
    # Fit objective
    # nfree is not used here but required since fitInter needs it
    #--------------------------------------------------------------------------#
    def fitObj(self,params,dbg=0,nfree=0):
        self.excit.updateFromParam(params)
        try:
            ffit = self.getFitFunc(self.excit,self.rc).flatten()
        except:
            raise
        fdat = self.ftrc.flatten()
        if self.breakMinimization==True:
            obj  = np.zeros(len(ffit),dtype=float)
        else:
            obj  = np.subtract(fdat,ffit)
        return obj

    #--------------------------------------------------------------------------#
    # Fit inter - is called in every minimizer iteration
    #--------------------------------------------------------------------------#
    def fitInter(self,params,iter,resid,dbg=0,breakmod=5,nfree=0):
        if self.breakMinimization: return
        currentError = np.sum(np.abs(resid))/self.ftrcnorm
        if iter==-1:
            self.bestParams   = params
            self.bestError    = currentError
            self.runningError = currentError
        else:
            if currentError<self.bestError:
                self.bestError  = currentError
                self.bestParams = params
        if dbg>0: print("DEBUG: Objective Norm: ", iter, currentError)
        if iter%(breakmod*len(params.valuesdict()))==0 and iter>0:
            if abs(currentError/self.runningError)>0.99:
                print("WARNING: Minimization stuck; abort and fall back to best parameter set",file=sys.stderr)
                self.breakMinimization = True
            if iter>-1: self.runningError = currentError

    #--------------------------------------------------------------------------#
    # Fit Atomic
    #--------------------------------------------------------------------------#
    def fitAtomic(self,dbg=0):
        params = self.excit.toParams(noPhase=self.exttype=="boost",noTmod=True) #fix phase if boost and fix time always
        #params = self.excit.toParams() #free phases and time modifier
        nfree  = self.excit.countFree()
        try:
            fitres = minimize(self.fitObj,params,iter_cb=self.fitInter,kws={"dbg":dbg,"nfree":nfree})
        except:
            raise #self.excit is not updated in this case
        if self.breakMinimization: #If minimization stuck and was force-aborted, reset the parameters to the best ones reached
            params = self.bestParams
            self.breakMinimization = False
        if dbg>0: print(fit_report(fitres))
        self.excit.updateFromParam(fitres.params)

    #--------------------------------------------------------------------------#
    # Add new excitation
    #--------------------------------------------------------------------------#
    def addEx(self,dbg=0):
        fdat   = self.pw   .flatten()
        ffit   = self.pwfit.flatten()
        #diff   = butter_lowpass_filter(np.subtract(fdat,ffit),0.1/self.fwhm,1./self.dw,order=5)
        diff   = np.subtract(fdat,ffit)
        maxen  = self.freq[np.argmax(np.abs(diff))%self.Nf]
        phase, dipoles = self.guessExcit(maxen)
        #######
        #print(maxen, phase, dipoles)
        #plt.plot(np.array([self.freq]*3).flatten(),diff)
        #plt.show()
        #######
        self.excit.add(energy=maxen,phase=phase,dipoles=dipoles)
        if dbg>0: self.excit.print()

    #--------------------------------------------------------------------------#
    # Fit Wrapper
    #--------------------------------------------------------------------------#
    def fit(self,dbg=0,maxex=10,tol=0.05,skipfirst=False):
        if (self.exttype=="boost"): self.excit.zeroPhase()
        self.update(dbg=dbg)
        if not skipfirst:
            self.fitAtomic(dbg=dbg)   # Fit existing excitations
            self.update(dbg=dbg)      # Update some self.components
            self.reportFit(dbg=dbg)
        nex = len(self.excit.exlist)
        for it in range(nex+1,maxex+1):
            self.excit.fix()          # Temporarily fix all existing excitations
            self.addEx(dbg=dbg)       # Add new excitation
            try:
                self.fitAtomic(dbg=dbg)   # Fit new excitation alone
                self.update(dbg=dbg)      # Update some self.components
                self.excit.release()      # Release all temporarily fixed excitations
            except:
                break #self.excit is not updated in this case
            try:
                self.fitAtomic(dbg=dbg)   # Fit all non-permanently fixed excitations
                self.update(dbg=dbg)      # Update some self.components
            except:
                break #self.excit is not updated in this case
            self.reportFit(dbg=dbg)
            if dbg>1: self.plotFitDebug(it)
            if self.fiterr<tol: break
        return self.excit

    #--------------------------------------------------------------------------#
    # Update ftfit and fiterr based on self.excit and report
    #--------------------------------------------------------------------------#
    def update(self,dbg=0):
        self.ftfit   = self.getFitFunc(self.excit)
        self.pwfit   = self.getPw     (self.ftfit)
        self.ftfitrc = self.getFitFunc(self.excit,self.rc)
        self.fiterr  = self.getError(self.ftrc,self.ftfitrc,datnorm=self.ftrcnorm)

    #--------------------------------------------------------------------------#
    # Report fit
    #--------------------------------------------------------------------------#
    def reportFit(self,dbg=0):
        print("Nex | Error: ",len(self.excit.exlist),self.fiterr)
        if dbg>0:
            print("")
            print("Excitations")
            self.excit.print()
            print("")

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
    # Write fit functions
    #--------------------------------------------------------------------------#
    def writeFit(self):
        head = "Energy (Ry) | real(Fitf) | imag(Fitf)"
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                dip   = self.dip[icalc][iarea]
                dname = dip.dipname
                for icomp in range(self.ncomp):
                    fname = os.path.splitext(dname)[0]+'_fit_'+dip.descript[icomp]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq,self.ftfit[icalc][iarea][icomp][0],self.ftfit[icalc][iarea][icomp][1])),header=head)
