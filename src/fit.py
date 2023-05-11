#==============================================================================#
# Fitting Class
#==============================================================================#
#intrinsic
import os
import numpy as np
import math
from lmfit import Parameters, minimize, fit_report
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#own
import errorHandler as err
import dipole
import excitations

class Fit:
    def __init__(self,dip,ext,excit,fitrange):
        self.minw     = fitrange[0]
        self.maxw     = fitrange[1]
        self.dip      = dip
        self.ext      = ext
        self.exttype  = self.dip[0][0].ext
        self.ncalc    = len(self.dip)
        self.narea    = len(self.dip[0])
        self.ncomp    = len(self.dip[0][0].ft)
        self.excit    = excit
        self.freq     = np.array([dip[0][0].freq    [i] for i in range(len(dip[0][0].freq    )) if dip[0][0].freq    [i]>=self.minw and dip[0][0].freq    [i]<=self.maxw])
        self.freqPade = np.array([dip[0][0].freqPade[i] for i in range(len(dip[0][0].freqPade)) if dip[0][0].freqPade[i]>=self.minw and dip[0][0].freqPade[i]<=self.maxw])
        self.Nf       = len(self.freq)
        self.NfPade   = len(self.freqPade)
        self.ft       = []
        self.pw       = []
        self.pade     = []
        for icalc in range(self.ncalc):
            self.ft  .append([])
            self.pw  .append([])
            self.pade.append([])
            for iarea in range(self.narea):
                d = self.dip[icalc][iarea]
                self.ft  [icalc].append([])
                self.pw  [icalc].append([])
                self.pade[icalc].append([])
                for n in range(len(d.ft)):
                    tmpr         =                     np.array([d.ft  [n][0][i] for i in range(len(d.ft  [n][0])) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw])
                    tmpi         =                     np.array([d.ft  [n][1][i] for i in range(len(d.ft  [n][1])) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw])
                    self.ft      [icalc][iarea].append([tmpr,tmpi])
                    self.pw      [icalc][iarea].append(np.array([d.pw  [n]   [i] for i in range(len(d.pw  [n]   )) if d.freq    [i]>=self.minw and d.freq    [i]<=self.maxw]))
                    tmpr         =                     np.array([d.pade[n][0][i] for i in range(len(d.pade[n][0])) if d.freqPade[i]>=self.minw and d.freqPade[i]<=self.maxw])
                    tmpi         =                     np.array([d.pade[n][1][i] for i in range(len(d.pade[n][1])) if d.freqPade[i]>=self.minw and d.freqPade[i]<=self.maxw])
                    self.pade    [icalc][iarea].append([tmpr,tmpi])
        self.fitf   = self.getFitFunc(self.excit)
        self.fiterr = self.getError(self.ft,self.fitf)
        self.minerr = np.amin   (self.fiterr)
        self.maxerr = np.amax   (self.fiterr)
        self.avgerr = np.average(self.fiterr)
        self.breakMinimization = False
        self.runningError = 1000000.

    #--------------------------------------------------------------------------#
    # Get the spectral line function
    #--------------------------------------------------------------------------#
    def lineFunc(self,T,w):
        sinc  = np.sinc(w*T/np.pi) #np.sinc is defined as sin(pi*x)/(pi*x)
        cosc  = (1.-math.cos(w*T))/(w*T) if abs(w)>0. else 0.
        return T*(cosc + 1.j*sinc)

    #--------------------------------------------------------------------------#
    # Get the fit function for a given set of frequencies and excitations
    #--------------------------------------------------------------------------#
    def getFitFunc(self,excit):
        f  = []
        for icalc in range(self.ncalc):
            f.append([])
            for iarea in range(self.narea):
                f[icalc].append([])
                for n in range(self.ncomp):
                    fr = np.zeros(self.Nf,dtype=float)
                    fi = np.zeros(self.Nf,dtype=float)
                    f[icalc][iarea].append([fr,fi])
        T  = self.dip[0][0].tprop
        Ef = self.dip[0][0].efield
        for ex in excit.exlist:
            Hw = self.ext.getVal([ex.energy])
            for icalc in range(self.ncalc):
                Ep = self.dip[icalc][0].epol #Polarization can differ between calculations
                for iarea in range(self.narea):
                    mua= ex.dipoles[iarea]
                    mu = ex.dipole
                    for n in range(self.ncomp):
                        ampl   = Ef*np.abs(np.dot(Ep,mu))*np.abs(Hw[0])*mua[n]
                        for i in range(self.Nf):
                            wm     = self.freq[i]-ex.energy
                            wp     = self.freq[i]+ex.energy
                            tmp    = ampl*(\
                                np.exp(-1.0j*ex.phase)*self.lineFunc(T,wm) -\
                                np.exp(+1.0j*ex.phase)*self.lineFunc(T,wp) )
                            f[icalc][iarea][n][0][i] += np.real(tmp)
                            f[icalc][iarea][n][1][i] += np.imag(tmp)
        return f


    #--------------------------------------------------------------------------#
    # Errors dat vs. fit functions
    #--------------------------------------------------------------------------#
    def getError(self,dat,fit,relative=True):
        e = []
        a = 0.
        for icalc in range(self.ncalc):
            e.append([])
            for iarea in range(self.narea):
                e[icalc].append([])
                for icomp in range(self.ncomp):
                    er = np.linalg.norm([dat[icalc][iarea][icomp][0][i]-fit[icalc][iarea][icomp][0][i] for i in range(len(dat[icalc][iarea][icomp][0]))])
                    ei = np.linalg.norm([dat[icalc][iarea][icomp][1][i]-fit[icalc][iarea][icomp][1][i] for i in range(len(dat[icalc][iarea][icomp][1]))])
                    a += np.linalg.norm([dat[icalc][iarea][icomp][0][i]                                for i in range(len(dat[icalc][iarea][icomp][0]))])
                    a += np.linalg.norm([dat[icalc][iarea][icomp][1][i]                                for i in range(len(dat[icalc][iarea][icomp][1]))])
                    e[icalc][iarea].append([er   ,ei   ])
        if relative:
            for icalc in range(self.ncalc):
                for iarea in range(self.narea):
                    for icomp in range(self.ncomp):
                        for rc in range(2):
                            e[icalc][iarea][icomp][rc] /= a
        return e

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
                    for rc in range(2):
                        for i in range(self.NfPade):
                            f[i] += fdat[icalc][iarea][n][rc][i]**2
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
        phase = (np.angle(Hw) - energy*t0)%(2.*math.pi)
        # Get the closest index in the FT array
        dw  = (self.freq[-1]-self.freq[0])/(len(self.freq)-1)
        idx = int(np.rint((energy-self.freq[0])/dw))
        # Find the calculation for which the excitation's peak is largest
        ampl        = [[[self.ft[icalc][iarea][icomp][0][idx]+1.0j*self.ft[icalc][iarea][icomp][1][idx] for icomp in range(self.ncomp)] for iarea in range(self.narea)] for icalc in range(self.ncalc)]
        jcalc       = np.argmax([np.linalg.norm([sum([ampl[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)]) for icalc in range(self.ncalc)]) #Use this calculation for the dipole-moment guess
        # Approx. dipole directions
        T       = self.dip[    0][0].tprop
        Ef      = self.dip[    0][0].efield
        Ep      = self.dip[jcalc][0].epol
        heightArea  = [[np.imag(np.exp(-1.0j*phase)*     ampl[jcalc][iarea][icomp]) for iarea in range(self.narea)]   for icomp in range(self.ncomp)] #Rotate the height by e^-i*phi and take its imag part to get the line heights (with sign)
        height      = [ np.imag(np.exp(-1.0j*phase)*sum([ampl[jcalc][iarea][icomp]  for iarea in range(self.narea)])) for icomp in range(self.ncomp)]
        dipdir      = height/np.linalg.norm(height) # Direction of the dipole moment
        eped        = np.dot(dipdir,Ep) # Scalar product e_pol.e_mu
        # Approx. dipoles
        dipabs  = math.sqrt(np.linalg.norm(height)/(T*Ef*np.abs(eped)*np.abs(Hw))) #Absolute value of the total dipole moment
        dipoles = []
        for iarea in range(self.narea):
            dipoles.append([])
            for icomp in range(self.ncomp):
                dipoles[iarea].append(heightArea[icomp][iarea]/(T*Ef*eped*np.abs(Hw)*dipabs)) # hightarea[icomp] = T*Ef*eped*abs(Hw)*abs(mu)*mu[icomp] -> solve for mu[icomp]

        return phase, dipoles
            
    #--------------------------------------------------------------------------#
    # Search for the highest peak in the given function
    #--------------------------------------------------------------------------#
    def findPeaks(self,freq,f,dbg=0,hf=0.05,mw=100):
        dw = (freq[-1]-freq[0])/(len(freq)-1)
        pos, prop = find_peaks(f,height=hf*np.amax(f),width=(0.,100))
      
        if dbg>1:
            plt.plot(freq,f)
            plt.plot(freq[pos], f[pos], "x")
            plt.show()

        return pos

    #--------------------------------------------------------------------------#
    # Fit objective
    #--------------------------------------------------------------------------#
    def fitObj(self,params,dbg=0):
        self.excit.updateFromParam(params)
        ffit = np.array(self.getFitFunc(self.excit)).flatten()
        fdat = np.array(self.ft).flatten()
        if self.breakMinimization==True:
            obj  = np.zeros(len(ffit),dtype=float)
        else:
            obj  = np.subtract(fdat,ffit)
        return obj

    #--------------------------------------------------------------------------#
    # Fit inter - is called in every minimizer iteration
    #--------------------------------------------------------------------------#
    def fitInter(self,params,iter,resid,dbg=0,breakmod=5):
        currentError = np.sum(np.abs(resid))
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
                print("WARNING: Minimization stuck; abort and fall back to best parameter set")
                self.breakMinimization = True
            if iter>-1: self.runningError = currentError
        return

    #--------------------------------------------------------------------------#
    # Fit Atomic
    #--------------------------------------------------------------------------#
    def fitAtomic(self,dbg=0):
        params = self.excit.toParams()
        fitres = minimize(self.fitObj,params,iter_cb=self.fitInter,kws={"dbg":dbg})
        if self.breakMinimization: #If minimization stuck and was force-aborted, reset the parameters to the best ones reached
            params = self.bestParams
            self.breakMinimization = False
        if dbg>0: print(fit_report(fitres))
        self.excit.updateFromParam(fitres.params)

    #--------------------------------------------------------------------------#
    # Add new excitation
    #--------------------------------------------------------------------------#
    def addEx(self):
        maxval = 0.
        maxen  = 0.
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    for rc in range(2):
                        fdat   = self.ft  [icalc][iarea][icomp][rc]
                        ffit   = self.fitf[icalc][iarea][icomp][rc]
                        r      = np.abs(np.subtract(fdat,ffit))
                        maxpos = np.argmax(r)
                        if r[maxpos]>maxval:
                            maxval = r[maxpos]
                            maxen  = self.freq[maxpos]
        phase, dipoles = self.guessExcit(maxen)
        self.excit.add(energy=maxen,phase=phase,dipoles=dipoles)
        self.excit.print()

    #--------------------------------------------------------------------------#
    # Fit Wrapper
    #--------------------------------------------------------------------------#
    def fit(self,dbg=0,maxit=10,tol=0.05,skipfirst=False):
        self.update(0,dbg=dbg)
        if not skipfirst:
            self.fitAtomic(dbg=dbg)   # Fit existing excitations
        for it in range(1,maxit+1):
            self.excit.fix()          # Temporarily fix all existing excitations
            self.addEx()              # Add new excitation
            self.fitAtomic(dbg=dbg)   # Fit new excitation alone
            self.update(it+1,dbg=dbg) # Update some self.components
            self.excit.release()      # Release all temporarily fixed excitations
            self.fitAtomic(dbg=dbg)   # Fit all non-permanently fixed excitations
            self.update(it+1,dbg=dbg) # Update some self.components
            if self.avgerr<tol: break
        return self.excit

    #--------------------------------------------------------------------------#
    # Update fitf, fiterr, minerr, maxerr, afgerr based on self.excit and report
    #--------------------------------------------------------------------------#
    def update(self,it,dbg=0):
        self.fitf   = self.getFitFunc(self.excit)
        self.fiterr = self.getError(self.ft,self.fitf)
        self.minerr = np.amin   (self.fiterr)
        self.maxerr = np.amax   (self.fiterr)
        self.avgerr = np.average(self.fiterr)
        if dbg>0: self.reportFit(it)
        if dbg>1: self.plotFitDebug(it)

    #--------------------------------------------------------------------------#
    # Report fit
    #--------------------------------------------------------------------------#
    def reportFit(self,it):
        print("Fit report: ", it)
        print("")
        print("Excitations")
        self.excit.print()
        print("")
        print("Error: Min | Max | Avg")
        print("  ",self.minerr,self.maxerr,self.avgerr)
        print("")
        print("Errors:")
        for icalc in range(self.ncalc):
            print("  ",icalc)
            for iarea in range(self.narea):
                print("    ",iarea)
                for icomp in range(self.ncomp):
                    print("",self.fiterr[icalc][iarea][icomp])

    #--------------------------------------------------------------------------#
    # Plot fit
    #--------------------------------------------------------------------------#
    def plotFitDebug(self,it=0):
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    for rc in range(2):
                        plt.title(f"it{it} calc{icalc} area{iarea} comp{icomp} rc{rc}")
                        plt.plot(self.freq,self.ft  [icalc][iarea][icomp][rc],label=f"Dat")
                        plt.plot(self.freq,self.fitf[icalc][iarea][icomp][rc],label=f"Fit")
                        plt.show()

    #--------------------------------------------------------------------------#
    # Write fit functions
    #--------------------------------------------------------------------------#
    def writeFit(self):
        head = 'Energy (Ry) | real(Fitf) | imag(Fitf)'
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                dip   = self.dip[icalc][iarea]
                dname = dip.dipname
                for icomp in range(self.ncomp):
                    fname = os.path.splitext(dname)[0]+'_fit_'+dip.descript[icomp]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq,self.fitf[icalc][iarea][icomp][0],self.fitf[icalc][iarea][icomp][1])),header=head)
