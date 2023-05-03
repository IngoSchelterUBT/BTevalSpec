#==============================================================================#
# Fitting Class
#==============================================================================#
#intrinsic
import numpy as np
import math
from lmfit import Parameters, minimize, report_fit
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
        self.averr  = np.average(self.fiterr)


    #--------------------------------------------------------------------------#
    # Get the fit function for a given set of frequencies and excitations
    #--------------------------------------------------------------------------#
    def getFitFunc(self,excit):
        f  = []
        T  = self.dip[0][0].tprop
        Ef = self.dip[0][0].efield
        Ep = self.dip[0][0].epol
        for icalc in range(excit.ncalc):
            f.append([])
            for iarea in range(excit.narea):
                f[icalc].append([])
                for n in range(excit.ncomp):
                    fr = np.zeros(self.Nf,dtype=float)
                    fi = np.zeros(self.Nf,dtype=float)
                    f[icalc][iarea].append([fr,fi])
        for ex in excit.exlist:
            Hw = self.ext.getVal([excit.exlist.energy])
            for icalc in range(excit.ncalc):
                for iarea in range(excit.narea):
                    for n in range(excit.ncomp):
                        for i in range(self.Nf):
                            wm     = self.freq[i]-ex.energy
                            wp     = self.freq[i]+ex.energy
                            mu     = ex.dipoles[iarea]

                            ampl   = -Ef*np.dot(Ep,mu)*Hw*mu[n] #a_j = -k.mu_j/e*1 mu_j = -E.mu_j/hbar H(w) mu_j -> Ry: -|E|e_E.mu_j H(w) mu_j
                            tmp    = ampl*(\
                                np.exp(+1.0j*ex.phase)*(math.cos(wm*T)-1.+1.j*math.sin(wm*T))/wm -\
                                np.exp(-1.0j*ex.phase)*(math.cos(wp*T)-1.+1.j*math.sin(wp*T))/wp )
                            f[icalc][iarea][n][0][i] += np.real(tmp)
                            f[icalc][iarea][n][1][i] += np.imag(tmp)
        return f

    #--------------------------------------------------------------------------#
    # Errors dat vs. fit functions
    #--------------------------------------------------------------------------#
    def getError(self,dat,fit):
        e  = []
        for icalc in range(len(dat)):
            e.append([])
            for iarea in range(len(dat[icalc])):
                e[icalc].append([])
                for n in range(len(dat[icalc][iarea])):
                    er = np.linalg.norm([dat[icalc][iarea][0][i]-fit[icalc][iarea][0][i] for i in range(len(dat[icalc][iarea][0]))])
                    ei = np.linalg.norm([dat[icalc][iarea][1][i]-fit[icalc][iarea][1][i] for i in range(len(dat[icalc][iarea][1]))])
                    e[icalc][iarea].append([er,ei])
        return e

    #--------------------------------------------------------------------------#
    # Make an initial guess
    #--------------------------------------------------------------------------#
    def newGuess(self,hf=0.05):
        excit = excitations.Excitations(ncalc=self.ncalc,narea=self.narea,ncomp=self.ncomp)
        fdat  = self.pade
        f     = np.zeros(len(self.freq),dtype=float)
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                for n in range(self.ncomp):
                    for rc in range(2):
                        for i in range(self.Nf):
                            f[i] += fdat[icalc][iarea][n][rc][i]**2
        peaks = self.findPeaks(np.sqrt(f),dbg=2,hf=hf)
        if self.exttype=="boost":
            phi = 0.
        else:
            phi = 0.1 #Todo: Make a better guess based on the external field
        for ipeak in range(len(peaks)):
            excit.add(energy=self.freqPade[peaks[ipeak]],phase=phi)

        return excit

    #--------------------------------------------------------------------------#
    # Search for the highest peak in the given function
    #--------------------------------------------------------------------------#
    def findPeaks(self,f,dbg=0,hf=0.05,mw=100):
        dw = (self.freqPade[-1]-self.freqPade[0])/(len(self.freqPade)-1)
        pos, prop = find_peaks(f,height=0.05*np.amax(f),width=(0.,100))
      
        if dbg>1:
            plt.plot(f)
            plt.plot(pos, f[pos], "x")
            plt.show()

        return pos
