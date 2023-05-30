import re
import numpy as np
from lmfit import Parameters
import matplotlib.pyplot as plt
#==============================================================================#
# Defines an excitation
#==============================================================================#
class Excitation:
    def __init__(self,narea=1,ncomp=3,exdict={},name="S",energy=0.01,erange=0.01,phase=0.01,tmod=1.,dipoles=None,fix=False):
        # Book keeping
        self.narea    = narea
        self.ncomp    = ncomp
        if not isinstance(dipoles,list):
            dipoles = [[0.01]*self.ncomp for i in range(self.narea)]
        # Fundamental
        self.name     = exdict.get("name"     ,   name     )
        self.fix      = exdict.get("fix"      ,   fix      ) #Global fix (by configuration file)
        self.fixtmp   = exdict.get("fix"      ,   fix      ) #Temporary fix (that can be released and is used for generating fit parameters)
        self.energy   = exdict.get("energy"   ,   energy   )
        if isinstance(erange,float):
            erange = [self.energy*(1.-erange),self.energy*(1.+erange)]
        self.erange   = exdict.get("erange"   ,   erange   ) #Energy range [min/max] (usually +/- 1% (i.e., x (1 +/- 0.01)) of the excitation energy)
        self.phase    = exdict.get("phase"    ,   phase    )
        self.tmod     = exdict.get("tmod"     ,   tmod     ) #T_eff = T*tmod (effective time for fitting)
        self.dipoles  = exdict.get("dipoles"  ,   dipoles  )
        # Derived
        self.dipole   = exdict.get("dipole"   ,  [0.]*ncomp)
        self.strength = exdict.get("strength" ,   0.       )
        self.strengths= exdict.get("strengths",  [0.]*narea) #Oscillator strength equivalent derived from the area's dipole (does NOT sum up to the total oscillator strength)
        # Update derived (overwrite the latter)
        self.derived()

        # Error checking
        if len(self.dipoles)!=self.narea or len(self.dipoles[0])!=self.ncomp: err(1,"Wrong shape of field 'dipoles'")
        if                                  len(self.dipole    )!=self.ncomp: err(1,"Wrong shape of field 'dipole'" )

    def todict(self):
        exdict = {}
        exdict["name"      ] = self.name
        exdict["fix"       ] = self.fix
        exdict["energy"    ] = float(self.energy) #Explicit casting to float is necesary for ruaml.dump
        exdict["erange"    ] = [float(self.erange[i]) for i in range(2)]
        exdict["phase"     ] = float(self.phase)
        exdict["tmod"      ] = float(self.tmod)
        exdict["dipoles"   ] = [[float(self.dipoles[i][j]) for j in range(len(self.dipoles[i]))] for i in range(len(self.dipoles))]
        exdict["dipole"    ] = [ float(self.dipole[i])     for i in range(len(self.dipole))]
        exdict["strength"  ] = float(self.strength)
        exdict["strengths" ] = [float(self.strengths[i])   for i in range(len(self.dipoles))]
        return exdict

    # Updates derived components
    def derived(self):
        #Define constants in Rydberg units
        m    = 0.5
        e2   = 2. #e^2
        hbar = 1.
        self.dipole   = np.array([sum([self.dipoles[iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)])
        self.strength = 2.*m*self.energy*np.linalg.norm(self.dipole)**2/(3.*e2*hbar)

        self.strengths= [2.*m*self.energy*np.linalg.norm(self.dipoles[iarea])**2/(3.*e2*hbar) for iarea in range(self.narea)]

    # Temporarily fixes the excitation
    def fixMe(self):
        self.fixtmp = True

    # Release temporary fix
    def releaseMe(self):
        self.fixtmp = self.fix

#============================================================================#
# Defines a list of excitations
#============================================================================#
class Excitations:
    def __init__(self,narea,ncomp,exlist=None):
        self.exlist = []
        self.narea  = narea
        self.ncomp  = ncomp
        try:
            for ex in exlist:
                self.add(exdict=ex,sort=False)
            self.sort()
            #Todo: Check if excitations in exlist match narea, ncomp
        except:
            pass

    #-------------------------------------------------------------------------
    # Add an excitation to the list
    #-------------------------------------------------------------------------
    def add(self,exdict=None,name="S",energy=0.,erange=0.01,phase=0.,tmod=1.,dipoles=None,fix=False,sort=True):
        #Todo: Check if new excitation matches narea, ncomp
        if isinstance(exdict,dict):
            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,exdict=exdict))
        else:
            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,name=name,energy=energy,erange=erange,phase=phase,tmod=tmod,dipoles=dipoles,fix=fix))
        if sort: self.sort()

    #-------------------------------------------------------------------------
    # Sorts the list of excitations and rename those with generic names "S[1234..]"
    #-------------------------------------------------------------------------
    def sort(self):
        self.exlist = sorted(self.exlist, key=lambda x: x.energy)
        pattern     = re.compile("^S[0-9]*$") #Matches and name "S<number>" or "S"
        for iex, ex in enumerate(self.exlist):
            if pattern.match(ex.name):
                ex.name = "S"+str(iex)

    #-------------------------------------------------------------------------
    # Set phases to zero (usually if boost excitation)
    #-------------------------------------------------------------------------
    def zeroPhase(self):
        for ex in self.exlist:
            ex.phase = 0.

    #-------------------------------------------------------------------------
    # Set tmod to one
    #-------------------------------------------------------------------------
    def unityTmod(self):
        for ex in self.exlist:
            ex.tmod = 1.

    #-------------------------------------------------------------------------
    # Generate lm fit parameters
    # Uses the fixtmp variable instead of fix
    #-------------------------------------------------------------------------
    def toParams(self,noEnergy=False,noPhase=False,noTmod=False):
        params = Parameters()
        for iex, ex in enumerate(self.exlist):
            params.add(f"w{iex}",vary=not ex.fixtmp and not noEnergy,value=ex.energy,min=ex.erange[0],max=ex.erange[1])
            params.add(f"p{iex}",vary=not ex.fixtmp and not noPhase ,value=ex.phase )
            params.add(f"t{iex}",vary=not ex.fixtmp and not noTmod  ,value=ex.tmod  )
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    params.add(f"d{iex}_{iarea}_{icomp}",vary=not ex.fixtmp,value=ex.dipoles[iarea][icomp])
        return params

    #-------------------------------------------------------------------------
    # Update from fit parameters
    #-------------------------------------------------------------------------
    def updateFromParam(self,params):
        pdict = params.valuesdict()
        for iex, ex in enumerate(self.exlist):
            ex.energy = pdict[f"w{iex}"]
            ex.phase  = pdict[f"p{iex}"]
            ex.tmod   = pdict[f"t{iex}"]
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    ex.dipoles[iarea][icomp] = pdict[f"d{iex}_{iarea}_{icomp}"]
            ex.derived()

    #-------------------------------------------------------------------------
    # Temporarily fix all or specific excitations
    #-------------------------------------------------------------------------
    def fix(self,which=None):
        for iex, ex in enumerate(self.exlist):
            if isinstance(which,list):
                if not iex in which: continue
            ex.fixMe()

    #-------------------------------------------------------------------------
    # Release all or specific excitations
    #-------------------------------------------------------------------------
    def release(self,which=None):
        for iex, ex in enumerate(self.exlist):
            if isinstance(which,list):
                if not iex in which: continue
            ex.releaseMe()

    #-------------------------------------------------------------------------
    # Count free excitations
    #-------------------------------------------------------------------------
    def countFree(self):
        nfree = 0
        for ex in self.exlist:
            if not ex.fixtmp: nfree +=1
        return nfree

    #-------------------------------------------------------------------------
    # Print excitations
    #-------------------------------------------------------------------------
    def print(self):
        for iex, ex in enumerate(self.exlist):
            ex.derived()
            print("Ex ",iex)
            print("  Name:     ", ex.name    )
            print("  Energy:   ", ex.energy  )
            print("  Erange:   ", ex.erange  )
            print("  Phase:    ", ex.phase   )
            print("  Tmod:     ", ex.tmod    )
            print("  Strength: ", ex.strength)
            print("  Dipole:   ", ex.dipole  )
            print("  Dipoles per area:")
            for iarea in range(self.narea):
                print("      ", iarea, ex.dipoles[iarea])

    #-------------------------------------------------------------------------
    # Return Lorentz function as well as energies, strengths, and labels
    #  - jex is a list of excitation indices. If None: Take all indices
    #  - jarea is an integer specifying the area. If None: Sum up all areas
    #  - jcomp is 0 (x), 1 (y), or 2 (z). if None: Use an oscillator-strength equivalent derived from the area's dipole
    #  - if jarea and jcomp are None: Use the oscillator strengths as height
    #-------------------------------------------------------------------------
    def lorentz(self,w,gamma,jarea=None,jcomp=None,jex=None):
        n  = len(w)
        f  = np.zeros(n,dtype=float)
        w0 = []
        f0 = []
        l0 = []
        for iex, ex in enumerate(self.exlist):
            if isinstance(jex,list):
                if not jex in iex: cycle
            ex.derived()
            w0.append(ex.energy)
            if isinstance(jarea,int) and isinstance(jcomp,int):
                f0.append(ex.dipoles[jarea][jcomp])
            elif isinstance(jarea,int):
                f0.append(ex.strengths[jarea])
            elif isinstance(jcomp,int):
                f0.append(ex.dipole[jcomp])
            else:
                f0.append(ex.strength)
            l0.append(ex.name)
            for i in range(n):
                gw  = gamma*w0[-1]
                gw2 = gw*gw
                dw  = w[i]-w0[-1]
                dw2 = dw*dw
                f[i] += f0[-1]*gw2/(dw2+gw2)
        return f, w0, f0, l0

    #-------------------------------------------------------------------------
    # Plot excitations
    #-------------------------------------------------------------------------
    def plot(self,wb,dw,gamma,jarea=None,jcomp=None,jex=None):
        w  = np.arange(wb[0],wb[1],dw)
        n  = len(w)
        f, w0, f0, l0 = self.lorentz(w,gamma,jarea,jcomp,jex)
        plt.title("Fitted Spectrum")
        plt.xlabel("Energy [Ry]")
        plt.ylabel("Oscillator Strength (Peak Height)")
        plt.plot(w,f)
        fmin = [0.]*len(f0)
        plt.vlines(w0,fmin,f0,colors="k")
        for iex in range(len(self.exlist)):
            plt.text(w0[iex],f0[iex]+0.01*np.amax(f0),l0[iex])
        plt.show()

    #-------------------------------------------------------------------------
    # Plot all excitations into panels
    #-------------------------------------------------------------------------
    def plotPanels(self,wb,dw,gamma,jex=None):
        w  = np.arange(wb[0],wb[1],dw)
        n  = len(w)
        fig, axs = plt.subplots(self.ncomp,self.narea,squeeze=False,sharex=True,sharey=True)
        fig.suptitle("Fitted Spectrum")
        for ax in axs.flat:
            ax.set(xlabel="Energy [Ry]", ylabel="Oscillator Strength (Peak Height)")
        for ax in axs.flat:# Hide x labels and tick labels for top plots and y ticks for right plots.
            ax.label_outer()
        for iarea in range(self.narea):
            for icomp in range(self.ncomp):
                f, w0, f0, l0 = self.lorentz(w,gamma,iarea,icomp,jex)
                axs[icomp][iarea].plot(w,f)
                fmin = [0.]*len(f0)
                axs[icomp][iarea].vlines(w0,fmin,f0,colors="k")
                for iex in range(len(self.exlist)):
                    axs[icomp][iarea].text(w0[iex],f0[iex]+0.01*np.amax(f0),l0[iex])
        plt.show()
