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
            erange = [self.energy-erange,self.energy+erange]
        self.erange   = exdict.get("erange"   ,   erange   ) #Energy range [min/max] (usually +/- 1% (i.e., x (1 +/- 0.01)) of the excitation energy)
        self.phase    = exdict.get("phase"    ,   phase    )
        self.tmod     = exdict.get("tmod"     ,   tmod     ) #T_eff = T*tmod (effective time for fitting)
        self.dipoles  = exdict.get("dipoles"  ,   dipoles  )
        # Derived
        self.dipole   = exdict.get("dipole"   ,  [0.]*ncomp)
        self.strength = exdict.get("strength" ,   0.       )
        self.strengths= exdict.get("strengths",  [0.]*narea) #Oscillator strength equivalent derived from the area's dipole (does NOT sum up to the total oscillator strength)
        #Significance
        self.signif   = exdict.get("signif"   ,   0.       ) #Significance of the excitation
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
        exdict["signif"    ] = float(self.signif)
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

    # Set significance
    def setSignificance(self,signif):
        self.signif = signif

#============================================================================#
# Defines a list of excitations
#============================================================================#
class Excitations:
    def __init__(self,narea,ncomp,exlist=None):
        self.exlist = []
        self.narea  = narea
        self.ncomp  = ncomp
        if isinstance(exlist,list):
            for ex in exlist:
                if isinstance(ex,Excitation):
                    self.add(exdict=ex.todict(),sort=False)
                elif isinstance(ex,dict):
                    self.add(exdict=ex,sort=False)
            self.sort()
            #Todo: Check if excitations in exlist match narea, ncomp

    #-------------------------------------------------------------------------
    # Add an excitation to the list
    # Returns the index of the new excitation inside self.exlist
    #-------------------------------------------------------------------------
    def add(self,exdict={},name="S",energy=0.,erange=0.01,phase=0.,tmod=1.,dipoles=None,fix=False,sort=True):
        #Todo: Check if new excitation matches narea, ncomp
#        if isinstance(exdict,dict):
#            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,exdict=exdict))
#        else:
#            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,name=name,energy=energy,erange=erange,phase=phase,tmod=tmod,dipoles=dipoles,fix=fix))
        ex = Excitation(narea=self.narea,ncomp=self.ncomp,exdict=exdict,name=name,energy=energy,erange=erange,phase=phase,tmod=tmod,dipoles=dipoles,fix=fix)
        self.exlist.append(ex)
        if sort: self.sort()
        #return self.exlist.index(ex)

    #-------------------------------------------------------------------------
    # Copy excitations object (potentially remove excitations)
    #-------------------------------------------------------------------------
    def copy(self):
        return Excitations(self.narea,self.ncomp,exlist=self.exlist)

    #-------------------------------------------------------------------------
    # Remove excitation by index or name
    #-------------------------------------------------------------------------
    def remove(self,rmidx=[],rmname=[]):
        #Remove by indices first and extend rmname list
        for iex, ex in enumerate(self.exlist):
            if iex in rmidx: rmname.append(ex.name)
            #Don't remove yet since this would change the original indexing
        #Then remove by names
        for iex, ex in enumerate(self.exlist):
            #if ex.name in rmname: del self.exlist[iex]
            if ex.name in rmname: self.exlist.remove(ex)

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
    # Temporarily fix excitations in a given interval (or the inverse interval)
    #-------------------------------------------------------------------------
    def fixErange(self,erange,inverse=False):
        for iex, ex in enumerate(self.exlist):
            if ex.energy < erange[1] and ex.energy >erange[0]: # inside inverval
                if not inverse: ex.fixMe()                     # fix if not inverse
            else:                                              #outside inverval
                if     inverse: ex.fixMe()                     # fix if     inverse

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
    # Restrict range of energy parameter
    # If neigh=True: Additionally restrict the energy range to half the
    # distance to the next lower and upper neighbor excitation, respectively.
    #-------------------------------------------------------------------------
    def restrict(self,erange=0.01,neigh=False,which=None): #arange=0.,prange=np.pi/4.1,
        #Assume that excitations are sorted by energy
        for iex, ex in enumerate(self.exlist):
            if isinstance(which,list):
                if not iex in which: continue
            lb = ex.energy*(1.-erange)
            rb = ex.energy*(1.+erange)
            if neigh:
                if iex >                  0: lb = max(lb,0.5*(ex.energy+self.exlist[iex-1].energy))
                if iex < len(self.exlist)-1: rb = min(rb,0.5*(ex.energy+self.exlist[iex+1].energy))
            ex.erange = [lb,rb]

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
        return f, np.array(w0), np.array(f0), l0

    #-------------------------------------------------------------------------
    # Plot excitations
    #-------------------------------------------------------------------------
    def plot(self,wb,dw,gamma,jarea=None,jcomp=None,jex=None,units="eV",fname=""):
        if units=="eV":
            wscal=13.605684958
        elif units=="Ry":
            wscal=1.
        else:
            err(1,"Unknown units "+units)
        w  = np.arange(wb[0],wb[1],dw)
        n  = len(w)
        f, w0, f0, l0 = self.lorentz(w,gamma,jarea,jcomp,jex)
        plt.title("Fitted Spectrum")
        plt.xlabel("Energy ["+units+"]")
        plt.ylabel("Oscillator Strength (Peak Height)")
        plt.plot(w*wscal,f)
        fmin = [0.]*len(f0)
        plt.vlines(w0*wscal,fmin,f0,colors="k")
        for iex in range(len(self.exlist)):
            plt.text(w0[iex]*wscal,f0[iex]+0.01*np.amax(f0),l0[iex])
        if fname=="":
            plt.show()
        else:
            plt.savefig(fname)

    #-------------------------------------------------------------------------
    # Plot all excitations into panels
    #-------------------------------------------------------------------------
    def plotPanels(self,wb,dw,gamma,jex=None,units="eV",fname=""):
        if units=="eV":
            wscal=13.605684958
        elif units=="Ry":
            wscal=1.
        else:
            err(1,"Unknown units "+units)
        w  = np.arange(wb[0],wb[1],dw)
        n  = len(w)
        fig, axs = plt.subplots(self.ncomp,self.narea,squeeze=False,sharex=True,sharey=True)
        fig.suptitle("Fitted Spectrum")
        for ax in axs.flat:
            ax.set(xlabel="Energy ["+units+"]", ylabel="Oscillator Strength (Peak Height)")
        for ax in axs.flat:# Hide x labels and tick labels for top plots and y ticks for right plots.
            ax.label_outer()
        for iarea in range(self.narea):
            for icomp in range(self.ncomp):
                f, w0, f0, l0 = self.lorentz(w,gamma,iarea,icomp,jex)
                axs[icomp][iarea].plot(w*wscal,f)
                fmin = [0.]*len(f0)
                axs[icomp][iarea].vlines(w0*wscal,fmin,f0,colors="k")
                for iex in range(len(self.exlist)):
                    axs[icomp][iarea].text(w0[iex]*wscal,f0[iex]+0.01*np.amax(f0),l0[iex])
        if fname=="":
            plt.show()
        else:
            plt.savefig(fname)
