import re
import numpy as np
from lmfit import Parameters
import matplotlib.pyplot as plt
#==============================================================================#
# Defines an excitation
#==============================================================================#
class Excitation:
    def __init__(self,ncalc=1,narea=1,ncomp=3,exdict={},name="S",energy=0.01,erange=0.01,phase=0.01,tmod=1.,dipoles=None,fix=False,ext=None):
        # Book keeping
        self.ncalc        = ncalc
        self.narea        = narea
        self.ncomp        = ncomp
        if not isinstance(dipoles,list):
            dipoles       = [[0.01]*self.ncomp for i in range(self.narea)]
        # Fundamental
        self.name         = exdict.get("name"        ,   name     )
        self.fix          = exdict.get("fix"         ,   fix      ) #Global fix (by configuration file)
        self.fixtmp       = exdict.get("fix"         ,   fix      ) #Temporary fix (that can be released and is used for generating fit parameters)
        # External
        self.ext          = ext
        # Fitting
        self.energy       = exdict.get("energy"      ,   energy   )
        self.energyErr    = exdict.get("energyErr"   ,   0.       )
        if isinstance(erange,float):
            erange        = [self.energy-erange,self.energy+erange]
        self.erange       =          exdict.get("erange"      ,   erange   ) #Energy range [min/max] (usually +/- 1% (i.e., x (1 +/- 0.01)) of the excitation energy)
        self.phase        =          exdict.get("phase"       ,   phase    )
        self.phaseErr     =          exdict.get("phaseErr"    ,   0.       )
        self.tmod         =          exdict.get("tmod"        ,   tmod     ) #T_eff = T*tmod (effective time for fitting)
        self.dipoles      = np.array(exdict.get("dipoles"     ,   dipoles  ))
        self.dipolesErr   = np.array(exdict.get("dipolesErr"  ,   np.zeros(self.dipoles.shape)))
        # Derived
        self.ampl         = np.array(exdict.get("ampl"        ,  np.zeros((self.ncalc,)+self.dipoles.shape))) #Note: tuples are concatenated. (n,) makes a tuple from the int n
        self.dipole       = np.array(exdict.get("dipole"      ,  [0.]*ncomp))
        self.dipoleErr    = np.array(exdict.get("dipoleErr"   ,  [0.]*ncomp))
        self.strength     =          exdict.get("strength"    ,   0.          )
        self.strengthErr  =          exdict.get("strengthErr" ,   0.          )
        self.strengths    = np.array(exdict.get("strengths"   ,  [0.]*narea)) #Oscillator strength equivalent derived from the area's dipole (does NOT sum up to the total oscillator strength)
        self.strengthsErr = np.array(exdict.get("strengthsErr",  [0.]*narea)) 
        #Significance
        self.signifFit    = exdict.get("signifFit"   ,   0.       ) #Significance: Can excit be replaced by other excitations
        self.signifErr    = exdict.get("signifErr"   ,   0.       ) #Significance: Does excitation's contribution to the spectrum match the fraction of its strength to the sum of all strengths
        self.signifAng    = exdict.get("signifAng"   ,   0.       ) #Significance: sqrt(angle between dipole and ext. excitation)
        self.signifExc    = exdict.get("signifExc"   ,   0.       ) #Significance: Similar to signifFit but re-fit dipole moment of single excitation alone
        self.signifRng    = exdict.get("signifRng"   ,   0.       ) #Significance: close to 1 if energy is close to the center of the energy range
        # Update derived (overwrite the latter)
        self.derived(errors=False)

        # Error checking
        if len(self.dipoles)!=self.narea or len(self.dipoles[0])!=self.ncomp: err(1,"Wrong shape of field 'dipoles'")
        if                                  len(self.dipole    )!=self.ncomp: err(1,"Wrong shape of field 'dipole'" )

    def todict(self):
        exdict = {}
        exdict["name"        ] = self.name
        exdict["fix"         ] = self.fix
        exdict["energy"      ] =    float(self.energy) #Explicit casting to float is necesary for ruaml.dump
        exdict["energyErr"   ] =    float(self.energyErr) #Explicit casting to float is necesary for ruaml.dump
        exdict["erange"      ] = [  float(self.erange[i]) for i in range(2)]
        exdict["phase"       ] =    float(self.phase)
        exdict["phaseErr"    ] =    float(self.phaseErr)
        exdict["tmod"        ] =    float(self.tmod)
        exdict["dipoles"     ] = [[ float(self.dipoles    [i][j])    for j in range(len(self.dipoles   [i]))] for i in range(len(self.dipoles   ))]
        exdict["dipolesErr"  ] = [[ float(self.dipolesErr [i][j])    for j in range(len(self.dipolesErr[i]))] for i in range(len(self.dipolesErr))]
        exdict["ampl"        ] = [[[float(self.ampl       [i][j][k]) for k in range(len(self.ampl[i][j]   ))] for j in range(len(self.ampl[i]   ))] for i in range(len(self.ampl))]
        exdict["dipole"      ] = [  float(self.dipole     [i])       for i in range(len(self.dipole       ))]
        exdict["dipoleErr"   ] = [  float(self.dipoleErr  [i])       for i in range(len(self.dipoleErr    ))]
        exdict["strength"    ] =    float(self.strength)
        exdict["strengthErr" ] =    float(self.strengthErr)
        exdict["strengths"   ] = [  float(self.strengths   [i])     for i in range(len(self.strengths    ))]
        exdict["strengthsErr"] = [  float(self.strengthsErr[i])     for i in range(len(self.strengthsErr ))]
        exdict["signifFit"   ] =    float(self.signifFit)
        exdict["signifErr"   ] =    float(self.signifErr)
        exdict["signifAng"   ] =    float(self.signifAng)
        exdict["signifExc"   ] =    float(self.signifExc)
        exdict["signifRng"   ] =    float(self.signifRng)
        return exdict

    # Updates derived components
    def derived(self,errors=True):
        #Define constants in Rydberg units
        m    = 0.5
        e2   = 2. #e^2
        hbar = 1.
        self.dipole       = np.array([sum([self.dipoles[iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)])
        self.strength     =  2.*m/(3.*e2*hbar)*self.energy*np.linalg.norm(self.dipole        )**2
        self.strengths    = [2.*m/(3.*e2*hbar)*self.energy*np.linalg.norm(self.dipoles[iarea])**2      for iarea in range(self.narea)]
        Hw = self.ext.getVal([self.energy])[0]
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                self.ampl[icalc][iarea] = 1./hbar * self.ext.efield * np.dot(self.ext.epol[icalc],self.dipole) * np.abs(Hw) * self.dipoles[iarea]
        if errors:
            self.dipoleErr    = np.array([np.sqrt(sum([self.dipolesErr[iarea][icomp]**2 for iarea in range(self.narea)])) for icomp in range(self.ncomp)])
            self.strengthErr  =  2.*m/(3.*e2*hbar)*np.sqrt((self.energyErr*self.energy*np.linalg.norm(self.dipole        )**2)**2 + (2.*self.energy*np.dot(self.dipole        ,self.dipoleErr        ))**2)
            self.strengthsErr = [2.*m/(3.*e2*hbar)*np.sqrt((self.energyErr*self.energy*np.linalg.norm(self.dipoles[iarea])**2)**2 + (2.*self.energy*np.dot(self.dipoles[iarea],self.dipolesErr[iarea]))**2) for iarea in range(self.narea)]

    # Temporarily fixes the excitation
    def fixMe(self):
        self.fixtmp = True

    # Release temporary fix
    def releaseMe(self):
        self.fixtmp = self.fix

    # Reset erange
    def resetErange(self,erange=0.01):
        self.erange = np.array([self.energy-erange,self.energy+erange])

    # Set significance
    def setSignificance(self,signifFit,signifErr,signifAng,signifExc,signifRng):
        self.signifFit = signifFit
        self.signifErr = signifErr
        self.signifAng = signifAng
        self.signifExc = signifExc
        self.signifRng = signifRng

#============================================================================#
# Defines a list of excitations
#============================================================================#
class Excitations:
    def __init__(self,ncalc,narea,ncomp,ext,exlist=None):
        self.exlist = []
        self.ncalc  = ncalc
        self.narea  = narea
        self.ncomp  = ncomp
        self.ext    = ext
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
        ex = Excitation(ncalc=self.ncalc,narea=self.narea,ncomp=self.ncomp,exdict=exdict,name=name,energy=energy,erange=erange,phase=phase,tmod=tmod,dipoles=dipoles,fix=fix,ext=self.ext)
        self.exlist.append(ex)
        if sort: self.sort()
        #return self.exlist.index(ex)

    #-------------------------------------------------------------------------
    # Copy excitations object (potentially remove excitations)
    #-------------------------------------------------------------------------
    def copy(self):
        return Excitations(self.ncalc,self.narea,self.ncomp,self.ext,exlist=self.exlist)

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
            ex.phase    = 0.
            ex.phaseErr = 0.

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
    def updateFromParam(self,params,errors=True):
        for iex, ex in enumerate(self.exlist):
            if ex.fixtmp: continue #Continue if values of this excitation need no update
            ex.energy    = params[f"w{iex}"].value
            if errors and isinstance(params[f"w{iex}"].stderr,float):
                ex.energyErr = params[f"w{iex}"].stderr
            ex.phase     = params[f"p{iex}"].value
            if errors and isinstance(params[f"p{iex}"].stderr,float):
                ex.phaseErr  = params[f"p{iex}"].stderr
            ex.tmod      = params[f"t{iex}"].value
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    ex.dipoles   [iarea][icomp] = params[f"d{iex}_{iarea}_{icomp}"].value
                    if errors and isinstance(params[f"d{iex}_{iarea}_{icomp}"].stderr,float):
                        ex.dipolesErr[iarea][icomp] = params[f"d{iex}_{iarea}_{icomp}"].stderr
            ex.derived(errors=errors)

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
    # Reset erange
    #-------------------------------------------------------------------------
    def resetErange(self,erange=0.01):
        for ex in self.exlist:
            ex.resetErange(erange=erange)

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
            ex.derived(errors=False)
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
        we = []
        f0 = []
        fe = []
        l0 = []
        for iex, ex in enumerate(self.exlist):
            if isinstance(jex,list):
                if not jex in iex: cycle
            ex.derived()
            w0.append(ex.energy)
            we.append(ex.energyErr)
            if isinstance(jarea,int) and isinstance(jcomp,int):
                f0.append(ex.dipoles   [jarea][jcomp])
                fe.append(ex.dipolesErr[jarea][jcomp])
            elif isinstance(jarea,int):
                f0.append(ex.strengths   [jarea])
                fe.append(ex.strengthsErr[jarea])
            elif isinstance(jcomp,int):
                f0.append(ex.dipole   [jcomp])
                fe.append(ex.dipoleErr[jcomp])
            else:
                f0.append(ex.strength   )
                fe.append(ex.strengthErr)
            l0.append(ex.name)
            for i in range(n):
                gw  = gamma*w0[-1]
                gw2 = gw*gw
                dw  = w[i]-w0[-1]
                dw2 = dw*dw
                f[i] += f0[-1]*gw2/(dw2+gw2)
        return f, np.array(w0), np.array(we), np.array(f0), np.array(fe), l0

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
        f, w0, dw, f0, df, l0 = self.lorentz(w,gamma,jarea,jcomp,jex)
        plt.title("Fitted Spectrum")
        plt.xlabel("Energy ["+units+"]")
        plt.ylabel("Oscillator Strength (Peak Height)")
        plt.plot(w*wscal,f)
        fmin = [0.]*len(f0)
        plt.vlines  (w0*wscal,fmin,f0                      ,colors="k")
        plt.errorbar(w0*wscal     ,f0,xerr=dw*wscal,yerr=df,color ="k",markersize=2.,ecolor ="r",fmt="o",capsize=3.)
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
                f, w0, dw, f0, df, l0 = self.lorentz(w,gamma,iarea,icomp,jex)
                axs[icomp][iarea].plot(w*wscal,f)
                fmin = [0.]*len(f0)
                axs[icomp][iarea].vlines  (w0*wscal,fmin,f0                      ,colors="k")
                axs[icomp][iarea].errorbar(w0*wscal     ,f0,xerr=dw*wscal,yerr=df,color ="k",markersize=2.,ecolor ="r",fmt="o",capsize=3.)
                for iex in range(len(self.exlist)):
                    axs[icomp][iarea].text(w0[iex]*wscal,f0[iex]+0.01*np.amax(f0),l0[iex])
        if fname=="":
            plt.show()
        else:
            plt.savefig(fname)
