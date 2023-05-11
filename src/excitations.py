import re
import numpy as np
from lmfit import Parameters
#==============================================================================#
# Defines an excitation
#==============================================================================#
class Excitation:
    def __init__(self,narea=1,ncomp=3,exdict={},name="S",energy="0.01",phase="0.01",dipoles=None,fix=False):
        # Book keeping
        self.narea    = narea
        self.ncomp    = ncomp
        if not isinstance(dipoles,list):
            dipoles = [[0.01]*self.ncomp for i in range(self.narea)]
        # Fundamental
        self.name     = exdict.get("name"    ,   name     )
        self.fix      = exdict.get("fix"     ,   fix      ) #Global fix (by configuration file)
        self.fixtmp   = exdict.get("fix"     ,   fix      ) #Temporary fix (that can be released and is used for generating fit parameters)
        self.energy   = exdict.get("energy"  ,   energy   )
        self.phase    = exdict.get("phase"   ,   phase    )
        self.dipoles  = exdict.get("dipoles" ,   dipoles  )
        # Derived
        self.dipole   = exdict.get("dipole"  ,  [0.]*ncomp)
        self.strength = exdict.get("strength",   0.       )
        # Update derived (overwrite the latter)
        self.derived()

        # Error checking
        if len(self.dipoles)!=self.narea or len(self.dipoles[0])!=self.ncomp: err(1,"Wrong shape of field 'dipoles'")
        if                                  len(self.dipole    )!=self.ncomp: err(1,"Wrong shape of field 'dipole'" )

    def todict(self):
        exdict = {}
        exdict["name"    ] = self.name
        exdict["fix"     ] = self.fix
        exdict["energy"  ] = float(self.energy) #Explicit casting to float is necesary for ruaml.dump
        exdict["phase"   ] = float(self.phase)
        exdict["dipoles" ] = [[float(self.dipoles[i][j]) for j in range(len(self.dipoles[i]))] for i in range(len(self.dipoles))]
        exdict["dipole"  ] = [ float(self.dipole[i])     for i in range(len(self.dipole))]
        exdict["strength"] = float(self.strength)
        return exdict

    # Updates derived components
    def derived(self):
        #Define constants in Rydberg units
        m    = 0.5
        e2   = 2. #e^2
        hbar = 1.
        self.dipole   = np.array([sum([self.dipoles[iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)])
        self.strength = 2.*m*self.energy*np.linalg.norm(self.dipole)**2/(3.*e2*hbar)

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
    def add(self,exdict=None,name="S",energy=0.,phase=0.,dipoles=None,fix=False,sort=True):
        #Todo: Check if new excitation matches narea, ncomp
        if isinstance(exdict,dict):
            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,exdict=exdict))
        else:
            self.exlist.append(Excitation(narea=self.narea,ncomp=self.ncomp,name=name,energy=energy,phase=phase,dipoles=dipoles,fix=fix))
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
    # Generate lm fit parameters
    # Uses the fixtmp variable instead of fix
    #-------------------------------------------------------------------------
    def toParams(self):
        params = Parameters()
        for iex, ex in enumerate(self.exlist):
            params.add(f"w{iex}",vary=not ex.fixtmp,value=ex.energy)
            params.add(f"p{iex}",vary=not ex.fixtmp,value=ex.phase)
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
    # Print excitations
    #-------------------------------------------------------------------------
    def print(self):
        for iex, ex in enumerate(self.exlist):
            ex.derived()
            print("Ex ",iex)
            print("  Name:     ", ex.name    )
            print("  Energy:   ", ex.energy  )
            print("  Phase:    ", ex.phase   )
            print("  Strength: ", ex.strength)
            print("  Dipole:   ", ex.dipole  )
            print("  Dipoles per area:")
            for iarea in range(self.narea):
                print("      ", iarea, ex.dipoles[iarea])

