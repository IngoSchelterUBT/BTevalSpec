import re
import numpy as np
#==============================================================================#
# Defines an excitation
#==============================================================================#
class Excitation:
    def __init__(self,ncalc=1,narea=1,ncomp=3,exdict={},name="S",energy="0.",phase="0.01",dipole=0.01):
        # Book keeping
        self.ncalc    = ncalc
        self.narea    = narea
        self.ncomp    = ncomp
        # Fundamental
        self.name     = exdict.get("name"    ,   name                                                     )
        self.fix      = exdict.get("fix"     ,   False                                                    ) 
        self.energy   = exdict.get("energy"  ,   energy                                                   )
        self.phase    = exdict.get("phase"   ,   phase                                                    )
        self.dipoles  = exdict.get("dipoles" , [[dipole]*ncomp   for i in range(narea)]                  )
#        # Derived
        self.dipole   = exdict.get("dipole"  ,  [0.]*ncomp                                                )
        self.strength = exdict.get("strength",   0.                                                       )
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
        self.dipole   = np.array([sum(np.transpose(self.dipoles)[i]) for i in range(self.ncomp)])
        self.strength = 2.*m*self.energy*np.linalg.norm(self.dipole)**2/(3.*e2*hbar)

#============================================================================#
# Defines a list of excitations
#============================================================================#
class Excitations:
    def __init__(self,ncalc,narea,ncomp,exlist=None):
        self.exlist = []
        self.ncalc  = ncalc
        self.narea  = narea
        self.ncomp  = ncomp
        try:
            for ex in exlist:
                self.add(exdict=ex,sort=False)
            self.sort()
            #Todo: Check if excitations in exlist match ncalc, narea, ncomp
        except:
            pass

    #-------------------------------------------------------------------------
    # Add an excitation to the list
    #-------------------------------------------------------------------------
    def add(self,exdict=None,name="S",energy=0.,phase=0.,dipole=0.1,sort=True):
        #Todo: Check if new excitation matches ncalc, narea, ncomp
        if isinstance(exdict,dict):
            self.exlist.append(Excitation(ncalc=self.ncalc,narea=self.narea,ncomp=self.ncomp,exdict=exdict))
        else:
            self.exlist.append(Excitation(ncalc=self.ncalc,narea=self.narea,ncomp=self.ncomp,name=name,energy=energy,phase=phase,dipole=dipole))
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
