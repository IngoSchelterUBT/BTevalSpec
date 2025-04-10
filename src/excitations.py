import re
import numpy as np
from lmfit import Parameters
import extern
import matplotlib.pyplot as plt
#==============================================================================#
# Defines an excitation
#==============================================================================#
class Excitation:
    def __init__(self,ncalc=1,narea=1,ncomp=3,exdict={},name="S",energy=0.01,erange=0.01,phase=0.01,tmod=1.,dipoles=None,fix=False,ext=None):
        # Book keeping
        self.ncalc           = ncalc
        self.narea           = narea
        self.ncomp           = ncomp
        if not isinstance(dipoles,list):
            dipoles          = [[0.01]*self.ncomp for i in range(self.narea)]
        # Fundamental
        self.name            = exdict.get("name"        ,   name     )
        self.fix             = exdict.get("fix"         ,   fix      ) #Global fix (by configuration file)
        self.fixtmp          = exdict.get("fix"         ,   fix      ) #Temporary fix (that can be released and is used for generating fit parameters)
        # External
        self.ext             = ext
        # Fitting
        self.energy          = exdict.get("energy"      ,   energy   )
        self.energyErr       = exdict.get("energyErr"   ,   0.       )
        if isinstance(erange,float):
            erange           = [self.energy-erange,self.energy+erange]
        self.erange          =          exdict.get("erange"         ,   erange   ) #Energy range [min/max] (usually +/- 1% (i.e., x (1 +/- 0.01)) of the excitation energy)
        self.phase           =          exdict.get("phase"          ,   phase    )
        self.phaseErr        =          exdict.get("phaseErr"       ,   0.       )
        self.tmod            =          exdict.get("tmod"           ,   tmod     ) #T_eff = T*tmod (effective time for fitting)
        self.dipoles         = np.array(exdict.get("dipoles"        ,   dipoles  ))
        self.dipolesErr      = np.array(exdict.get("dipolesErr"     ,   np.zeros(self.dipoles.shape)))
        # Derived
        self.ampl            = np.array(exdict.get("ampl"           ,  np.zeros((self.ncalc,)+self.dipoles.shape))) #Note: tuples are concatenated. (n,) makes a tuple from the int n
        self.amplErr         = np.array(exdict.get("amplErr"        ,  np.zeros((self.ncalc,)+self.dipoles.shape))) #Note: -"-
        self.dipole          = np.array(exdict.get("dipole"         ,  [0.]*self.ncomp))
        self.dipoleErr       = np.array(exdict.get("dipoleErr"      ,  [0.]*self.ncomp))
        self.transdensDipole = np.array(exdict.get("transdensDipole",  [0.]*self.ncomp))
        self.strength        =          exdict.get("strength"       ,   0.          )
        self.strengthErr     =          exdict.get("strengthErr"    ,   0.          )
        self.strengthEped    = np.array(exdict.get("strengthEped"   ,  [0.]*self.ncalc)) #Oscillator strength times eped => Projection of oscillator strength on external field polarization direction
        self.strengthEpedErr = np.array(exdict.get("strengthEpedErr",  [0.]*self.ncalc))
        self.strengths       = np.array(exdict.get("strengths"      ,  [0.]*self.narea)) #Oscillator strength equivalent derived from the area's dipole (does NOT sum up to the total oscillator strength)
        self.strengthsErr    = np.array(exdict.get("strengthsErr"   ,  [0.]*self.narea)) 
        self.eped            = np.array(exdict.get("eped"           ,  [0.]*self.ncalc)) #Angle factor <e_dip,e_pol> = cos(angle between excitation's dipole and external polarization vector)
        self.epedErr         = np.array(exdict.get("epedErr"        ,  [0.]*self.ncalc))
        #Significance
        self.signifFit       =          exdict.get("signifFit"      ,   0.       ) #Significance: Can excit be replaced by other excitations
        self.signifErr       =          exdict.get("signifErr"      ,   0.       ) #Significance: Does excitation's contribution to the spectrum match the fraction of its strength to the sum of all strengths
        self.signifAng       =          exdict.get("signifAng"      ,   0.       ) #Significance: sqrt(angle between dipole and ext. excitation)
        self.signifExc       =          exdict.get("signifExc"      ,   0.       ) #Significance: Similar to signifFit but re-fit dipole moment of single excitation alone
        self.signifRng       =          exdict.get("signifRng"      ,   0.       ) #Significance: close to 1 if energy is close to the center of the energy range
        self.signifPha       =          exdict.get("signifPha"      ,   0.       ) #Significance: close to 1 if energy is close to the center of the energy range
        # Update derived (overwrite the latter)
        self.derived(errors=False)

        # Error checking
        if len(self.dipoles)!=self.narea or len(self.dipoles[0])!=self.ncomp: err(1,"Wrong shape of field 'dipoles'")
        if                                  len(self.dipole    )!=self.ncomp: err(1,"Wrong shape of field 'dipole'" )

    def todict(self):
        exdict = {}
        exdict["name"           ] = self.name
        exdict["fix"            ] = self.fix
        exdict["energy"         ] =    float(self.energy) #Explicit casting to float is necesary for ruaml.dump
        exdict["energyErr"      ] =    float(self.energyErr) #Explicit casting to float is necesary for ruaml.dump
        exdict["erange"         ] = [  float(self.erange[i]) for i in range(2)]
        exdict["phase"          ] =    float(self.phase)
        exdict["phaseErr"       ] =    float(self.phaseErr)
        exdict["tmod"           ] =    float(self.tmod)
        exdict["dipoles"        ] = [[ float(self.dipoles        [i][j])    for j in range(len(self.dipoles     [i]))] for i in range(len(self.dipoles   ))]
        exdict["dipolesErr"     ] = [[ float(self.dipolesErr     [i][j])    for j in range(len(self.dipolesErr  [i]))] for i in range(len(self.dipolesErr))]
        exdict["ampl"           ] = [[[float(self.ampl           [i][j][k]) for k in range(len(self.ampl     [i][j]))] for j in range(len(self.ampl   [i]))] for i in range(len(self.ampl   ))]
        exdict["amplErr"        ] = [[[float(self.amplErr        [i][j][k]) for k in range(len(self.amplErr  [i][j]))] for j in range(len(self.amplErr[i]))] for i in range(len(self.amplErr))]
        exdict["dipole"         ] = [  float(self.dipole         [i])       for i in range(len(self.dipole         ))]
        exdict["dipoleErr"      ] = [  float(self.dipoleErr      [i])       for i in range(len(self.dipoleErr      ))]
        exdict["transdensDipole"] = [  float(self.transdensDipole[i])       for i in range(len(self.transdensDipole))]
        exdict["strength"       ] =    float(self.strength)
        exdict["strengthErr"    ] =    float(self.strengthErr)
        exdict["strengthEped"   ] = [  float(self.strengthEped[i])         for i in range(len(self.strengthEped   ))]
        exdict["strengthErrEped"] = [  float(self.strengthEpedErr[i])      for i in range(len(self.strengthEpedErr))]
        exdict["strengths"      ] = [  float(self.strengths   [i])      for i in range(len(self.strengths      ))]
        exdict["strengthsErr"   ] = [  float(self.strengthsErr[i])      for i in range(len(self.strengthsErr   ))]
        exdict["eped"           ] = [  float(self.eped[i])              for i in range(len(self.eped           ))]
        exdict["epedErr"        ] = [  float(self.epedErr[i])           for i in range(len(self.epedErr        ))]
        exdict["signifFit"      ] =    float(self.signifFit)
        exdict["signifErr"      ] =    float(self.signifErr)
        exdict["signifAng"      ] =    float(self.signifAng)
        exdict["signifExc"      ] =    float(self.signifExc)
        exdict["signifRng"      ] =    float(self.signifRng)
        exdict["signifPha"      ] =    float(self.signifPha)
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
            epol                     = self.ext.epol[icalc]/np.linalg.norm(self.ext.epol[icalc])
            self.eped        [icalc] = np.dot(epol,self.dipole)/np.linalg.norm(self.dipole)
            self.strengthEped[icalc] = self.strength * self.eped[icalc]
            for iarea in range(self.narea):
                self.ampl   [icalc][iarea] = -1./hbar * self.ext.efield * np.abs(Hw) *          np.dot(self.ext.epol[icalc],self.dipole   ) * self.dipoles[iarea]
                self.amplErr[icalc][iarea] =  1./hbar * self.ext.efield * np.abs(Hw) * np.sqrt((np.dot(self.ext.epol[icalc],self.dipoleErr) * self.dipoles[iarea])**2 + (np.dot(self.ext.epol[icalc],self.dipole) * self.dipolesErr[iarea])**2)
        if errors:
            self.dipoleErr    = np.array([np.sqrt(sum([self.dipolesErr[iarea][icomp]**2 for iarea in range(self.narea)])) for icomp in range(self.ncomp)])
            self.strengthErr  =  2.*m/(3.*e2*hbar)*np.sqrt((self.energyErr*self.energy*np.linalg.norm(self.dipole        )**2)**2 + (2.*self.energy*np.dot(self.dipole        ,self.dipoleErr        ))**2)
            self.strengthsErr = [2.*m/(3.*e2*hbar)*np.sqrt((self.energyErr*self.energy*np.linalg.norm(self.dipoles[iarea])**2)**2 + (2.*self.energy*np.dot(self.dipoles[iarea],self.dipolesErr[iarea]))**2) for iarea in range(self.narea)]
            tmp2 = np.sqrt(np.abs(np.dot(self.dipole,self.dipoleErr)))
            tmp  = [np.sqrt((self.dipoleErr[i]/np.linalg.norm(self.dipole))**2 + (self.dipole[i]*tmp2/np.linalg.norm(self.dipole)**2)**2) for i in range(self.ncomp)]
            for icalc in range(self.ncalc):
                epol              = self.ext.epol[icalc]/np.linalg.norm(self.ext.epol[icalc])
                self.epedErr        [icalc] = np.sqrt(sum([(tmp[i]*epol[i])**2 for i in range(self.ncomp)]))
                self.strengthEpedErr[icalc] = np.sqrt((self.strength*self.epedErr[icalc])**2 + (self.strengthErr*self.eped[icalc])**2)

    # Set dipoles from transition densities
    def setTransdensDipole(self,tdd):
        self.transdensDipole = tdd

    # Temporarily fixes the excitation
    def fixMe(self,permanent=False):
        self.fixtmp = True
        if permanent: self.fix = True

    # Release temporary fix
    def releaseMe(self,permanent=False):
        if permanent: self.fix = False
        self.fixtmp = self.fix

    # Reset erange
    def resetErange(self,erange=0.01):
        self.erange = np.array([self.energy-erange,self.energy+erange])

    # Set significance
    def setSignificance(self,signifFit,signifErr,signifAng,signifExc,signifRng,signifPha):
        self.signifFit = signifFit
        self.signifErr = signifErr
        self.signifAng = signifAng
        self.signifExc = signifExc
        self.signifRng = signifRng
        self.signifPha = signifPha

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
    # Set dipoles from transition densities
    #-------------------------------------------------------------------------
    def setTransdensDipoles(self,tdd):
        for iex, ex in enumerate(self.exlist):
            ex.setTransdensDipole(tdd[iex])

    #-------------------------------------------------------------------------
    # Copy excitations object (potentially remove excitations)
    #-------------------------------------------------------------------------
    def copy(self):
        return Excitations(self.ncalc,self.narea,self.ncomp,self.ext,exlist=self.exlist)

    #-------------------------------------------------------------------------
    # Remove excitation by index or name
    #-------------------------------------------------------------------------
    def remove(self,rmidx=[],rmname=[],dbg=0):
        # Collect excitations to remove
        #Don't remove yet since this would change the original indexing
        rmex = []
        for iex, ex in enumerate(self.exlist):
            if iex     in rmidx or ex.name in rmname:  rmex.append(ex)
        # Remove excitations
        for ex in rmex:
            if dbg>0: print(f"  Remove {ex.name}")
            self.exlist.remove(ex)

    #-------------------------------------------------------------------------
    # Sorts the list of excitations and rename those with generic names "S[1234..]"
    #-------------------------------------------------------------------------
    def sort(self):
        self.exlist = sorted(self.exlist, key=lambda x: x.energy)
        pattern     = re.compile("^S[0-9]*$") #Matches name "S<number>" or "S"
        for iex, ex in enumerate(self.exlist):
            if pattern.match(ex.name):
                ex.name = "S"+str(iex+1)

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
    def toParams(self,ext,noEnergy=False,noPhase=False,noTmod=False):
        params = Parameters()
        for iex, ex in enumerate(self.exlist):
            params.add(f"w{iex}",vary=not ex.fixtmp and not noEnergy,value=ex.energy,min=ex.erange[0],max=ex.erange[1])
            if not noPhase: #Otherwise, the phase is taken from the external-field profile
                #Phase boundaries are (i) not necessary and (ii) can spoil the fit
                #phasemm = np.angle(ext.getVal(ex.erange))%(2*np.pi)
                #if phasemm[0]>phasemm[1]: phasemm[1] += 2.*np.pi
                params.add(f"p{iex}",vary=not ex.fixtmp,value=ex.phase)#,min=phasemm[0],max=phasemm[1])
            params.add(f"t{iex}",vary=not ex.fixtmp and not noTmod  ,value=ex.tmod  )
            for iarea in range(self.narea):
                for icomp in range(self.ncomp):
                    params.add(f"d{iex}_{iarea}_{icomp}",vary=not ex.fixtmp,value=ex.dipoles[iarea][icomp])
        return params

    #-------------------------------------------------------------------------
    # Update from fit parameters
    #-------------------------------------------------------------------------
    def updateFromParam(self,params,noPhase,ext,errors=True):
        for iex, ex in enumerate(self.exlist):
            if ex.fixtmp: continue #Continue if values of this excitation need no update
            ex.energy    = params[f"w{iex}"].value
            if errors and isinstance(params[f"w{iex}"].stderr,float):
                ex.energyErr = params[f"w{iex}"].stderr
            if noPhase:
                ex.phase     = np.angle(ext.getVal([ex.energy])[0])%(2*np.pi)
                ex.phaseErr  = 0.
            else:
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
    def fix(self,whichidx=None,whichname=None,dbg=0,inverse=False,permanent=False):
        for iex, ex in enumerate(self.exlist):
            if inverse:
                if isinstance(whichidx,list):
                    if iex in whichidx: continue
                if isinstance(whichname,list):
                    if ex.name in whichname: continue
            else:
                if isinstance(whichidx,list):
                    if not iex in whichidx: continue
                if isinstance(whichname,list):
                    if not (ex.name in whichname or "all" in whichname): continue
            ex.fixMe(permanent=permanent)
            if dbg>0: print(f"Fix excitation {ex.name}")

    #-------------------------------------------------------------------------
    # Temporarily fix excitations in a given interval (or the inverse interval)
    #-------------------------------------------------------------------------
    def fixErange(self,erange,inverse=False,permanent=False):
        for iex, ex in enumerate(self.exlist):
            if ex.energy < erange[1] and ex.energy >erange[0]: # inside inverval
                if not inverse: ex.fixMe(permanent=permanent)  # fix if not inverse
            else:                                              #outside inverval
                if     inverse: ex.fixMe(permanent=permanent)  # fix if     inverse

    #-------------------------------------------------------------------------
    # Release all or specific excitations
    #-------------------------------------------------------------------------
    def release(self,whichidx=None,whichname=None,dbg=0,inverse=False,permanent=False):
        for iex, ex in enumerate(self.exlist):
            if inverse:
                if isinstance(whichidx,list):
                    if iex in whichidx: continue
                if isinstance(whichname,list):
                    if ex.name in whichname: continue
            else:
                if isinstance(whichidx,list):
                    if not iex in whichidx: continue
                if isinstance(whichname,list):
                    if not (ex.name in whichname or "all" in whichname): continue
            ex.releaseMe(permanent=permanent)
            if dbg>0: print(f"Release excitation {ex.name}")

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
    def print(self,long=False):
        print("")
        print("Excitations")
        for iex, ex in enumerate(self.exlist):
            ex.derived(errors=False)
            if long:
                print("Ex ",iex)
                print("  Name:     ", ex.name    )
                print("  Fix:      ", ex.fix     )
                print("  Energy:   ", ex.energy  )
                print("  Erange:   ", ex.erange  )
                print("  Phase:    ", ex.phase   )
                print("  Tmod:     ", ex.tmod    )
                print("  Strength: ", ex.strength)
                print("  Dipole:   ", ex.dipole  )
                print("  Dipoles per area:")
                print("")
                for iarea in range(self.narea):
                    print("      ", iarea, ex.dipoles[iarea])
            else:
                pistr=u"\u03c0"
                print(f"{ex.name:5s} | fix:{ex.fix!s:^5} | E={ex.energy:7.4f}Ry | phi={ex.phase/np.pi:5.2f}{pistr} | f={ex.strength:6.3f}",end="")
                for icalc in range(self.ncalc):
                    ampl = np.linalg.norm([sum([ex.ampl[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)])
                    print(f"| |a|({icalc})={ampl:6.3e}",end="")
                print("")

    #-------------------------------------------------------------------------
    # Create amplitude files
    #-------------------------------------------------------------------------
    def excitFiles(self,tprop):

        for icalc in range(self.ncalc):
            with open(f"excit_{icalc+1:1d}.dat","w") as fh:
                fh.write("# name  | energy[Ry] |  strength  | phase[pi]  |  eped      |strengthEped|   energyErr  | strengthError|   phaseErr   |   epedErr   |strengthEpedErr|signifFit|signifAng|signifExc|signifErr|signifRng|signifPha\n")
                for iex, ex in enumerate(self.exlist):
                    pha = (ex.phase/np.pi)%(2.*np.pi)
                    if pha >   np.pi: pha -= 2.*np.pi
                    if pha <= -np.pi: pha += 2.*np.pi
                    fh.write(f"{ex.name:8s} {ex.energy:12.5f} {ex.strength:12.5f} {pha:12.5f} {ex.eped[icalc]:12.5f} {ex.strengthEped[icalc]:12.5f} {ex.energyErr:14.7f} {ex.strengthErr:14.7f} {ex.phaseErr/np.pi:14.7f} {ex.epedErr[icalc]:14.7f} {ex.strengthEpedErr[icalc]:14.7f} {ex.signifFit:9.2f} {ex.signifAng:9.2f} {ex.signifExc:9.2f} {ex.signifErr:9.2f} {ex.signifRng:9.2f} {ex.signifPha:9.2f}\n")

        xyz = ["x","y","z"]
        for icalc in range(self.ncalc):
            for iarea in range(self.narea):
                with open(f"ampl_{icalc+1:1d}_{str(iarea+1).zfill(2)}.dat","w") as fh:
                    fh.write("# name  |   energy   |   ampl x   |   ampl y   |   ampl z    |  abs(ampl)  | energy err | ampl err x | ampl err y | ampl err z | abs(ampl) err\n")
                    for iex, ex in enumerate(self.exlist):
                        ampl    = np.linalg.norm([ex.ampl   [icalc][iarea][icomp] for icomp in range(self.ncomp)])
                        amplErr = np.linalg.norm([ex.amplErr[icalc][iarea][icomp] for icomp in range(self.ncomp)])
                        fh.write(f"{ex.name:8s} {ex.energy:12.5f} {ex.ampl[icalc][iarea][0]:12.5e} {ex.ampl[icalc][iarea][1]:12.5e} {ex.ampl[icalc][iarea][2]:12.5e} {ampl:12.5e} {ex.energyErr:12.5f} {ex.amplErr[icalc][iarea][0]:12.5e} {ex.amplErr[icalc][iarea][1]:12.5e} {ex.amplErr[icalc][iarea][2]:12.5e} {amplErr:12.5e}\n")
            with open(f"ampl_{icalc+1:1d}_glob.dat","w") as fh:
                fh.write("# name  |   energy   |   ampl x   |   ampl y   |   ampl z    |  abs(ampl)  | energy err | ampl err x | ampl err y | ampl err z | abs(ampl) err\n")
                for iex, ex in enumerate(self.exlist):
                    amplGlob    = [           sum([ex.ampl   [icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)]
                    amplGlobErr = [np.linalg.norm([ex.amplErr[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)]
                    ampl    = np.linalg.norm(amplGlob   )
                    amplErr = np.linalg.norm(amplGlobErr)
                    fh.write(f"{ex.name:8s} {ex.energy:12.5f} {amplGlob[0]:12.5e} {amplGlob[1]:12.5e} {amplGlob[2]:12.5e} {ampl:12.5e} {ex.energyErr:12.5f} {amplGlobErr[0]:12.5e} {amplGlobErr[1]:12.5e} {amplGlobErr[2]:12.5e} {amplErr:12.5e}\n")

    #-------------------------------------------------------------------------
    # Return Lorentz function as well as energies, strengths, and labels
    #  - jex is a list of excitation indices. If None: Take all indices
    #  - jarea is an integer specifying the area. If None: Sum up all areas
    #  - jcomp is 0 (x), 1 (y), or 2 (z). if None: Use an oscillator-strength equivalent derived from the area's dipole
    #  - if jcalc, jarea, and jcomp are None: Use the oscillator strengths as height
    #  - if jarea and jcomp are None but jcalc is an integer, use the absolute amplitude as height (error is still zero)
    #-------------------------------------------------------------------------
    def lorenz(self,w,gamma,jcalc=None,jarea=None,jcomp=None,jex=None):
        n  = len(w)
        f  = np.zeros(n,dtype=float)
        w0 = []
        we = []
        f0 = []
        fe = []
        l0 = []
        pattern = re.compile("^S[0-9]*$") #Matches name "S<number>" or "S"
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
                if isinstance(jcalc,int):
#                    f0.append(ex.strengthEped   [jcalc])
#                    fe.append(ex.strengthEpedErr[jcalc])
                    ampl = np.linalg.norm([sum([ex.ampl[jcalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)])
                    f0.append(ampl)
                    fe.append(0.)
                else:
                    f0.append(ex.strength   )
                    fe.append(ex.strengthErr)
            if pattern.match(ex.name):
                l0.append(u"\u03C3"+str(iex+1))
            else:
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
    def plot(self,wb=[None,None],dw=0.00001,gamma=np.pi/1000.,jcalc=None,jarea=None,jcomp=None,jex=None,units="eV",fname=""):
        if units=="eV":
            wscal=13.605684958
        elif units=="Ry":
            wscal=1.
        else:
            err(1,"Unknown units "+units)
        wrange = self.exlist[-1].energy-self.exlist[0].energy
        wb0 = self.exlist[ 0].energy - 0.1*wrange
        wb1 = self.exlist[-1].energy + 0.1*wrange
        if wb[0]: wb0 = wb[0]
        if wb[1]: wb1 = wb[1]
        w  = np.arange(wb0,wb1,dw)
        n  = len(w)
        f, w0, dw, f0, df, l0 = self.lorenz(w,gamma,jcalc=jcalc,jarea=jarea,jcomp=jcomp,jex=jex)
        plt.clf()
        #plt.title("Fitted Spectrum")
        plt.xlabel("Energy ["+units+"]")
        if isinstance(jcalc,int):
            plt.ylabel("absolute amplitude |a_0j| (Peak Height)")
        else:
            plt.ylabel("isotropic oscillator strength f_0j (Peak Height)")
        plt.plot(w*wscal,f,color="#C0C0C3")
        fmin = [0.]*len(f0)
        plt.vlines  (w0*wscal,fmin,f0                      ,colors="#36454F")
        plt.errorbar(w0*wscal     ,f0,xerr=dw*wscal,yerr=df,color ="#B22222",markersize=2.,ecolor ="r",fmt="o",capsize=3.)
        for iex in range(len(self.exlist)):
            plt.text(w0[iex]*wscal,f0[iex]+0.01*np.amax(f0),l0[iex])
        if fname=="":
            plt.show()
        else:
            plt.savefig(fname,dpi=300)

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
                f, w0, dw, f0, df, l0 = self.lorenz(w,gamma,jarea=iarea,jcomp=icomp,jex=jex)
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

    #-------------------------------------------------------------------------
    # Output tex table
    #-------------------------------------------------------------------------
    def latexTable(self,fname=""):
        if fname=="": fname="excit_LatexTable.dat"
        ev=13.605684958
        lines = []
        lines.append(fr"\begin{{table}}[ht]")
        lines.append(fr"\centering")
        lines.append(fr"\fontsize{{8}}{{10}}\selectfont")
        lines.append(fr"\begin{{tabular}}{{l"+len(self.exlist)*"r"+"}}")

        lines.append(fr"Label                                      ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"& {self.exlist[iex].name}"
        lines[-1] += fr"\\"

        lines.append(fr"Energy $[\mathrm{{Ry}}]$                   ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].energy:8.3f}       "
        lines[-1] += fr"\\"

        lines.append(fr"Energy $[\mathrm{{eV}}]$                   ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].energy*ev:8.3f}       "
        lines[-1] += fr"\\"

        lines.append(fr"Strength                                   ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].strength:8.3f}        "
        lines[-1] += fr"\\"

        lines.append(fr"Phase $[\pi]$                              ")
        for iex in range(len(self.exlist)):
            pha = (self.exlist[iex].phase/np.pi)%(2.*np.pi)
            if pha >   np.pi: pha-=2.*np.pi
            if pha <= -np.pi: pha+=2.*np.pi
            lines[-1] += fr"&    {pha:8.3f}        "
        lines[-1] += fr"\\"

        lines.append(fr"Area-transition dipoles [a.u.]\\")
        
        for iarea in range(self.narea):
            lines.append(fr"\hline     area {iarea}  -$x$          ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {self.exlist[iex].dipoles[iarea][0]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$y$          ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {self.exlist[iex].dipoles[iarea][1]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$z$          ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {self.exlist[iex].dipoles[iarea][2]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$\|\ldots\|$   ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {np.linalg.norm(self.exlist[iex].dipoles[iarea])*100.:8.3f}        "
            lines[-1] += fr"\\"

        lines.append(fr"Global transition dipoles [a.u.]\\")
        lines.append(fr"\hline                   -$x$          ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].dipole[0]*100.:8.3f}        "
        lines[-1] += fr"\\"
        lines.append(fr"                         -$y$          ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].dipole[1]*100.:8.3f}        "
        lines[-1] += fr"\\"
        lines.append(fr"                         -$z$          ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {self.exlist[iex].dipole[2]*100.:8.3f}        "
        lines[-1] += fr"\\"
        lines.append(fr"                         -$\|\ldots\|$   ")
        for iex in range(len(self.exlist)):
            lines[-1] += fr"&    {np.linalg.norm(self.exlist[iex].dipole)*100.:8.3f}        "
        lines[-1] += fr"\\"

        
        for icalc in range(self.ncalc):
            lines.append(fr"\hline     calc {icalc}  \\")
            lines.append(fr"Area amplitudes [$10^{{-2}}$ a.u.]\\")
            for iarea in range(self.narea):
                lines.append(fr"\hline     area {iarea}  -$x$          ")
                for iex in range(len(self.exlist)):
                    lines[-1] += fr"&    {self.exlist[iex].ampl[icalc][iarea][0]*100.:8.3f}        "
                lines[-1] += fr"\\"
                lines.append(fr"\hphantom{{area {iarea}}}-$y$          ")
                for iex in range(len(self.exlist)):
                    lines[-1] += fr"&    {self.exlist[iex].ampl[icalc][iarea][1]*100.:8.3f}        "
                lines[-1] += fr"\\"
                lines.append(fr"\hphantom{{area {iarea}}}-$z$          ")
                for iex in range(len(self.exlist)):
                    lines[-1] += fr"&    {self.exlist[iex].ampl[icalc][iarea][2]*100.:8.3f}        "
                lines[-1] += fr"\\"
                lines.append(fr"\hphantom{{area {iarea}}}-$\|\ldots\|$   ")
                for iex in range(len(self.exlist)):
                    lines[-1] += fr"&    {np.linalg.norm(self.exlist[iex].ampl[icalc][iarea])*100.:8.3f}        "
                lines[-1] += fr"\\"
            lines.append(fr"Global amplitudes [$10^{{-2}}$ a.u.]\\")
            lines.append(fr"\hline                   -$x$          ")
            ampl = [[sum([self.exlist[iex].ampl[icalc][iarea][icomp] for iarea in range(self.narea)]) for icomp in range(self.ncomp)] for iex in range(len(self.exlist))]
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {ampl[iex][0]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$y$          ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {ampl[iex][1]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$z$          ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {ampl[iex][2]*100.:8.3f}        "
            lines[-1] += fr"\\"
            lines.append(fr"\hphantom{{area {iarea}}}-$\|\ldots\|$   ")
            for iex in range(len(self.exlist)):
                lines[-1] += fr"&    {np.linalg.norm(ampl[iex])*100.:8.3f}        "
            lines[-1] += fr"\\"
        lines.append(fr"\end{{tabular}}")
        lines.append(fr"\caption{{'fit data'}}")
        lines.append(fr"\label{{tab:fit_data}}")
        lines.append(fr"\end{{table}}")

        with open(fname,"w") as fh:
            for line in lines:
                fh.write(line+"\n")

    #-------------------------------------------------------------------------
    # Output amplitude table for gnuplot script
    #-------------------------------------------------------------------------
    def gnuTable(self):
        comp = ["x","y","z"]
        for icalc in range(self.ncalc):
            fname=f"excit_GnuTable_{icalc}.dat"
            lines = []
            for iex, ex in enumerate(self.exlist):
                lines.append(f"w{iex+1} = {ex.energy}")
                lines.append(f"p{iex+1} = {ex.phase}")
                for iarea in range(self.narea):
                    for icomp in range(self.ncomp):
                        lines.append(f"a{iex+1}{comp[icomp]}{iarea+1} = {ex.ampl[icalc][iarea][icomp]}")

            with open(fname,"w") as fh:
                for line in lines:
                    fh.write(line+"\n")
