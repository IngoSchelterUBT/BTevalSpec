#===============================================================================
# Defines an excitation
#===============================================================================
class Excitation:
    def __init__(self,exdict=None):
        try:
            self.name     = exdict[name]    
            self.energy   = exdict[energy]  
            self.strength = exdict[strength]
            self.dipole   = exdict[dipole]  
            self.phase    = exdict[phase]   
            self.fix      = exdict[fix]     
        except:
            self.name     = "S"
            self.energy   = 0.
            self.strength = 0.
            self.dipole   = [0.,0.,0.]
            self.phase    = 0.
            self.fix      = False

    def todict(self):
        exdict = {}
        exdict[name]     = self.name
        exdict[energy]   = self.energy
        exdict[strength] = self.strength
        exdict[dipole]   = self.dipole
        exdict[phase]    = self.phase
        exdict[fix]      = self.fix
        return exdict

#=============================================================================
# Defines a list of excitations
#=============================================================================
class Excitations:
    def __init__(self,exdict=None):
        self.exlist = []
        try:
            for ex in exdict["SPEC"]:
                self.add(Excitation(ex))
        except:
            pass

    #=========================================================================
    # Add an excitation to the list
    #=========================================================================
    def add(self,ex):
        self.exlist.append(ex)
        self.sort()

    #=========================================================================
    # Sorts the list of excitations
    #=========================================================================
    def sort(self):
        self.exlist = sorted(self.exlist, key=lambda x: x['energy'])

    #=========================================================================
    # Converts excitation list to dictionary
    #=========================================================================
    def excit2dict(self):
        exdict = {"SPEC": []}
        for ex in exlist:
            exdict["SPEC"].append(ex.todict())
        return exdict
