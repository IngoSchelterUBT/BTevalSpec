#==============================================================================#
# Configuration
#==============================================================================#
#intrinsic
import numpy as np
import ruamel.yaml
import sys

#custom
import errorHandler as err

#==============================================================================#
# Class Config
#==============================================================================#
class Config:
    def __init__(self,ifile):
        yaml = ruamel.yaml.YAML()
        with open(ifile) as fh:
            self.conf = yaml.load(fh)
        self.descript = self.conf.get("DESCRIPTION","")
        self.dipfiles = self.conf.get("DIPOLE"     ,[])
        self.ext      = self.conf.get("EXT"        ,{})
        self.densft   = self.conf.get("DENSFT"     ,{})
        self.opt      = self.conf.get("OPT"        ,{})
        self.excit    = self.conf.get("SPEC"       ,[])
        tmp, self.ind, self.bsi = ruamel.yaml.util.load_yaml_guess_indent(open(ifile))

    #--------------------------------------------------------------------------#
    # Write configuration
    #--------------------------------------------------------------------------#
    def write(self,ofile):
        yaml = ruamel.yaml.YAML()
        self.conf["DESCRIPTION"] = self.descript
        self.conf["DIPOLE"     ] = self.dipfiles
        self.conf["EXT"        ] = self.ext
        self.conf["DENSFT"     ] = self.densft
        self.conf["OPT"        ] = self.opt
        self.conf["SPEC"       ] = self.excit
        yaml.indent(mapping=self.ind,sequence=self.ind,offset=self.bsi)
        with open(ofile,'w') as fh:
            yaml.dump(self.conf,fh)

#==========================================================================#
# Write configuration template
#==========================================================================#
def writeTemplate(ofile):
    conf = """\
    DESCRIPTION:                    # Description of evaluation 
    DIPOLE:                         # List of dipole moment files; different calculations are only allowed to differ in the boost direction or laser polarization
      - [dipole_calc1_area1.dat, dipole_calc1_area2.dat]
      - [dipole_calc2_area1.dat, dipole_calc2_area2.dat]
    EXT:
      profile: laser_profile.dat    # Profile of excitation
      invertPhase: False            # Before BTDFT v3.5.2, the laser profile missed a factor "-1". InvertPhase==True compensates this error.
    DENSFT:
      densft: []                    #if cube files: must contain a 2-list with real/imag part cube files per energy / if compact files: must contain one complex-valued compact file per energy
      densen: []
    OPT:
      FT:
        calc: true
      Pade:
        calc: true
      Fit:
        range:                      # Range of spectrum which should be fitted in Ry
          - 0.10
          - 0.40
    """
    yaml = ruamel.yaml.YAML()
    code = yaml.load(conf)
    with open(ofile,'w') as fh:
        yaml.dump(code,fh)
