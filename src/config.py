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
      - [dipole_calc1_area1.dat, dipole_calc1_area2]
      - [dipole_calc2_area1.dat, dipole_calc2_area2]
    EXT:
      profile: laser_profile.dat  # Profile of excitation
    OPT:
      FT:
        calc: True                  # Turn fourier transformation on/off
        minpw: 17                   # Min. length of fourier transform vector (2^n). Can increase sampling rate.
        smooth: 0                   # Artificial decay rate
        window: 0                   # Kaiser-Bessel windowing parameter (>0)
      Pade:
        calc: True                  # Turn Pade Approximation on/off
        wmax: 0.6                   # Maximum energy for Pade Approximation in Ry
        dw: 1.0e-05                 # Step of Pade Approximation
        smooth: 0.                  # Smooth for Pade Approximation (0:auto)
        thin: 0                     # Only keep every 2^n data point
      Fit:
        calc: True                  # Turn fit on/off
        skipfirst: False            # Skip first fit of existing excitations
        guess: True                 # Turn guess for fit via Pade Approximation on/off
        plot_result: True           # If True: The fit results are plotted without fitting again
        gnuplot_spectrum: False     # If True: A gnuplot script for plotting the resulting spectrum is created.
        dat_spectrum: False         # If True: A .dat file for the spectrum is created
        guess_thres: 0.05           # Relative height of line in Pade Approximation compared to highest line 
                                    #   which should be identified as a line for fitting (only relevent, if 
                                    #   fit_guess == True).
        relerr_crit: 0.05           # Criterium for relative error between fit and data (only relevant, if 
                                    #   fit_guess == True).
        max_excit: 10               # Maximum numer of iterations used to reach relative error between fit and 
                                    #   data (only relevant, if fit_guess == True)
        relspacing_lines: 0.01      # Threshold for relative error between two lines which should be identified 
                                    #   as one in fit of trace (only relevant, if number of dipole files == 3).
        range:                      # Range of spectrum which should be fitted in Ry
          - 0.10
          - 0.40
        nsigma: 2.                  # New excitations are guessed at energies where the peak height is larger than
                                    # the mean value plus nsigma times the standard deviation of the distribution of
                                    # maxima of Fourier Transform minus current fit
                                    # (and the latter scaled for error compensation around existing lines)
        significances: False
        fiterr: 0.                  # Current fit error
    SPEC:
      - name:         "S1"                    #Identifier
        fix:          False                   #Set to true if the excitation shall be unchanged in a restart run
        energy:       0.                      #Energy [Ry]            can be changed manually to help fit
        energyErr:    0.                      #Energy Error (directly from fit)
        erange:       [0.1,0.5]               #Energy fit range
        phase:        0.                      #Phase (==0. for boost) can be changed manually to help fit
        phaseErr:     0.                      #Phas Error (directly from fit)
        dipoles:      [[0.,0.,0.],[0.,0.,0.]] #Areas' transition-dipole contributions
        dipolesErr:   [[0.,0.,0.],[0.,0.,0.]] #Dipoles Error (directly from fit)
        dipole:       [0.,0.,0.]              #Total Transition dipole
        dipoleErr:    [0.,0.,0.]              #Dipole Error (derived from dipolesErr)
        strength:     0.                      #Oscillator strength
        strengthErr:  0.                      #Strength Error (derived from dipolesErr and energyErr)
        signifFit:    0.                      #Significanc: Fit
        signifAng:    0.                      #Significanc: sqrt(cos(angle between ext polarization and dipole))
        signifExc:    0.                      #Significanc: Fitting a single line again
        signifErr:    0.                      #Significanc: From Fitting-Error
    """
    yaml = ruamel.yaml.YAML()
    code = yaml.load(conf)
    with open(ofile,'w') as fh:
        yaml.dump(code,fh)
