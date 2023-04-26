#==============================================================================#
# Load modules
#==============================================================================#
#intrinsic
import numpy as np
import ruamel.yaml

#custom
import inout
import errorHandler as err
import excitations

#==============================================================================#
# Class Config
#==============================================================================#
class Config:
    def __init__(self,ifile):
        yaml = ruamel.yaml.YAML()
        with open(ifile) as fh:
            conf = yaml.load(fh)
        self.descript = conf.get("DESCRIPTION","")
        self.dipfiles = conf.get("DIPOLE"     ,[])
        self.extfile  = conf.get("EXCIT"      ,"")
        self.opt      = conf.get("OPT"        ,{})
        self.excit    = conf.get("SPEC"       ,[])

    #==========================================================================#
    # Write configuration
    #==========================================================================#
    def write(self,ofile):
        yaml = ruamel.yaml.YAML()
        conf = {"DESCRIPTION": self.descript, "DIPOLE": self.dipfiles, "EXCIT": self.extfiles, "OPT": self.opt, "SPEC": self.excit}
        with open(ofile,'w') as fh:
            yaml.dump(conf,fh)

#==========================================================================#
# Write configuration template
#==========================================================================#
def writeTemplate(ofile):
    conf = """\
    DESCRIPTION:                    # Description of evaluation 
    DIPOLE:                         # List of dipole moment files
      - [dipole_calc1_area1.dat, dipole_calc1_area2]
      - [dipole_calc2_area1.dat, dipole_calc2_area2]
    EXCIT: laser_profile.dat        # Profile of excitation
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
        smooth: 1.0e-07             # Smooth for Pade Approximation
        thin: 0                     # Only keep every 2^n data point
      Fit:
        calc: True                  # Turn fit on/off
        guess: True                 # Turn guess for fit via Pade Approximation on/off
        plot_result: False          # If True: The fit results are plotted without fitting again
        gnuplot_spectrum: False     # If True: A gnuplot script for plotting the resulting spectrum is created.
        dat_spectrum: False         # If True: A .dat file for the spectrum is created
        guess_thres: 0.1            # Relative height of line in Pade Approximation compared to highest line 
                                    #   which should be identified as a line for fitting (only relevent, if 
                                    #   fit_guess == True).
        relerr_crit: 0.05           # Criterium for relative error between fit and data (only relevant, if 
                                    #   fit_guess == True).
        max_iter: 10                # Maximum numer of iterations used to reach relative error between fit and 
                                    #   data (only relevant, if fit_guess == True)
        relspacing_lines: 0.01      # Threshold for relative error between two lines which should be identified 
                                    #   as one in fit of trace (only relevant, if number of dipole files == 3).
        range:                      # Range of spectrum which should be fitted in Ry
          - 0.10
          - 0.40
    SPEC:
      - name:     "Dummy"
        energy:   0.
        strength: 0.
        dipole:   [0.,0.,0.]
        phase:    0.
        fix:      False
    """
    yaml = ruamel.yaml.YAML()
    code = yaml.load(conf)
    with open(ofile,'w') as fh:
        yaml.dump(code,fh)
