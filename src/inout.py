#import modules
import numpy as np
import shutil
import os
import ruamel.yaml
import re

#own modules
import errorHandler as err

def writeEmptyConfig():
  yaml = ruamel.yaml.YAML()
  
  data = {  'DESCRIPTION': '>',
            'DIPOLE': ['dipole_k100.dat', 'dipole_k010.dat', 'dipole_k001.dat'],
            'EXCIT': 'laser_profile.dat',
            'OPT': { 'FourierTransform':  {  'fourier': True, 'window': 0, 'smooth': 0, 'min_pw': 17
                                            },
                      'PadeApprox' : {  'pade': False, 'pade_wmax': 0.6, 'pade_dw': 0.00001, 'pade_smooth': 0.0000001,
                                        'pade_thin': 0
                                     },
                      'FitSpectrum': { 'fit': False, 'fit_guess': False, 'fit_relerr_crit': 0.1, 'fit_max_iter' : 10,'fit_range': [0.1, 0.6], 'guess_thres': 0.1, 'fit_relspacing_lines': 0.01},
                    },
            'SPEC': { 'nex': 2,
                       'excitations': [{'name': 'S1', 'fix': False, 'energy': 0.2, 'strength': 0.1, 'phase': 0.05,
                                            'transdip': [0.0, 0.1, 0.2]},
                                        {'name': 'S2', 'fix': True, 'energy': 0.15, 'strength': 0.12, 'phase': 0.04,
                                            'transdip': [0.0, 0.3, 0.4]}]
                     },
          }
  #in excitations is an array of dictionaries

  data_str = """\
  DESCRIPTION:                    # Description of evaluation 
  DIPOLE:                         # List of dipole moment files
    - dipole_k100.dat
    - dipole_k010.dat
    - dipole_k001.dat
  EXCIT: laser_profile.dat        # profile of excitation
  OPT:
    FourierTransform:
      fourier: True               # Turn fourier transformation on/off
      min_pw: 17                  # Min. length of fourier transform vector (2^n). Can increase sampling rate.
      smooth: 0                   # Artificial decay rate
      window: 0                   # Kaiser-Bessel windowing parameter (>0)
    PadeApprox:
      pade: True                  # Turn Pade Approximation on/off
      pade_wmax: 0.4              # Maximum energy for Pade Approximation in Ry
      pade_dw: 1.0e-05            # Step of Pade Approximation
      pade_smooth: 1.0e-07        # Smooth for Pade Approximation
      pade_thin: 0                # Only keep every 2^n data point
    FitSpectrum:
      fit: True                   # Turn fit on/off
      fit_guess: True             # Turn guess for fit via Pade Approximation on/off
      guess_thres: 0.1            # Relative height of line in Pade Approximation compared to highest line which should be identified as a line for fitting (only relevent, if fit_guess == True).
      fit_relerr_crit: 0.05       # Criterium for relative error between fit and data (only relevant, if fit_guess == True)
      fit_max_iter: 10            # Maximum numer of iterations used to reach relative error between fit and data (only relevant, if fit_guess == True)
      fit_relspacing_lines: 0.01  # threshold for relative error between two lines which should be identified as one in fit of trace (only relevant, if number of dipole files == 3).
      fit_range:                  # Range of spectrum which should be fitted in Ry
        - 0.10
        - 0.40
  """
  code = yaml.load(data_str)
  with open('eval.yaml','w') as file:
    yaml.dump(code,file)

def readConfig():
  yaml = ruamel.yaml.YAML()
  with open('eval.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

  return config

#Routine for writing excitation information
#fit is a list of all fits (x-, y- and z-direction and trace)
#TODO: Change boolean guess=True to guess=False
def writeExcitations(conf,fit): #attention conf is the object!!
  #read yaml file again
  yaml = ruamel.yaml.YAML()
  config, ind, bsi = ruamel.yaml.util.load_yaml_guess_indent(open('eval.yaml'))

  #create a list of dictionaries for the excitation lines
  try:
    config['SPEC']
  except Exception:
    config['SPEC'] = {}
  excitations = []
  if len(fit) == 4:
    config['SPEC']['nex'] = len(fit[3].fit_result[:,0])
    for i in range(len(fit[3].fit_result[:,0])):
      #create empty dictionary
      exc = {}
      try:
        exc['name'] = str(conf.excitations.names[i])
      except Exception:
        exc['name'] = 'S'
      try:
        exc['fix'] = boole(fit[3].fix[i])
      except Exception:
        exc['fix'] = False
      exc['energy'] = float(fit[3].fit_result[i,0])
      exc['strength'] = float(fit[3].fit_result[i,1])
      for j in range(3):
        index = np.where(fit[j].lineId == i)[0]
        if len(index) != 0:
          #ruamel.yaml only takes lists NOT numpy arrays
          exc['amplitude_file%i' % (j+1)] = np.ndarray.tolist(fit[j].fit_result[index,1:])[0]
        else:
          exc['amplitude_file%i' % (j+1)] = [0., 0., 0.]
      excitations.append(exc)
  elif len(fit) == 1:
    config['SPEC']['nex'] = len(fit[0].fit_result[:,0])
    for i in range(len(fit[0].fit_result[:,0])):
      exc = {}
      try:
        exc['name'] = str(conf.excitations.names[i])
      except Exception:
        exc['name'] = 'S'
      try:
        exc['fix'] = bool(fit[0].fix[i])
      except Exception:
        exc['fix'] = False
      exc['energy'] = float(fit[0].fit_result[i,0])
      exc['strength'] = float(fit[0].osciStrength[i])
      exc['amplitudes_file1'] = np.ndarray.tolist(fit[0].fit_result[i,1:])
      excitations.append(exc)

  excitations = sorted(excitations, key=lambda x: x['energy'])

  config['SPEC']['excitations'] = excitations

  yaml.indent(mapping=ind, sequence=ind, offset=bsi)
  with open('eval.yaml', 'w') as file:
    yaml.dump(config, file)

def cleanFT():
  if os.path.exists('Osci'):
    shutil.rmtree('Osci')
  if os.path.exists('PW'):
    shutil.rmtree('PW')

def cleanPade():
  if os.path.exists('PADE'):
    shutil.rmtree('PADE')
