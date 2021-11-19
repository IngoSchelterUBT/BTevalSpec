import yaml
import numpy as np
import shutil
import os

import errorHandler as err

def writeEmptyConfig():
  data = {  'DESCRIPTION' : '>',
            'DIPOLE' : ['dipole_k100.dat', 'dipole_k010.dat', 'dipole_k001.dat'],
            'EXCIT' : 'laser_profile.dat',
            'OPT' : { 'FourierTransform' :  {  'fourier' : True, 'window' : 0, 'smooth' : 0, 'min_pw' : 17 
                                            },
                      'PadeApprox' : {  'pade' : False, 'pade_wmax' : 0.6, 'pade_dw' : 0.00001, 'pade_smooth' : 0.0000001,
                                        'pade_thin' : 0
                                     },
                      'FitSpectrum' : { 'fit' : False, 'fit_guess' : False, 'fit_range' : [0.1, 0.6], 'guess_thres' : 10},
                    },
            'SPEC' : { 'nex' : 2, 'fit_abserr' : 1.0e-12, 'fit_referr' : 1.0e-4, 'fit_relerr' : 1.0e-8, 
                       'excitations' : [{'name' : 'S1', 'fix' : False, 'energy' : 0.2, 'strength' : 0.1, 'phase' : 0.05,
                                            'transdip' : [0.0, 0.1, 0.2]},
                                        {'name' : 'S2', 'fix' : True, 'energy' : 0.15, 'strength' : 0.12, 'phase' : 0.04,
                                            'transdip' : [0.0, 0.3, 0.4]}]
                     },
          }
  #in excitations is an array of dictionaries
 
  with open('eval_template.yaml','w') as file:
    yaml.dump(data,file)

def readConfig():
  with open('eval.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
   
  return config


def cleanFT():
  if os.path.exists('Osci'):
    shutil.rmtree('Osci')
  if os.path.exists('PW'):
    shutil.rmtree('PW')

def cleanPade():
  if os.path.exists('PADE'):
    shutil.rmtree('PADE')
