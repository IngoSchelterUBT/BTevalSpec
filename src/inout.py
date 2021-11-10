import yaml
import numpy as np

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
                      'FitSpectrum' : { 'fit' : False, 'fit_guess' : False },
                    },
            'SPEC' : { 'excitations' : { 'energies' : [0, 0.2, 0.5], 'strenghts' : [0, 0.2, 0.4]
                                      },
                    },
          }
 
  with open('eval_template.yaml','w') as file:
    yaml.dump(data,file)

def readConfig():
  with open('eval.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
    
  return config

