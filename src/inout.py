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
                      'FitSpectrum' : { 'fit' : False },
                    },
            'SPEC' : { 'excitations' : { 'energies' : [0, 0.2, 0.5], 'strenghts' : [0, 0.2, 0.4]
                                      },
                    },
          }
 
  with open('eval.yaml','w') as file:
    yaml.dump(data,file)

def readConfig():
  with open('eval.yaml') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
    
  return config

def readDipole(dipolefile):
  dipole = {}
  with open(dipolefile,'r') as f:
    for line in f:
      if '#!BT' in line:
        line_array = line.split()
        if len(line_array) == 3:
          (key, val) = (line_array[1],line_array[2])
          dipole[key] = [val]
        if len(line_array) == 4:
          (key, [val1, val2]) = (line_array[1],[line_array[3],line_array[2]])
          dipole[key] = [val1, val2]
  
  #Check units of Dipole File
  if dipole.get('DX')[-1] != 'a0':
    err.err(1,'Unit of DX unknown!')
  elif dipole.get('T_START')[-1] != 't_ry':
    err.err(1,'Unit of T_START unknown!')
  elif dipole.get('DT')[-1] != 't_ry':
    err.err(1,'Unit of DT unknown!')
  elif dipole.get('BOOSTENERGY')[-1] != 'Ry':
    err.err(1,'Unit of BOOSTENERGY unknown!')
  elif dipole.get('LASERFREQ')[-1] != 'Ry':
    err.err(1,'Unit of LASERFREQ unknown!')
  elif dipole.get('LASERINT')[-1] != 'W/cm^2':
    err.err(1,'Unit of LASERINT unknown!')
  elif dipole.get('LASEREFIELD')[-1] != 'E_ry':
    err.err(1,'Unit of LASEREFIELD unknown!')
  elif dipole.get('LASERSTART')[-1] != 't_ry':
    err.err(1,'Unit of LASERSTART unknown!')
  elif dipole.get('LASEREND')[-1] != 't_ry':
    err.err(1,'Unit of LASEREND unknown!')
  elif dipole.get('DIP0')[-1] != 'dip_ry':
    err.err(1,'unit of DIP0 unknown!')
  
  #To Do: Kommasepariertes Array als String umformatieren als numpy array!
  for k, v in dipole.items():
    if ',' in v[0]:
      dipole[k] = np.array([float(i) for i in v[0].split(',')]) #modifiziere v[0]: dort steht der value drinnen, d. h. wenn etwas Komma seperiertes ist dann aendere den value zu einer numpy list
    if isfloat(v[0]):
      v[0] = float(v[0])

  return dipole

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
