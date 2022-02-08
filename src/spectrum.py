#Import Libraries
import numpy as np
import ruamel.yaml

#==============================================================================#
# Class Spectrum
#==============================================================================#
# The Class consists of the following variables:
# - 
#------------------------------------------------------------------------------#

class Spectrum:
  def __init__(self,conf,fit):
    self.excs = self.createExcDict(conf,fit)



#-----------------------------------------------------------------------------#
#   Methods of the class Fit
#-----------------------------------------------------------------------------#
  def createExcDict(self,conf,fit):
    excitations = []
    if len(fit) == 4:
      nex = len(fit[3].fit_result[:,0])
      for i in range(len(fit[3].fit_result[:,0])):
        #create empty dictionary
        exc = {}
        try:
          exc['name'] = str(conf.excitations.names[i])
        except Exception:
          exc['name'] = 'S'
        try:
          exc['fix'] = bool(fit[3].fix[i])
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
      nex = len(fit[0].fit_result[:,0])
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

    excDict = {}
    excDict['SPEC'] = {}
    excDict['SPEC']['nex'] = nex
    excitations = sorted(excitations, key=lambda x: x['energy'])
    excDict['SPEC']['excitations'] = excitations


    return excDict
