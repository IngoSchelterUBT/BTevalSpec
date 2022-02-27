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
    if conf.gnuplot_spectrum: self.printGnuplotSpectrum()
    if conf.dat_spectrum: self.writeDatFile()



#-----------------------------------------------------------------------------#
#   Methods of the class Fit
#-----------------------------------------------------------------------------#

  #Routine for making list of dictionaries for excitations
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

  #Routine for writing a gnuplot file for plotting the spectrum
  def printGnuplotSpectrum(self):
    ev = 13.605693122994
    w = np.array([d['energy'] for d in self.excs['SPEC']['excitations']])
    f = np.array([d['strength'] for d in self.excs['SPEC']['excitations']])
    name = [d['name'] for d in self.excs['SPEC']['excitations']]

    sp = open('spectrum.plt','w')
    sp.write('reset')
    sp.write('\nset sample 1000')
    sp.write('\nset xrange [' + str(np.floor(np.amin(w)*ev)) + ':' + str(np.ceil(np.amax(w)*ev)) + ']')
    sp.write("\nset xlabel 'Energy (eV)'")
    sp.write("\nset ylabel 'Oscillator strength'")
    sp.write('\n\n###############################################################################')
    sp.write('\n\nev = 13.605693122994')
    sp.write('\nn = 0.025')
    sp.write('\n\n#Excitations')
    for i,val in enumerate(w):
      sp.write('\nw_' + str(i+1) + ' = ' + str(w[i]))
      sp.write('\nf_' + str(i+1) + ' = ' + str(f[i]))
      
    sp.write('\n\n#Single lines')
    for i,val in enumerate(w):
      sp.write('\ns_' + str(i+1) + '(x) = f_' + str(i+1) + '*exp(-((x/ev-w_' + str(i+1) + ')/(n/ev))**2)')
    sp.write('\n\n#Spectrum')
    sp.write('\nsp(x) = \\')
    sp.write('\n        ')
    for i,val in enumerate(w):
      if i != len(w)-1:
        sp.write('s_' + str(i+1) + '(x) + ')
      else:
        sp.write('s_' + str(i+1) + '(x)') 
    
    sp.write('\n\nunset arrow')


    sp.write('\n\n#Arrows')
    for i,val in enumerate(w):
      sp.write('\nset arrow from w_' + str(i+1) + '*ev,0 to w_' + str(i+1) + '*ev,f_' + str(i+1) + ' lw 2 head filled')
      sp.write(f"\nset label '{name[i]}' at w_" + str(i+1) + '*ev,f_' + str(i+1) + ' offset 0,0.5')

    sp.write('\n\n#Plot')
    sp.write('\nplot \\')
    sp.write("\n      sp(x) lc rgb 'black' dt 1 lw 2 ti 'Spectrum',\\")
    for i,val in enumerate(w):
      if i == 0:
        sp.write("\n      s_" + str(i+1) + "(x) lc 4 dt 2 lw 2 ti 'Excitations',\\")
      else:
        sp.write("\n      s_" + str(i+1) + "(x) lc 4 dt 2 lw 2 noti,\\")
  
    sp.write("\n      0 lc rgb 'black' lw 2 noti\n")

    sp.write('\n\nset output')
    
    sp.close()

  #Routine for writing spectrum data in dat file
  def writeDatFile(self):
    ev = 13.605693122994
    w = np.array([d['energy'] for d in self.excs['SPEC']['excitations']])
    f = np.array([d['strength'] for d in self.excs['SPEC']['excitations']])

    np.savetxt('spectrum.dat', np.column_stack([w,f]),header='Energy (eV) | Oscillator strength')
