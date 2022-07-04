#File for calculation and plot of fit

import numpy as np

#own modules
import errorHandler as err
import util

#TODO: maybe create a config-keyword for the realtive error threshold (0.01 at this moment)

#Routine for handling trace guess: The guess for trace consists only of the
#exciations lines out of the fit for the different dipole moments (x-, y- and
#z-direction)
#The routine also sets Ids for the excitation lines in the directions (for
#the trace the Id is straight forward, i.e. 0, 1, 2, 3, ...)
#input variable:
# - guessTrace object of class Trace (fourier transform of dipole moment 
#   in the fit range is saved there)
# - list of fit objects for fits of dipole moments with boost in x-, y- and
#   z-direction
def guessTrace(config, guessTrace, fit):
  #calculate relative error between the excitaion energies of the fits with boost
  #in x-, y- and z-direction
  w_trace = fit[0].fit_result[:,0]
  fix = fit[0].fix
  for i in range(1,len(fit)):
    for j in range(len(fit[i].fit_result[:,0])):
      err = np.abs(w_trace - fit[i].fit_result[j,0])/((w_trace + fit[i].fit_result[j,0])/2.)
      if np.any(err < config.fit_relspacing_lines): #default relspacing_lines: 0.01 equals a relative error of 1.0 %
        continue
      else:
        #print('line added:')
        #print(fit[i].fit_result[j,0])
        #print(i)
        w_trace = np.append(w_trace, fit[i].fit_result[j,0])
        fix = np.append(fix, fit[i].fix[j])

  arginds = w_trace.argsort()
  w_trace = w_trace[arginds]
  fix = fix[arginds]

  for i in range(len(fit)):
    fit[i].lineId = np.full(len(fit[i].fit_result[:,0]),-1,dtype=np.int32)
    for j, value in enumerate(w_trace):
      err = np.abs(value - fit[i].fit_result[:,0])/((value + fit[i].fit_result[:,0])/2.) 
      indices = np.argwhere(err < config.fit_relspacing_lines) #default relspacing_lines: 0.01 equals a relative error of 1.0 %
      fit[i].lineId[indices] = j

      
      

  #make a guess for the oscillator strengths out of the FT of the trace
  f_trace = guessTrace.searchAmplitude(w_trace)

  #save in guess for trace
  guessTrace.guess = np.column_stack([w_trace,f_trace])
  guessTrace.fix = fix
    
  
