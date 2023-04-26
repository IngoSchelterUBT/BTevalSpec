#Import Libraries
import numpy as np
from scipy import signal
import matplotlib
import matplotlib.pyplot as plt

#own modules
import dipole
import padeApprox
import errorHandler as err

#TODO: Delete making guess for trace! In this case only create self.guess = [] this is
#      handeled after the fit of the 3 directions.

#==============================================================================#
# Class Fit
#==============================================================================#
# The Class consists of the following variables:
# - guessId: Id of the Guess (starting with 0).
# - fit_range: Range which should be fitted.
# - propTime: Propagation time of the BTDFT calculation.
# - osci: Fourier transformation of dipole file (only for fit range).
# - padeOsci: Pade Approximation of dipole file (only for fit range).
# - guessThres: Relative threshold, when peak should be identified as 
#               excitation line
# - guess: Guess for fitting spectrum (energy | x-amp | y-amp | z-amp) for
#          trace (energy | oscistrength)
#------------------------------------------------------------------------------#

class Guess:
  def __init__(self,config,ft,pade,Id,calcFlag='no'):
    ################################
    #id of guess
    self.guessId = Id
    #read fit range oout of config-object
    self.fit_range = config.fit_range
    #read propagation time of ft-object
    self.propTime = ft.propTime
    #read relative guess threshold out of config-object
    if config.fit_guess: self.guess_thres = config.guess_thres

    ################################
    #Filter padeOsci and Osci
    self.osci = []
    self.padeOsci = []
    for i in range(len(ft.osci)):
      self.osci.append(ft.osci[i][(ft.osci[i][:,0] >= self.fit_range[0]) & (ft.osci[i][:,0] <= self.fit_range[1]),:])
      if config.fit_guess:
        self.padeOsci.append(pade.padeOsci[i][(pade.padeOsci[i][:,0] >= self.fit_range[0]) & (pade.padeOsci[i][:,0] <= self.fit_range[1]),:])

    #run acutal guess
    self.giveGuess(config,calcFlag) #saves guess in self.guess


#-----------------------------------------------------------------------------#
#   Methods of the class Guess
#-----------------------------------------------------------------------------#
  #Head-Routine for making Guess
  def giveGuess(self,config,calcFlag='no'):
    #Make a list for the guess,
    # - if three files are fitted (calcFlag=='no'), then the list contains the
    #   guess for the x-, y- and z-direction
    # - if only one file is fitted (calcFlag=='trace'), then the list only contains
    #   one guess.

    #Make guess or read guess out of config
    if config.fit_guess:
      self.guess = self.makeGuess(config,calcFlag)
    else:
      #In object config.excitations are the names, energies, osciStrengths, phases and transdips of
      #all the excitations saved
      self.guess = np.zeros((len(config.excitations.energies),1+len(self.osci)))
      self.guess[:,0] = config.excitations.energies
      self.osciStrengths = config.excitations.osciStrengths
      #read guess out of config-file
      if calcFlag == 'no':
        #- read a1 oder a2 oder a3 as guess for fit
        #- for new added line search in Osci for guess of amplitude
        for i in range(len(config.excitations.energies)):
          if True in np.isnan(config.excitations.amplitudes[self.guessId][i,:]):
            #new added line -> search in Osci for guess of amplitude
            #in searchAmplitude(w) w must be an array of energies, in this case of
            #length 1
            amp = self.searchAmplitude(np.array([config.excitations.energies[i]]))
            self.guess[i,1:] = amp
          else:
            self.guess[i,1:] = config.excitations.amplitudes[self.guessId][i,:]
      elif calcFlag == 'trace':
        #read oscillator strength as guess
        #- read strength
        #- for new added line search in Osci for guess of strength
        for i in range(len(config.excitations.energies)):
          if config.excitations.osciStrengths[i] == None:
            strength = self.searchAmplitude(np.array([config.excitations.energies[i]]))
            self.guess[i,1:] = strength
          else:
            self.guess[i,1:] = config.excitations.osciStrengths[i]


  #Routine for making guess out of pade-Approximation
  def makeGuess(self,config,calcFlag='no'):
    # 1) Create Pade-Approximation depending on all three files should be fitted or just the trace
    padeSum = self.createPade(config,calcFlag)

    # 2) look for peaks in the Pade Approximation (padeSum)
    w = self.searchPeaks(padeSum,calcFlag)

    # 3) Make a guess for the oscillator strength out of the imag part of Osci (imag part of ft of dipole moment)
    #    i.e. search for the corresponding energy and oscillator strength in the fourier transformation
    f = self.searchAmplitude(w)

    return np.column_stack([w,f])

  #Routine for creating padeSum, which is the Pade-Approximation from which the fit is made of
  def createPade(self,config,calcFlag):
    if calcFlag == 'no':
      #if x-, y- and z-file are fitted, then calculate a pade-Approximation of the sum of x-,y- and z-file
      #combine the three guesses to one guess, with
      #energy | a1 (x-direction) | a2 (y-direction) | a3 (z-direction)
      #create a dipole object which contains the sum of x-, y- and z-component of the fitted dipole moment
      dip = dipole.Dipole(config.dipoleFiles[self.guessId],-1,calcFlag='guess')
      #calculate pade-Approximation out of this summed up dipole File
      pade = padeApprox.Pade(config,dip,-1,calcFlag='guess')
      #set padeSum as this padeApproximation
      padeSum = pade.padeOsci[0][(pade.padeOsci[0][:,0] >= self.fit_range[0]) & (pade.padeOsci[0][:,0] <= self.fit_range[1]),:]
    elif calcFlag == 'trace':
      #only the trace is fitted, i.e. only one guess has to be made and amplitude is oscillator strength
      padeSum = self.padeOsci[0]

    return padeSum

  #Routine for searching Peaks in Pade-Approximation
  def searchPeaks(self,padeSum,calcFlag):
    maxPadeOsci = self.findMaxPeak(padeSum,calcFlag)
    pos_peaks, _ = signal.find_peaks(padeSum[:,1],height=self.guess_thres*maxPadeOsci)
    #only look for negative peaks, if the trace is NOT fitted (i.e. calcFlag==no)
    if calcFlag == 'no':
      neg_peaks, _ = signal.find_peaks(-padeSum[:,1],height=self.guess_thres*maxPadeOsci)
      w_pade = np.sort(padeSum[np.append(pos_peaks,neg_peaks),0])
    elif calcFlag == 'trace':
      w_pade = np.sort(padeSum[pos_peaks,0])

    return w_pade

  #Routine for searching Amplitudes in Osci, i.e. make Guess for amplitudes
  def searchAmplitude(self,w):
    f_temp = [np.array([])]*len(self.osci)
    for i in range(len(self.osci)):
      for j in range(len(w)):
        abs_diff = np.abs(self.osci[i][:,0]-w[j])
        index_smallest_diff = abs_diff.argmin()
        #f as oscillator strength of imaginary part (third column of osci)
        f_temp[i] = np.append(f_temp[i],self.osci[i][index_smallest_diff,2])

    f = np.empty((len(w),0))
    for i in range(len(self.osci)):
      f = np.column_stack([f,f_temp[i]])

    #calculate values for guess of fit-function
    for i in range(len(f[0,:])):
      f[:,i] = f[:,i]*w[:]/self.propTime

    return f


  #Routine for finding maximum peak in PadeOsci
  def findMaxPeak(self,padeSum,calcFlag):
    pos_peaks, _ = signal.find_peaks(padeSum[:,1])
    #only look for negative peaks,if the trace is NOT fitted (i.e. calcFlag==no)
    if calcFlag == 'no':
      neg_peaks, _ = signal.find_peaks(-padeSum[:,1])
      f_pade = np.sort(padeSum[np.append(pos_peaks,neg_peaks),1])
    elif calcFlag == 'trace':
      f_pade = np.sort(padeSum[pos_peaks,1])

    return np.amax(f_pade)
