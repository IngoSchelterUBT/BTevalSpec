#File for calculation the fit

import numpy as np
from scipy.fft import fft,fftfreq
from scipy import signal

import dipole
import padeApprox
import errorHandler as err

#For tests
from pprint import pprint

class Fit:
  
  def __init__(self,config,ft,pade,Id,calcFlag='no'):
    ################################
    #id of fit
    self.fitId = Id
    #read fit range out of config-File
    self.fit_range = config.fit_range
    #read propataion time out of dipole-File
    self.propTime = ft.propTime
    #read guess threshold out of config-File
    if config.fit_guess: self.guess_thres = config.guess_thres

    ################################
    #Filter padeOsci and Osci
    self.osci = []
    self.padeOsci = []
    for i in range(len(ft.osci)):
      self.osci.append(ft.osci[i][(ft.osci[i][:,0] >= self.fit_range[0]) & (ft.osci[i][:,0] <= self.fit_range[1]),:])
      if config.fit_guess:
        self.padeOsci.append(pade.padeOsci[i][(pade.padeOsci[i][:,0] >= self.fit_range[0]) & (pade.padeOsci[i][:,0] <= self.fit_range[1]),:])

    #Make Fit
    self.makeFit(config,calcFlag)


#-----------------------------------------------------------------------------#
#   Methods of the class Fit
#-----------------------------------------------------------------------------#
  #Routine for making fit
  def makeFit(self,config,calcFlag='no'):
    #Make an list for the guess, 
    # - if three files are fitted (calcFlag=='no'), then the list contains the
    #   guess for the x-, y- and z-direction
    # - if only one file is fitted (calcFlag=='trace'), then the list only contains
    #   one guess.
    
    #make guess or read guess out of config
    if config.fit_guess:
      #Make guess out of Pade-Approximation
      if calcFlag == 'no':
        #if x-, y- and z-file are fitted, then calculate a pade-Approximation of the sum of x-,y- and z-file
        #combine the three guesses to one guess, with
        #energy | a1 (x-direction) | a2 (y-direction) | a3 (z-direction)
        self.guess =  self.makeGuess(config,calcFlag)
      elif calcFlag == 'trace':
        #only the trace is fitted, i.e. only one guess has to be made and amplitude is oscillator strength
        self.guess = self.makeGuess(config,calcFlag)
    else:
      #In object config.excitations are the names, energies, osciStrengths, phases and transdips of 
      #all the excitations saved
      self.excitations = config.excitations
      #read guess out of config-file
      if calcFlag == 'no':
        #read a1, a2, a3 as guess for fit
        self.guess = np.column_stack([config.excitations.energies,config.excitations.transdips])
      elif calcFlag == 'trace':
        #read oscillator strength as guess
        self.guess = np.column_stack([config.excitations.energies,config.excitations.osciStrengths])


  #Routine for making guess out of pade-Approximation
  def makeGuess(self,config,calcFlag='no'):
    # 1) Create Pade-Approximation depending on all three files should be fitted or just the trace
    padeSum = self.createPade(config,calcFlag)
    
    # 2) look for peaks in the Pade Approximation (padeSum)
    w = self.searchPeaks(padeSum,calcFlag)

    # 3) Make a guess for the oscillator strength out of the imag part of Osci (imag part of ft of dipole moment)
    #    i.e. search for the corresponding energy and oscillator strength in the fourier transformation
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
    f = f*w[:,np.newaxis]/self.propTime

    return np.column_stack([w,f])

  #Routine for creating padeSum, which is the Pade-Approximation from which the fit is made of
  def createPade(self,config,calcFlag):
    if calcFlag == 'no':
      #create a dipole object which contains the sum of x-, y- and z-component of the fitted dipole moment
      dip = dipole.Dipole(config.dipoleFiles[self.fitId],-1,calcFlag='guess')
      #calculate pade-Approximation out of this summed up dipole File
      pade = padeApprox.Pade(config,dip,-1,calcFlag='guess')
      #set padeSum as this padeApproximation
      padeSum = pade.padeOsci[0][(pade.padeOsci[0][:,0] >= self.fit_range[0]) & (pade.padeOsci[0][:,0] <= self.fit_range[1]),:]
    elif calcFlag == 'trace':
      padeSum = self.padeOsci[0]

    return padeSum

  #Routine for searching Peaks in Pade-Approximation
  def searchPeaks(self,padeSum,calcFlag):
    pos_peaks, _ = signal.find_peaks(padeSum[:,1],height=self.guess_thres)
    #only look for negative peaks, if the trace is NOT fitted (i.e. calcFlag==no)
    if calcFlag == 'no':
      neg_peaks, _ = signal.find_peaks(-padeSum[:,1],height=self.guess_thres)
      w_pade = np.sort(padeSum[np.append(pos_peaks,neg_peaks),0])
    elif calcFlag == 'trace':
      w_pade = np.sort(padeSum[pos_peaks,0])

    return w_pade

  #Routine for whole fit-function
  #x: are the x-values of the fit function
  #guess: is a matrix of first column of energy of peaks and second column of oscillator strength of peaks
  def fit_sinc(x,guess):
    s = np.zeros(len(x))
    for j in range(len(x)):
      s[j] = 0.0
      for i in range(len(guess[:,0])):
        if (x[j]-guess[i,0]) != 0.0:
          s[j] += guess[i,1]/guess[i,0]*np.sin(self.propTime*(x[j]-guess[i,0]))/(x[j]-guess[i,0])
        else:
          s[j] = guess[i,1]/guess[i,0]*self.propTime
    
    return s


