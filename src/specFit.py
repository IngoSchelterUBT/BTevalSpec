#File for calculation the fit

import numpy as np
from scipy import signal

import dipole
import padeApprox
import errorHandler as err
from lmfit import Parameters, minimize, report_fit

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
    
    # 1.) Make guess or read guess out of config
    if config.fit_guess:
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

    # 2.) Make multifit or make single fit for trace
    if calcFlag == 'no':
      self.makeMultifit()
      #self.plotMultifit()
      return
    elif calcFlag == 'trace':
      #self.makeTracefit()
      #self.plotTracefit()
      return

#-----------------------------------------------------------------------------#
#   Methods for Making Guess
#-----------------------------------------------------------------------------#
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
    for i in range(len(f[0,:])):
      f[:,i] = f[:,i]*w[:]/self.propTime

    return np.column_stack([w,f])

  #Routine for creating padeSum, which is the Pade-Approximation from which the fit is made of
  def createPade(self,config,calcFlag):
    if calcFlag == 'no':
      #if x-, y- and z-file are fitted, then calculate a pade-Approximation of the sum of x-,y- and z-file
      #combine the three guesses to one guess, with
      #energy | a1 (x-direction) | a2 (y-direction) | a3 (z-direction)
      #create a dipole object which contains the sum of x-, y- and z-component of the fitted dipole moment
      dip = dipole.Dipole(config.dipoleFiles[self.fitId],-1,calcFlag='guess')
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
    pos_peaks, _ = signal.find_peaks(padeSum[:,1],height=self.guess_thres)
    #only look for negative peaks, if the trace is NOT fitted (i.e. calcFlag==no)
    if calcFlag == 'no':
      neg_peaks, _ = signal.find_peaks(-padeSum[:,1],height=self.guess_thres)
      w_pade = np.sort(padeSum[np.append(pos_peaks,neg_peaks),0])
    elif calcFlag == 'trace':
      w_pade = np.sort(padeSum[pos_peaks,0])

    return w_pade

#-----------------------------------------------------------------------------#
#   Methods for Making Multifit
#-----------------------------------------------------------------------------#
  #Routine for multifit
  def makeMultifit(self):
    #Make a list of w_i and f_i for all values of x-, y- and z-
    fit_params = Parameters()
    self.osci_range = np.abs(self.osci[0][0,0] - self.osci[0][len(self.osci[0][:,0])-1,0])
    for j in range(len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        fit_params.add('w_%i' % (i+1+off), value=self.guess[i,0])
        fit_params.add('f_%i' % (i+1+off), value=self.guess[i,j+1])
      
    #Make constraint that 'w_1' must equal 'w_(1+off)'
    for j in range(1,len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        fit_params['w_%i' % (i+1+off)].expr = 'w_%i' % (i+1)

    #calculate Fit
    #x is a list of numpy arrays with the energy axis of every direction
    x = []
    for j in range(len(self.guess[0,1:])):
      x.append(self.osci[j][:,0])
    #data is a list of the imaginary part of the fourier transform
    data = self.osci[0][:,2]
    for j in range(1,len(self.guess[0,1:])):
      data = np.append(data,self.osci[j][:,2])

    #Ausgabe zum Testen!
    #for j in range(len(self.guess[0,1:])):
    #  np.savetxt('guess_%i.dat' % (j),np.column_stack([self.guess[:,0],self.guess[:,j+1]]))
    #for j in range(len(x)):
    #  np.savetxt('test_%i.dat' % (j),np.column_stack([x[j],self.osci[j][:,2]]))
    #y = self.fit_sinc(fit_params,x)
    #omega = x[0]
    #for j in range(1,len(x)):
    #  omega = np.append(omega,j*self.osci_range+x[j])
    #np.savetxt('guess_test.dat',np.column_stack([omega,y]))
    #exit()

    #fuehre den Fit durch
    out = minimize(self.objective, fit_params, args=(x,data))
    report_fit(out.params)


  #Routine for whole fit-function
  #x is a list of numpy arrays with the energy axis of every direction
  def fit_sinc(self, w, f, x):
    sinc = np.array([])
    for k in range(len(x)):
      s = np.zeros(len(x[k]))
      off = int(k*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        s += f[i+off]/w[i+off]*\
             np.sinc(self.propTime*(x[k]-w[i+off]))*self.propTime
      
      sinc = np.append(sinc,s)
    
    return sinc


  #Routine for calculatin the total residual for fits of sinc-Function to osci_imag
  #x is an array appending self.osci[0][:,0], self.osci[1][:,0], self.osci[2][:,0]
  #data is an array appending self.osci[0][:,2], self.osci[1][:,2] and self.osci[2][:,2]. I. e. self.osci[:2][:,2]
  def objective(self, params, x, data):
    #read out of params w and f in arrays
    w = np.array([])
    f = np.array([])
    for j in range(len(self.guess[0,1:])):
      off = int(j*len(self.guess[:,0]))
      for i in range(len(self.guess[:,0])):
        w = np.append(w,params['w_%i' % (i+1+off)])
        f = np.append(f,params['f_%i' % (i+1+off)])
    
    #make residual per data set
    resid = data - self.fit_sinc(w, f, x)

    #now flatten this to a 1D array, as minimize() needs
    return resid
    
