#[1] A. Bruner, D. LaMaster, K. Lopata, "Accelerated Broadband Spectra Using i
#Transition Dipole Decomposition and Pade Approximants", JCTC 12, 3741-3750 (2016)



#File for calculation the fourier transformation of the input data (dipole file(s))

import yaml
import os, glob
import numpy as np
from scipy.fft import fft,fftfreq
from scipy import signal

#import own modules
import errorHandler as err
#import f90 module
import f90_tools.mathtools as mathtools

class Pade:
  
  def __init__(self,config,dipole,Id):
    #Save all the information needed for the fit, either from dipole files or from reading the fourier transformed dipole moment files
    #if ft should be transformed than the values for the fourier transformation are also needed

    ################################
    #information out of dipole file
    #id of Pade Approximation
    self.padeId = Id

    ################################
    if config.pade:
      #read propagation time and k-vecktor out of dipole
      self.propTime = dipole.dipData[len(dipole.dipData)-1,0]
      self.kvec = dipole.kvec

      #calculate fourier transformation and pw-spectrum with dipole object
      self.calcPade(config,dipole)

      #write Pade_Osci_* files in PADE directory
      self.writePade()
    elif not config.pade and config.fit_guess:
      #read Pade_Osci_* files in PADE directory
      self.readPadeOsci()



#-----------------------------------------------------------------------------#
#   Methods of the class Pade
#-----------------------------------------------------------------------------#
  #Routine for calculation the fourier transformation and the PW
  def calcPade(self,config,dipole):
    time = dipole.dipData[:,0] 
    func = dipole.dipData[:,1:]

    #Get number of frequencies
    freq = self.getPadeFreq(config.pade_wmax,config.pade_dw)

    #calculate the ft of the time for Osci
    self.padeOsci = [] #padeOsci as list of ft (freq,PadeAprox) for x-,y- and z-component

    #dipole moments
    #variable to transform Ry to Osci
    Ry2Osci = -np.sqrt(dipole.nelec[0])/np.sqrt(dipole.boostenergy)/np.sqrt(2)/3.

    #calculate Pade Approx
    for i in range(len(func[0,:])):
      mu = self.padeData(config,time,freq,func[:,i])
      mu = mu*dipole.dt*Ry2Osci
      mu_real = np.real(mu)
      mu_imag = np.imag(mu)
      self.padeOsci.append(np.column_stack([freq,mu_imag]))



      
  #Routine for calculation the frequencies in Pade
  def getPadeFreq(self,pade_wmax,pade_dw):
    w0 = 0.
    wn = int(np.ceil((float(pade_wmax)-float(w0))/pade_dw)+1)
    w = np.array([])
    for i in range(wn):
      w = np.append(w,w0 + i*pade_dw)
    return w


  #Routine for calculation of Pade Apprximation of the dipole moment
  def padeData(self,config,time,freq,dip):
    w = freq
    wn = len(freq)
    #Set dimension of the problem
    m = int(len(time))-1
    m = int(2*(m/2)) #ensure m to be even by integer division
    n = int(m/2)

    #Remove DC component of the dipole moment
    dip = dip - np.sum(dip)/(m+1)

    #Get time step and propagation time
    t_start = time[0]
    t_stop = time[m]
    t_prop = t_stop - t_start
    dt = t_prop/m

    #Get decay rate of the exponential function and apply to the dipole data
    #eta=-log(decayFraction)/t_prop
    if config.pade_smooth > 0.:
      for i in range(m+1):
        dip[i] = dip[i]*np.exp(-config.pade_smooth*(time[i]-t_start))

    #Get matrix G and vector d eq~(33) in [1]
    b0 = 1.0
    G = np.zeros((n,n))
    d = np.zeros(n)
    for k in range(n):
      for i in range(n):
        G[k,i] = dip[n-(i+1)+(k+1)]
      d[k] = -dip[n+(k+1)]

    #Solve G*b=d eq~(34) in [1] (resulting b are stored in d-array)
    b = np.linalg.solve(G,d) #is the solution of the LGS
    
    #Calculate a-coefficients eq~(35) in [1]
    a = np.zeros(n)
    a0 = dip[0]
    for k in range(n):
      a[k] = b0*dip[k+1]
      for i in range(k):
        a[k] = a[k] + b[i]*dip[(k+1)-(i+1)]

    #Pade approximation for all values of w (using fortran)
    mu = mathtools.padeseries(w, wn, n, dt, a0, b0, a, b)

    return mu

  #Routine for writing the fourier transformations (osci) and pw-spectrums to file
  def writePade(self):
    #PadeOsci
    headPadeOsci = 'Energy (Ry) | imag part of pade approximation'
    #create folder PADE, if it does not exist or delete all 'Pade_Osci_Id'-files in directory
    if not os.path.exists('PADE'):
      os.makedirs('PADE')

    #Save all the Pade_Osci-files in directory
    for i in range(len(self.padeOsci)):
      np.savetxt('PADE/Pade_Osci_' + str(self.padeId) + self.getDir(i),self.padeOsci[i],header=headPadeOsci)


  #Routine for reading the osci-files
  def readPadeOsci(self):
    try:
      PadeOsciFile = open('PADE/Pade_Osci_' + str(self.padeId) + 'x','r')
    except:
      err.err(1,('You are trying to do a fit_guess (with no PADE-Approximation) but there' +
                ' is no corresponding Pade-file in the PADE directory'))
    #Read the fourier transformation itself
    self.padeOsci = []

    for i in range(3):
      self.padeOsci.append(np.loadtxt('PADE/Pade_Osci_' + str(self.padeId) + self.getDir(i)))

  #Get direction
  def getDir(self,i):
    if i == 0:
      return 'x'
    elif i == 1:
      return 'y'
    elif i == 2:
      return 'z'
    else:
      err.err(1,'There can not be more directions than x, y and z!')
