#[1] A. Bruner, D. LaMaster, K. Lopata, "Accelerated Broadband Spectra Using i
#Transition Dipole Decomposition and Pade Approximants", JCTC 12, 3741-3750 (2016)

#Import Libraries
import os
import numpy as np
import numba as nb

#import own modules
import errorHandler as err
import util
import mathtools

#==============================================================================#
# Class Pade
#==============================================================================#
# The Class consists of the following variables:
# - padeId: Id of Pade Approximation (equals 3 for trace)
# - propTime: Propagation time of calculation
# - kvec: k-vector of excitation at beginning of calculation
# - padeOsci: List of Pade Approximations for x-, y- and z-direction of dipole
#             moment (energy | padeApprox)
#------------------------------------------------------------------------------#

class Pade:
  def __init__(self,config,dipole,Id,calcFlag='no'):
    #Save all the information needed for the fit, either from dipole files or from reading the fourier transformed dipole moment files
    #if ft should be transformed than the values for the fourier transformation are also needed

    ################################
    #information out of dipole file
    #id of Pade Approximation
    self.padeId = Id

    ################################
    if config.pade or calcFlag == 'guess':
      #read propagation time and k-vecktor out of dipole
      self.propTime = dipole.dipData[len(dipole.dipData)-1,0]
      self.kvec = dipole.kvec

      #calculate fourier transformation and pw-spectrum with dipole object
      self.calcPade(config,dipole)

      #write Pade_Osci_* files in PADE directory
      self.writePade(calcFlag)
    elif not config.pade and config.fit_guess:
      #read Pade_Osci_* files in PADE directory
      self.readPadeOsci(calcFlag)



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

    #Get time step and propagation time
    t_start = time[0]
    t_stop = time[m]
    t_prop = t_stop - t_start
    dt = t_prop/m

    #Thin: Remove every 2nd line from the file starting at line 2
    if config.pade_thin > 0:
      for i in range(config.pade_thin):
        time = np.delete(time, np.arange(0, time.size, 2))
        func = np.delete(func, np.arange(0, func.size, 2))
        dt = 2.*dt

    #Remove DC component of the dipole moment
    dip = dip - np.sum(dip)/(m+1)

    #Get decay rate of the exponential function and apply to the dipole data
    #eta=-log(decayFraction)/t_prop
    if config.pade_smooth > 0.:
      dip = dip*np.exp(-config.pade_smooth*(time-t_start))
      #for i in range(m+1):
      #  dip[i] = dip[i]*np.exp(-config.pade_smooth*(time[i]-t_start))
    
    mu = mathtools.numba_padeseries(w, wn, m, n, dt, dip)

    return mu


  #Routine for writing the fourier transformations (osci) and pw-spectrums to file
  def writePade(self,calcFlag='no'):
    #PadeOsci
    headPadeOsci = 'Energy (Ry) | imag part of pade approximation'
    #create folder PADE, if it does not exist or delete all 'Pade_Osci_Id'-files in directory
    if not os.path.exists('PADE'):
      os.makedirs('PADE')

    if calcFlag == 'no':
      #Save all the Pade_Osci-files in directory
      for i in range(len(self.padeOsci)):
        np.savetxt('PADE/Pade_Osci_' + str(self.padeId + 1) + util.getDir(i),self.padeOsci[i],header=headPadeOsci)
    elif calcFlag == 'trace':
      #Save trace Pade-file in directory
      for i in range(len(self.padeOsci)):
        np.savetxt('PADE/Pade_Osci',self.padeOsci[i],header=headPadeOsci)


  #Routine for reading the osci-files
  def readPadeOsci(self,calcFlag):
    try:
      if calcFlag == 'no':
        PadeOsciFile = open('PADE/Pade_Osci_' + str(self.padeId + 1) + 'x','r')
      elif calcFlag == 'trace':
        PadeOsciFile = open('PADE/Pade_Osci','r')
      elif calcFlag == 'guess':
        return
    except:
      err.err(1,('You are trying to do a fit_guess (with no PADE-Approximation) but there' +
                ' is no corresponding Pade-file in the PADE directory'))
    #Read the fourier transformation itself
    self.padeOsci = []

    if calcFlag == 'no':
      for i in range(3):
        self.padeOsci.append(np.loadtxt('PADE/Pade_Osci_' + str(self.padeId + 1) + util.getDir(i)))
    elif calcFlag == 'trace':
      self.padeOsci.append(np.loadtxt('PADE/Pade_Osci'))
