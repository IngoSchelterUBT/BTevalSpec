#This contains the Data class where all the information for the calculation is stored
import numpy as np
import configparser as cpar

class Data:

  def __init__(self):
    ##############################################################
    #Data of config-file
    #boolean for fourier transformation
    self.fourier = False
    #Number of calculations the fourier transform should be made of (1 or 3)
    self.numCalc = 0 
    #Kaiser-Bessel window for fourier transformation
    self.window = 0.
    #smooth parameter for fourier tranformation
    self.smooth = 0.
    #default value for min. length of fourier tranform vector. Can increase the sampling rate of the ft. (default = 17)
    self.min_pw = 17 

    #Data of dipole file
    #number of atoms
    self.natoms = 0
    #number of atom types
    self.ntypes = 0
    #array of number of atoms per type
    self.atomspertype = np.array([])
    #number of electrons for spin-up and spin-down
    self.nelec = np.zeros(2,dtype=int)
    #number of all valence electrons
    self.nval = 0
    #spacing of the used grid in bohr
    self.dx = 0
    #start time of propagation in rydberg units
    self.t_start = 0.
    #length of time step in propagation
    self.dt = 0.
    #boost energy of boost calculation in rydberg
    self.boostenergy = 0.
    #k-vector of boost excitation self.kvec = np.zeros(3,dtype=np.float64)

    ##############################################################
    #Data for Pade Approximation
    #Data of config file
    #boolean for Pade Approximation
    pade = False
    #maximum value of Pade-Approximation in Ry (default=0.6)
    pade_wmax = 0.6
    #PADE_DW (default=0.00001)
    pade_dw = 0.00001
    #smooth parameter for Pade-Approximation 
    pade_smooth = 0.0000001
    #only keep every 2^n th data point to avoid running out of memory in pade program
    pade_thin = 0
  
    ##############################################################
    #Data for Fit calculation


#---------------------------------------------------------------------#
#Routine for creating config file
#---------------------------------------------------------------------#

  def createConfig(self):
    config = cpar.ConfigParser()
    #############################################################
    #Config part for fourier transformation
    config['FourierTransform'] = {}
    ft = config['FourierTransform']
    ft['fourier'] = 'yes #should fourier tranformation be made?'
    ft['numCalc'] = '0 #one calculation or three?'
    ft['window'] = '0 #Kaiser-Bessel windowing parameter (>0)'
    ft['smooth'] = '0 #Artificial decay rate'
    ft['min_pw'] = '17 #Min. length of fourier transform vector. Can increase the sampling rate of the ft'
    
   
    #############################################################
    #Config part for fourier transformation
    config['PadeApproximation'] = {}
    pade = config['PadeApproximation']
    pade['pade'] = 'no #should Pade-Approximation be made?'
    pade['pade_wmax'] = '0.6 #maximum energy in ry of pade approximation'
    pade['pade_dw'] = '0.00001 #Pade_DW'
    pade['pade_smooth'] = '0.0000001 #Smooth Parameter for Pade-Approximation'
    pade['pade_thin'] = '0 #only keep every 2^n th data point to avoid running out of memory in pade program'
    
    #############################################################
    #Config part for fitting the spectrum
    config['Fit'] = {}
    fit = config['Fit']
    fit['fit'] = 'no #should the spectrum be fitted?'

    #Output of config-File
    with open('eval.ini', 'w') as configfile:
      config.write(configfile)

#---------------------------------------------------------------------#
#Routine for reading config file
#---------------------------------------------------------------------#

  def readConfig(self):
    config = cpar.ConfigParser(inline_comment_prefixes='#')
    config.read('eval.ini')

    #############################################################
    #Read part for fourier transformation
    ft = config['FourierTransform']
    self.fourier = ft.getboolean('fourier',fallback=self.fourier)
    self.numCalc = ft.getint('numCalc',fallback=self.numCalc)
    self.window = ft.getfloat('window',fallback=self.window)
    self.smooth = ft.getfloat('smooth',fallback=self.smooth)
    self.min_pw = ft.getint('min_pw',fallback=self.min_pw)

    #############################################################
    #Read part for Pade-Approximation
    pade = config['PadeApproximation']
    self.pade = pade.getboolean('pade',self.pade)
    self.pade_wmax = pade.getfloat('pade_wmax',fallback=self.pade_wmax)
    self.pade_dw = pade.getfloat('pade_dw',self.pade_dw)
    self.pade_smooth = pade.getfloat('pade_smooth',self.pade_smooth)
    self.pade_thin = pade.getint('pade_thin',self.pade_thin)

#---------------------------------------------------------------------#
#Routine for reading dipole file
#---------------------------------------------------------------------#

  def readDipole(self):
    
