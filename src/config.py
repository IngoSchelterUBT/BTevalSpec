#This contains the Data class where all the information for the calculation is stored
import numpy as np

#own modules
import inout

class Config:

  def __init__(self):
    #------------------------------------------------------------#
    ##############################################################
    #Read Config-File
    ##############################################################
    configFile = inout.readConfig()
    # - configFile: dictionary for the configFile variables
    self.dipoleFiles = configFile.get('DIPOLE','dipole.dat')
    if isinstance(self.dipoleFiles, list):
      self.numDipoleFiles = len(self.dipoleFiles)
    else:
      self.numDipoleFiles = 1
    self.exitFile = configFile.get('EXCIT')
    
    ##############################################################
    #Data of configFile-file
    #boolean for fourier transformation
    self.fourier = configFile.get('OPT').get('FourierTransform').get('fourier',False)
    #Kaiser-Bessel window for fourier transformation
    self.window = configFile.get('OPT').get('FourierTransform').get('window',0.)
    #smooth parameter for fourier tranformation
    self.smooth = configFile.get('OPT').get('FourierTransform').get('smooth',0.)
    #default value for min. length of fourier tranform vector. Can increase the sampling rate of the ft. (default = 17)
    self.min_pw = int(configFile.get('OPT').get('FourierTransform').get('min_pw',17))


    ##############################################################
    #Config for Pade Approximation
    #boolean for Pade Approximation
    self.pade = configFile.get('OPT').get('PadeApprox').get('pade',False)
    #maximum value of Pade-Approximation in Ry (default=0.6)
    self.pade_wmax = configFile.get('OPT').get('PadeApprox').get('pade_wmax',0.6)
    #PADE_DW (default=0.00001)
    self.pade_dw = configFile.get('OPT').get('PadeApprox').get('pade_dw',1.0e-5)
    #smooth parameter for Pade-Approximation 
    self.pade_smooth = configFile.get('OPT').get('PadeApprox').get('pade_smooth',1.0e-7)
    #only keep every 2^n th data point to avoid running out of memory in pade program
    self.pade_thin = int(configFile.get('OPT').get('PadeApprox').get('pade_thin',0))
  
    ##############################################################
    #Config for Fit calculation
    self.fit = configFile.get('OPT').get('FitSpectrum').get('fit',False)
    self.fit_guess = configFile.get('OPT').get('FitSpectrum').get('fit_guess',False)
