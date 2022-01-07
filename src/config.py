#This contains the Data class where all the information for the calculation is stored
import numpy as np

#own modules
import inout

#TODO: make errors for a false configuration!!

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
    #boolean if fit should be calculated
    self.fit = configFile.get('OPT').get('FitSpectrum').get('fit',False)
    #boolean if guess for fit out of Pade approximation should be done (or read out of config file if False)
    self.fit_guess = configFile.get('OPT').get('FitSpectrum').get('fit_guess',False)
    #criterium for fit, i.e. absolute deviation between raw data and fit
    self.fit_relerr_crit = configFile.get('OPT').get('FitSpectrum').get('fit_relerr_crit',0.1)
    #maximum numer of iterations for reaching the desired absolute deviation between raw data and fit
    self.fit_max_iter = configFile.get('OPT').get('FitSpectrum').get('fit_max_iter',1)
    #range of fit in spectrum
    self.fit_range = np.sort(configFile.get('OPT').get('FitSpectrum').get('fit_range',np.array([0])))
    #relative threshold (according to the maximum line) which line in the Pade approximation should be identified as an excitation line
    self.guess_thres = configFile.get('OPT').get('FitSpectrum').get('guess_thres',0.1)

    ##############################################################
    #Config for Fit output/guess
    #Read list of excitations as class
    self.excitations = self.Excitations(configFile,self.numDipoleFiles)

  class Excitations:

    def __init__(self,configFile,numDipoleFiles):
      excitations = configFile.get('SPEC','none')
      if excitations == 'none':
        return
      else:
        excitations = excitations.get('excitations',[])
      self.names = []
      self.energies = np.array([])
      self.osciStrengths = np.array([])
      self.phases = np.array([])
      self.amplitudes = [] #amplitudes as list of numpy arrays
      self.fix = np.array([])
      for i in range(len(excitations)):
        self.names.append(excitations[i].get('name','none'))
        self.energies = np.append(self.energies,excitations[i].get('energy'))
        self.osciStrengths = np.append(self.osciStrengths,excitations[i].get('strength',0.0))
        self.fix = np.append(self.fix,excitations[i].get('fix',False))
        #self.phases = np.append(self.phases,excitations[i].get('phase'))
        for j in range(numDipoleFiles):
          self.amplitudes.append(excitations[i].get('amplitude_file%i' % (j+1),np.array([])))
