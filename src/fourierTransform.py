#File for calculation the fourier transformation of the input data (dipole file(s))

import yaml
import os, glob
import numpy as np
from scipy.fft import fft,fftfreq
from scipy import signal

import errorHandler as err

class FT:
  
  def __init__(self,config,dipole,Id):
    #Save all the information needed for the fit, either from dipole files or from reading the fourier transformed dipole moment files
    #if ft should be transformed than the values for the fourier transformation are also needed

    ################################
    #information out of dipole file
    #id of fourier transform
    self.ftId = Id

    ################################
    if config.fourier:
      #read propagation time and k-vecktor out of dipole
      self.propTime = dipole.dipData[len(dipole.dipData)-1,0]
      self.kvec = dipole.kvec

      #calculate fourier transformation and pw-spectrum with dipole object
      self.calcFT(config,dipole)

      #write fourier transformation and pw-spectrum
      self.writeFT()
    elif not config.fourier and config.fit:
      #read propagation time and k-vector out of head Osci files
      #read fourier transformation out of Osci files, PW-Spectrum ist not needed
      self.readOsci()



#-----------------------------------------------------------------------------#
#   Methods of the class FT
#-----------------------------------------------------------------------------#
  #Routine for calculation the fourier transformation and the PW
  def calcFT(self,config,dipole):
    time = dipole.dipData[:,0] 
    func = dipole.dipData[:,1:]

    #calculating lenght of ft array as a power of 2
    pw2 = 0
    while 2**pw2 < len(time)-1:
      pw2 = pw2 + 1
    pw2 = max(pw2,config.min_pw)

    #calculate the ft of the time for Osci
    self.osci = [] #Osci as list of ft (freq,realFT,imagFT) for x-,y- and z-component

    #frequencies
    freq = self.getFreq(2**pw2,dipole.dt,False)
    freqPw = self.getFreq(2**(pw2-1)+1,dipole.dt,True)
    
    #smooth the signal if smooth parameter is greater than 0.
    if config.smooth > 0.:
      for i in range(func.shape[1]):
        for j in range(func.shape[0]):
          func[j,i] = func[j,i]*np.ext(-config.smooth*j*dipole.dt)

    #Kaiser-Bessel window if parameter is greater than 0.
    if config.window > 0.:
      for i in range(func.shape[1]):
        window = signal.windows.kaiser(len(func[:,i]),config.window)
        func[:,i] = func[:,i]*window

    #dipole moments
    #variable to transform Ry to Osci
    Ry2Osci = np.sqrt(dipole.nelec[0])/np.sqrt(dipole.boostenergy)/np.sqrt(2)/3.
    
    #Do fourier transformation for x-, y- and z-dipole moment
    for i in range(len(func[0,:])):
      ft = self.fft(func[:,i],2**pw2,dipole.dt)
      ft = ft*dipole.dt*Ry2Osci
      ft_real = np.real(ft)
      ft_imag = np.imag(ft)
      self.osci.append(np.column_stack([freq,ft_real,ft_imag]))

    #calculate the pw spectrum
    self.pw = []
    
    for i in range(func.shape[1]):
      powerSpec = self.pwSpec((self.osci[i][:,1]+1j*self.osci[i][:,2])/dipole.dt/Ry2Osci,
                              2**pw2,dipole.dt)
      self.pw.append(np.column_stack([freqPw,powerSpec]))

    np.savetxt('test_pw.dat',self.pw[0])


    #output of osci and pw with the header cosisting the propagation time and the kvector
      
      


  #Routine for calculating the fourier transformation of the time
  def getFreq(self,Nt,dt,pw):
    #Nt is the length of the vector
    #dt is the time step
    #pw (angular frequency), i.e. for Osci multiply with 2*pi or for PW not and
    #only positive frequencies
    if not pw:
      return fftfreq(Nt,d=dt)*2.*np.pi
    else:
      freqPW = np.zeros(int(Nt))
      for i in range(0,int(Nt)):
        freqPW[i] = float(i)/(Nt+1)/2/dt
      return freqPW


  #Routine for calculating the fourier transformation of the signal
  #it returns the real and imag value of the fft
  def fft(self,tfunc,N,dt):
    #tfunc is the time dependent function (i.e. dipole moment)
    #N is length of the ft-vector

    #drop DC
    dc = np.sum(tfunc)
    tfunc = tfunc - dc/float(len(tfunc))


    return fft(tfunc,n=N)

  #Routine for calculating the pw spectrum
  #it returns the pw spectrum
  def pwSpec(self,ft,N,dt):
    #ft is the fourier transformation (real and imag part)
    #N is the length of the ft array (2**pw2)
    #dt is the length of time step

    pw = np.zeros(int(N/2)+1)
    pw[0] = np.abs(ft[0])**2
    for i in range(1,int(N/2)+1):
      pw[i] = np.abs(ft[i])**2+np.abs(ft[int(N)-i])**2

    return pw

  #Routine for writing the fourier transformations (osci) and pw-spectrums to file
  def writeFT(self):
    #make Header
    head = self.makeHeader()

    #Osci
    headOsci = head + '\n' + 'Energy (Ry) | real part of FT | imag part of FT'
    #create folder Osci, if it does not exist or delete all 'Osci_Id'-files in directory
    if not os.path.exists('Osci'):
      os.makedirs('Osci')

    #Save all the Osci-files in directory
    for i in range(len(self.osci)):
      np.savetxt('Osci/Osci_' + str(self.ftId) + self.getDir(i),self.osci[i],header=headOsci)

    #PW
    headPW = head + '\n' + 'Energy (Ry) | power spectrum'
    #create folder PW, if it does not exist or delete all 'PW_Id'-files in directory
    if not os.path.exists('PW'):
      os.makedirs('PW')

    #Save all the PW-files in directory
    for i in range(len(self.pw)):
      np.savetxt('PW/PW_' + str(self.ftId) + self.getDir(i),self.pw[i],header=headPW)

  #Routine for making header as yaml-dictionary
  def makeHeader(self):
    #Save head information as dictionary
    headDict = {'PropTime(Ry)' : str(self.propTime), 'kVec' : self.kvec.astype(str).tolist()}
    
    headDictStr = yaml.dump(headDict)
    
    head = str()
    for i, line in enumerate(headDictStr.splitlines()):
      if i == 0:
        head += ('!BT ' + line)
      else:
        head += ('\n' + '!BT ' + line)
    
    return head



  #Routine for reading the osci-files
  def readOsci(self):
    #Read head of Osci-file in x-direction
    self.readHeaderOsci()

    #Read the fourier transformation itself
    self.osci = []

    for i in range(3):
      self.osci.append(np.loadtxt('Osci/Osci_' + str(self.ftId) + self.getDir(i)))

  #Routine for reading header information in osci-files with ftId (only x-direction)
  def readHeaderOsci(self):
    try:
      OsciFile = open('Osci/Osci_' + str(self.ftId) + 'x','r')
    except:
      err.err(1,('You are trying to do a fit (with no Fourier Trafo) but there' +
                ' is no corresponding Osci-file in the Osci directory'))
    head = str()
    for line in OsciFile.readlines():
      l = line.split()
      if l[1] == '!BT':
        head += (str(line))
    
    headYaml = str()
    for i, line in enumerate(head.splitlines()):
      if i == 0:
        headYaml += line[5:]
      else:
        headYaml += ('\n' + line[5:])

    headDict = yaml.load(headYaml,Loader=yaml.FullLoader)

    self.propTime = headDict.get('PropTime(Ry)',0.)
    if self.propTime == 0.:
      err(1,'There is no PropTime in Osci file!')
    self.kvec = headDict.get('kVec',0)
    if self.kvec == 0:
      err(1,'There is no kvec in Osci file!')
    self.kvec = np.array(self.kvec,dtype=float)

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
