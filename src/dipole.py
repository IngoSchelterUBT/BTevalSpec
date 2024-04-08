#Import Libraries
import numpy as np
import os
from scipy.fft import fft,fftfreq,fftshift
from scipy import signal

#own module
import errorHandler as err
import mathtools
import util

#==============================================================================#
# Class Dipole
#==============================================================================#
# The Class consists of the following variables:
#
# General
# - dipname     file name of input dipole file
# - descript    description list for each component, e.g., "['x','y','z']"
# Time-dependent data
# - time        numpy array containing the free-propagation-time axis
# - dat         contains dipole moment(s) on during the free-propagation time
# Frequency-dependent data (may be added)
# - freq        Frequencies of Fourier transform
# - ft          Fourier transform of dat
# - osci        Dipole strength of dat
# - pw          Power spectrum of dat
# - pade...
# Further components
# - meta        Meta-information dictionary
# - nelec       Number of valence electrons
# - tprop       Free propagation time
# - dt          Time step of propagation
# - ext         "boost" or "laser"
# - efield      electric field strength in Ry (boost: from boost energy)
# - epol        electric field polarization vector
# - text        End-time of ext. laser (or 0. in case of boost)
#               excitation period is removed from the dipole moment
#------------------------------------------------------------------------------#
class Dipole:
    #--------------------------------------------------------------------------#
    # Read Diople-File for creating dipole object
    #--------------------------------------------------------------------------#
    def __init__(self,fname,descript=[]):

        #Read and interpred meta data
        self.dipname = fname
        self.meta    = self.readDipoleHeader(fname)
        self.nelec   = sum(self.meta.get('NELEC',[0.,0.])) #number of electrons
        if self.meta.get("BOOSTMASK",["no"])[0] != "no":
            self.ext    = "boost"
            energy      = self.meta.get('BOOSTENERGY',[0.])[0] # En = N*hbar^2*k^2/2/m => k = sqrt(2*m*En/N/hbar^2) -> Ry: k=sqrt(En/N)
            self.efield = np.sqrt(energy/self.nelec/2.) # identify hbar*k*1 = e*E*H(omega) => E = hbar*k/e -> Ry: E=k/sqrt(2) = sqrt(En/2/N)
            self.epol   = self.meta.get('BOOSTKVEC',['empty']) #k-vector of boost excitation
            self.text   = 0.
        elif self.meta.get("EXCITATION",["no"])[0] == "laser_mono":
            if self.meta.get("NLASER",[1])[0] > 1: err(1,"More than one laser not supported")
            self.ext    = "laser"
            self.efield =  self.meta.get('LASEREFIELD',[0.,""])[0]
            self.epol   = [self.meta.get('LASERPOLX',['empty'])[0],self.meta.get('LASERPOLY','empty')[0],self.meta.get('LASERPOLZ','empty')[0]]
            self.text   = self.meta.get('LASEREND',0.)[0]
        else:
            err(1,"Unknown excitation")

        #Read dipole moment data from dipole file
        dat = np.transpose(np.loadtxt(fname,comments='#'))
        if len(descript)==0:
            self.descript = [i for i in range(len(dat)-1)]
        elif len(descript)==len(dat)-1:
            self.descript = descript
        else:
            err(1,"Invalid length of description list")
        self.time  = np.array([dat[0][i]-self.text for i in range(len(dat[0])) if dat[0][i]>=self.text])
        self.tprop = self.time[-1]-self.time[0]
        self.dt    = self.tprop/(len(self.time)-1)
        self.dat   = []
        for j in range(1,len(dat)):
            self.dat.append(np.array([dat[j][i] for i in range(len(dat[0])) if dat[0][i]>=self.text]))

        #Init derived components that are not initially computed
        self.ft       = []
        self.pw       = []
        self.freq     = []
        self.pade     = []
        self.freqPade = []

    #--------------------------------------------------------------------------#
    # Read header information of the dipole file and returns it as dictionary
    #--------------------------------------------------------------------------#
    def readDipoleHeader(self,dipolefile):
        meta = {}
        with open(dipolefile,'r') as f:
            for line in f:
                if '#!BT' in line:
                    line_array = line.split()
                    if len(line_array) == 3:
                        (key, val) = (line_array[1],line_array[2])
                        meta[key] = [val]
                    if len(line_array) == 4:
                        (key, [val1, val2]) = (line_array[1],[line_array[3],line_array[2]])
                        meta[key] = [val1, val2]

            #Check units of Dipole File
            if meta.get('DX')[-1] != 'a0':
                    err.warn('Unit of DX unknown!')
            elif meta.get('T_START')[-1] != 't_ry':
                    err.warn('Unit of T_START unknown!')
            elif meta.get('DT')[-1] != 't_ry':
                    err.warn('Unit of DT unknown!')
            elif meta.get('BOOSTENERGY',["Ry"])[-1] != 'Ry':
                    err.warn('Unit of BOOSTENERGY unknown!')
            elif meta.get('LASERFREQ',["Ry"])[-1] != 'Ry':
                    err.warn('Unit of LASERFREQ unknown!')
            elif meta.get('LASERINT',["W/cm^2"])[-1] != 'W/cm^2':
                    err.warn('Unit of LASERINT unknown!')
            elif meta.get('LASEREFIELD',["E_ry"])[-1] != 'E_ry':
                    err.warn('Unit of LASEREFIELD unknown!')
            elif meta.get('LASERSTART',["t_ry"])[-1] != 't_ry':
                    err.warn('Unit of LASERSTART unknown!')
            elif meta.get('LASEREND',["t_ry"])[-1] != 't_ry':
                    err.warn('Unit of LASEREND unknown!')
            elif meta.get('DIP0')[-1] != 'dip_ry':
                    err.warn('unit of DIP0 unknown!')

            for k, v in meta.items():
                    if ',' in v[0]:
                            meta[k] = np.array([float(i) for i in v[0].split(',')]) #modifiziere v[0]: dort steht der value drinnen, d. h. wenn etwas Komma seperiertes ist dann aendere den value zu einer numpy list
                    if util.isFloat(v[0]):
                            v[0] = float(v[0])

            return meta

    #----------------------------------------------------------------------------#
    # Compute Fourier transform and power spectrum
    #----------------------------------------------------------------------------#
    def getFt(self,minpw=0,smooth=0.,window=0.,rmDC=True):

        # Frequencies
        pw = np.ceil(np.log2(len(self.time)))
        if minpw>0: pw = max(pw,minpw)
        Nf = 2**pw
        self.freq = fftfreq(Nf,d=self.dt)*2.*np.pi

        # Generate signal to transform
        tfunc = np.copy(self.dat)
        # Remove DC component
        if rmDC:
            for i in range(tfunc.shape[0]):
                tfunc[i] -= np.average(tfunc[i]) 
        # Smoothing
        if smooth>0.:
            for i in range(tfunc.shape[0]):
                tfunc[i] *= np.exp(-smooth*(self.time-self.time[0]))
        # Windowing
        if window>0.:
            for i in range(tfunc.shape[0]):
                tfunc[i] *= signal.windows.kaiser(len(tfunc[i]),window)

        # Fourier transform
        self.ft = []
        self.pw = []
        for i in range(len(tfunc)):
            tmp = self.dt*fft(tfunc[i],n=Nf )
            self.ft.append([np.real(tmp),np.imag(tmp)])
            self.pw.append(np.zeros(Nf//2,dtype=float))
            for j in range(len(self.pw[i])):
                self.pw[i][j] = np.abs(tmp[j])*np.abs(tmp[j]) + np.abs(tmp[-j])*np.abs(tmp[-j])

    #----------------------------------------------------------------------------#
    # Compute Pade approximation
    #----------------------------------------------------------------------------#
    def getPade(self,wmax,dw,thin=0,smooth=0.,window=0.,rmDC=True):

        # Frequencies
        wmax          = np.ceil(wmax/dw)*dw
        Nf            = int(np.rint(wmax/dw)+1)
        self.freqPade = np.linspace(0.,wmax,Nf)

        # Generate signal to transform
        tfunc = np.copy(self.dat)
        time  = np.copy(self.time)
        dt    = self.dt
        # Remove DC component
        if rmDC:
            for i in range(tfunc.shape[0]):
                tfunc[i] -= np.average(tfunc[i]) 
        # Smoothing
        if smooth>0.:
            for i in range(tfunc.shape[0]):
                tfunc[i] *= np.exp(-smooth*(self.time-self.time[0]))
        # Windowing
        if window>0.:
            for i in range(tfunc.shape[0]):
                tfunc[i] *= signal.windows.kaiser(len(tfunc[i]),window)
        # Thin out
        if thin>0:
            for j in range(thin):
                for i in range(tfunc.shape[0]):
                    tfunc = np.delete(tfunc[i],np.arange(0,tfunc[i].size, 2))
                time      = np.delete( time   , np.arange(0, time.size, 2))
                dt       *= 2.

        # Pade approximation
        m = ((len(time)-1)//2)*2 #ensure even m
        n = m//2
        self.pade = []
        for i in range(len(tfunc)):
            tmp = self.dt*mathtools.numba_padeseries(self.freqPade,m,n,dt,tfunc[i])
            self.pade.append([np.real(tmp),np.imag(tmp)])

    #----------------------------------------------------------------------------#
    # Write spectra
    #----------------------------------------------------------------------------#
    def writeSpectra(self,what=["ft","pw","pade"]):
        #FT
        if "ft" in what:
            head = 'Energy (Ry) | real(FT) | imag(FT)'
            try:
                for i in range(len(self.ft)):
                    fname = os.path.splitext(self.dipname)[0]+'_ft_'+self.descript[i]+'.dat'
                    np.savetxt(fname,np.column_stack((fftshift(self.freq),fftshift(self.ft[i][0]),fftshift(self.ft[i][1]))),header=head)
            except:
                err.warn("No FT data to write")

        #PW
        if "pw" in what:
            head = 'Energy (Ry) | Power Spectrum'
            try:
                for i in range(len(self.pw)):
                    fname = os.path.splitext(self.dipname)[0]+'_pw_'+self.descript[i]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq[:len(self.freq)//2],self.pw[i])),header=head)
                    
            except:
                err.warn("No PW data to write")

        #Pade
        if "pade" in what:
            head = 'Energy (Ry) | real(Pade) | imag(Pade) '
            try:
                for i in range(len(self.pade)):
                    fname = os.path.splitext(self.dipname)[0]+'_pade_'+self.descript[i]+'.dat'
                    np.savetxt(fname,np.column_stack((self.freqPade,self.pade[i][0],self.pade[i][1])),header=head)
            except:
                err.warn("No Pade data to write")

    #----------------------------------------------------------------------------#
    # Read spectra
    #----------------------------------------------------------------------------#
    def readSpectra(self,what=["ft","pw","pade"]):
        #FT
        if "ft" in what:
            try:
                self.ft = []
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_ft_'+self.descript[i]+'.dat'
                    tmp   = np.transpose(np.loadtxt(fname))
                    self.ft.append([fftshift(tmp[1]),fftshift(tmp[2])])
                    if i==0: self.freq = fftshift(tmp[0])
            except:
                err.err(1,"No FT files to read")

        #PW
        if "pw" in what:
            try:
                self.pw = []
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_pw_'+self.descript[i]+'.dat'
                    self.pw.append(np.transpose(np.loadtxt(fname))[1])
            except:
                err.err(1,"No PW files to read")

        #Pade
        if "pade" in what:
            try:
                self.pade = []
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_pade_'+self.descript[i]+'.dat'
                    tmp = np.transpose(np.loadtxt(fname))
                    self.pade.append([tmp[1],tmp[2]])
                    if i==0: self.freqPade = tmp[0]
            except:
                err.err(1,"No Pade files to read")
