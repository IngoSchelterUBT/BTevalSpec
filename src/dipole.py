#Import Libraries
import numpy as np
import math
import os
from scipy.fft import fft,fftfreq,fftshift
from scipy import signal

#own module
import errorHandler as err
import mathtools

#==============================================================================#
# Class Dipole
#==============================================================================#
# The Class consists of the following variables:
#
# General
# - dipname     file name of input dipole file
# Time-dependent data
# - time        numpy array containing the time axis
# - dat         contains dipole moment(s)
# Frequency-dependent data (may be added)
# - freq        Frequencies of Fourier transform
# - ft          Fourier transform of dat
# - osci        Dipole strength of dat
# - pw          Power spectrum of dat
# - pade...
# Further components
# - kvec:       k-vector of excitation
# - natoms      Number of atoms of molecule
# - ntypes      Number of types of atoms of molecule
# - atomspertype Number of atoms per type
# - nelec       Number of valence electrons [total, spin-up, spin-down]
# - dxyz        Grid spacing
# - t_start     Start time of propagation
# - dt          Time step of propagation
# - t_prop      Propagation time
# - boostenergy Energy of boost excitation in Ry
#------------------------------------------------------------------------------#
class Dipole:
    #--------------------------------------------------------------------------#
    # Read Diople-File for creating dipole object
    #--------------------------------------------------------------------------#
    def __init__(self,fname): #,dipoleId,calcFlag='no'):

        #Read and interpred meta data
        self.dipname = fname
        self.setHeadInformation(self.readDipoleHeader(fname))

        #Read data dipole moment for dipole file
        dat         = np.transpose(np.loadtxt(fname,comments='#'))
        self.time   = dat[0]
        self.dat    = dat[1:]
        self.t_prop = self.dt * (len(dat[0])-1)

    #--------------------------------------------------------------------------#
    # Routine checks if it is a float in a String
    #--------------------------------------------------------------------------#
    def isfloat(self,value):
        try:
            float(value)
            return True
        except ValueError:
            return False

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
                    err.warning(1,'Unit of DX unknown!')
            elif meta.get('T_START')[-1] != 't_ry':
                    err.warning(1,'Unit of T_START unknown!')
            elif meta.get('DT')[-1] != 't_ry':
                    err.warning(1,'Unit of DT unknown!')
            elif meta.get('BOOSTENERGY')[-1] != 'Ry':
                    err.warning(1,'Unit of BOOSTENERGY unknown!')
            elif meta.get('LASERFREQ')[-1] != 'Ry':
                    err.warning(1,'Unit of LASERFREQ unknown!')
            elif meta.get('LASERINT')[-1] != 'W/cm^2':
                    err.warning(1,'Unit of LASERINT unknown!')
            elif meta.get('LASEREFIELD')[-1] != 'E_ry':
                    err.warning(1,'Unit of LASEREFIELD unknown!')
            elif meta.get('LASERSTART')[-1] != 't_ry':
                    err.warning(1,'Unit of LASERSTART unknown!')
            elif meta.get('LASEREND')[-1] != 't_ry':
                    err.warning(1,'Unit of LASEREND unknown!')
            elif meta.get('DIP0')[-1] != 'dip_ry':
                    err.warning(1,'unit of DIP0 unknown!')

            for k, v in meta.items():
                    if ',' in v[0]:
                            meta[k] = np.array([float(i) for i in v[0].split(',')]) #modifiziere v[0]: dort steht der value drinnen, d. h. wenn etwas Komma seperiertes ist dann aendere den value zu einer numpy list
                    if self.isfloat(v[0]):
                            v[0] = float(v[0])

            return meta

    #--------------------------------------------------------------------------#
    # Sets the information of the header of the dipole file into variables
    # meta is a dictionary
    #--------------------------------------------------------------------------#
    def setHeadInformation(self,meta):
        #Data of dipole file
        #number of atoms
        self.natoms = int(meta.get('NATOMS',[0])[0])
        if self.natoms == 0:
            err.warning(1,'There are no atoms in dipole file!')
        #number of atom types
        self.ntypes = int(meta.get('NTYPES',[0])[0])
        if self.ntypes == 0:
            err.warning(1,'There are no atom types in dipole file!')
        #array of number of atoms per type
        self.atomspertype = meta.get('ATOMSPERTYPE','empty')
        if 'empty' in str(self.atomspertype):
            err.warning(1,'There are no atomtypes in dipole file!')
        #total number of electrons, electrons for spin-up and spin-down (1., 2. and 3. position)
        self.nelec = meta.get('NELEC','empty')
        if 'empty' in str(self.nelec):
            err.err(1, 'There are no electrons in spin-up and spin-down in dipole file!')
        totelec = sum(self.nelec)
        self.nelec = np.insert(self.nelec, 0, totelec, axis=0)
        #spacing of the used grid in bohr
        self.dxyz = meta.get('DX',[0.])[0]
        if self.dxyz == 0.:
            err.warning(1, 'There is no grid spacing in dipole file!')
        #start time of propagation in rydberg units
        self.t_start = meta.get('T_START',[0.])[0]
        #length of time step in propagation
        self.dt = meta.get('DT',[0.])[0]
        if self.dt == 0.:
            err.warning(1, 'There is no time step in dipole file!')
        #boost energy of boost calculation in rydberg
        self.boostenergy = meta.get('BOOSTENERGY',[0.])[0]
        #k-vector of boost excitation
        self.kvec = meta.get('BOOSTKVEC','empty')
        if 'empty' in str(self.kvec):
            err.err(1, 'There was no k-vector in dipole file!')

    #----------------------------------------------------------------------------#
    # Compute Fourier transform and power spectrum
    #----------------------------------------------------------------------------#
    def ft(self,minpw=0,smooth=0.,window=0.,rmDC=True):

        # Frequencies
        pw = math.ceil(math.log2(len(self.time)))
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
                tfunc[i] *= np.exp(-smooth*(self.time-self.t_start))
        # Windowing
        if window>0.:
            for i in range(tfunc.shape[0]):
                tfunc[i] *= signal.windows.kaiser(len(tfunc[i]),window)

        # Fourier transform
        self.ft = [] #np.zeros((len(self.dat),Nf   ),dtype=complex)
        self.pw = [] #np.zeros((len(self.dat),Nf//2),dtype=float  )
        for i in range(len(tfunc)):
            self.ft.append(self.dt*fft(tfunc[i],n=Nf ))
            self.pw.append(np.zeros(Nf//2,dtype=float))
            for j in range(len(self.pw[i])):
                self.pw[i][j] = np.abs(self.ft[i][j])*np.abs(self.ft[i][j]) + np.abs(self.ft[i][-j])*np.abs(self.ft[i][-j])

    #----------------------------------------------------------------------------#
    # Compute Pade approximation
    #----------------------------------------------------------------------------#
    def pade(self,wmax,dw,thin=0,smooth=0.,window=0.,rmDC=True):

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
                tfunc[i] *= np.exp(-smooth*(self.time-self.t_start))
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
        self.pade = [] #np.zeros((len(tfunc),Nf),dtype=complex)
        for i in range(len(tfunc)):
            self.pade.append(self.dt*mathtools.numba_padeseries(self.freqPade,m,n,dt,tfunc[i]))

    #----------------------------------------------------------------------------#
    # Write spectra
    #----------------------------------------------------------------------------#
    def writeSpectra(self,what=["ft","pw","pade"]):
        #FT
        if "ft" in what:
            head = 'Energy (Ry) | real(FT) | imag(FT)'
            try:
                for i in range(len(self.ft)):
                    fname = os.path.splitext(self.dipname)[0]+'_ft_'+str(i)+'.dat'
                    np.savetxt(fname,np.column_stack((fftshift(self.freq),np.real(fftshift(self.ft[i])),np.imag(fftshift(self.ft[i])))),header=head)
            except:
                err.warn("No FT data to write")

        #PW
        if "pw" in what:
            head = 'Energy (Ry) | Power Spectrum'
            try:
                for i in range(len(self.pw)):
                    fname = os.path.splitext(self.dipname)[0]+'_pw_'+str(i)+'.dat'
                    np.savetxt(fname,np.column_stack((self.freq[:len(self.freq)//2],self.pw[i])),header=head)
                    
            except:
                err.warn("No PW data to write")

        #Pade
        if "pade" in what:
            head = 'Energy (Ry) | real(Pade) | imag(Pade) '
            try:
                for i in range(len(self.pade)):
                    fname = os.path.splitext(self.dipname)[0]+'_pade_'+str(i)+'.dat'
                    np.savetxt(fname,np.column_stack((self.freqPade,np.real(self.pade[i]),np.imag(self.pade[i]))),header=head)
            except:
                err.warn("No Pade data to write")

    #----------------------------------------------------------------------------#
    # Read spectra
    #----------------------------------------------------------------------------#
    def readSpectra(self,what=["ft","pw","pade"]):
        #FT
        self.ft = []
        if "ft" in what:
            try:
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_ft_'+str(i)+'.dat'
                    tmp = np.transpose(np.loadtxt(fname))
                    self.ft.append(fftshift(tmp[1]) + 1.j*fftshift(tmp[2]))
                    if i==0: self.freq = fftshift(tmp[0])
            except:
                err.warn("No FT files to read")

        #PW
        if "pw" in what:
            try:
                self.pw = []
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_pw_'+str(i)+'.dat'
                    self.pw.append(np.transpose(np.loadtxt(fname))[1])
            except:
                err.warn("No PW files to read")

        #Pade
        self.pade = []
        if "pade" in what:
            try:
                for i in range(len(self.dat)):
                    fname = os.path.splitext(self.dipname)[0]+'_pade_'+str(i)+'.dat'
                    tmp = np.transpose(np.loadtxt(fname))
                    self.pade.append(tmp[1] + 1.j*tmp[2])
                    if i==0: self.freqPade = tmp[0]
            except:
                err.warn("No Pade files to read")
