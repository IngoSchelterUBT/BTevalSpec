#Import Libraries
import numpy as np
import os
from scipy.fft import fft,fftfreq,fftshift
from scipy     import interpolate

#own module
import errorHandler as err

#==============================================================================#
# External excitation time profile: h(t) -> H(w)
#   maxdw: Minimum energy resolution (choice)
#==============================================================================#
class Extern:
    def __init__(self,fname,efield,epol,maxdw=0.0001):
        self.fname  = fname
        self.efield = efield
        self.epol   = epol #List of efield-polarizations for each calculation
        if (self.fname!=""): #laser
            dat         = np.transpose(np.loadtxt(fname,comments='#'))
            self.time   = dat[0] # Keep original time frame
            self.ext    = dat[1]
            self.dt     = (self.time[-1]-self.time[0])/(len(self.time)-1)
            self.ncalc  = len(self.epol)
            pw          = np.int(np.ceil(np.log2(2.*np.pi/self.dt/maxdw)))
            Nf          = 2**pw

            self.freq   = fftfreq(Nf,d=self.dt)*2.*np.pi
            self.ft     = self.dt*fft(self.ext,n=Nf)
            self.ftint  = interpolate.interp1d(self.freq,self.ft,kind="cubic")

    #----------------------------------------------------------------------------#
    # Write spectra
    #----------------------------------------------------------------------------#
    def write(self):
        if (self.fname!=""): #laser
            head = 'Energy (Ry) | real(FT) | imag(FT)'
            fname = os.path.splitext(self.fname)[0]+'_ft.dat'
            np.savetxt(fname,np.column_stack((fftshift(self.freq),np.real(fftshift(self.ft)),np.imag(fftshift(self.ft)))),header=head)

    #----------------------------------------------------------------------------#
    # Return interpolated value
    #----------------------------------------------------------------------------#
    def getVal(self,w):
        if self.fname!="": #laser
            try:
                ft = [self.ftint(w[i]) for i in range(len(w))]
            except ValueError:
                print("Value Error in extern interpolation (e.g., energy out of range)")
                raise
        else: #boost
            ft = [1.+0.j]*len(w)
        return ft
