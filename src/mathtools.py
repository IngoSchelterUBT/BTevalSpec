#[1] A. Bruner, D. LaMaster, K. Lopata, "Accelerated Broadband Spectra Using i
#Transition Dipole Decomposition and Pade Approximants", JCTC 12, 3741-3750 (2016)

#Mathtools module for calculating for loops fast

#Import Libraries
import numpy as np
import numba as nb
from scipy.signal import butter, lfilter, freqz

#------------------------------------------------------------------------------#
#Routine for calculating the PADE-series
#------------------------------------------------------------------------------#
@nb.njit
def numba_padeseries(w, m, n, dt, dip):
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
 
  wn = len(w)
  z = np.full(n,0.0+0.0j)
  p = np.full(wn,0.0+0.0j)
  q = np.full(wn,0.0+0.0j)
  mu = np.full(wn,0.0+0.0j)
  for wi in range(wn):
    #Get powers of the complex exponential function
    z0 = 1.0
    z1 = np.exp(1.0j*w[wi]*dt)
    z[0] = z1
    
    for i in range(1,n):
      z[i] = z[i-1]*z1

    #Create Pade approximation eq~(29) in [1]
    p[wi] = a0
    q[wi] = b0
    for k in range(n):
      p[wi] = p[wi]+a[k]*z[k]
      q[wi] = q[wi]+b[k]*z[k]

  mu = p/q
  return mu

#------------------------------------------------------------------------------#
# Return Fit spectrum 
# T                             Eff. Propagation time
# w[nf]                         Frequency samples
# wi[nex]                       Excitation Frequencies
# p[nex]                        Excitation Phases
# t[nf]                         Time modifier (should be very close to 1.)
# a[ncalc][narea][ncomp][nex]   Amplitudes
#------------------------------------------------------------------------------#
@nb.njit
def fspectrum(ncalc,narea,ncomp,rc,nf,T,w,wi,p,tm,a):
    nrc = len(rc)
    f = np.zeros((ncalc,narea,ncomp,nrc,nf)) #,dtype=float)
    for icalc in range(ncalc):
        for iarea in range(narea):
            for n in range(ncomp):
                for i in range(nf):
                    for iex in range(len(wi)):
                        wm     = w[i]-wi[iex]
                        wp     = w[i]+wi[iex]
                        t      = T*tm[iex]
                        sincm  = np.sinc(wm*t/np.pi) #np.sinc is defined as sin(pi*x)/(pi*x)
                        coscm  = (1.-np.cos(wm*t))/(wm*t) if abs(wm)>0. else 0.
                        sincp  = np.sinc(wp*t/np.pi)
                        coscp  = (1.-np.cos(wp*t))/(wp*t) if abs(wp)>0. else 0.
                        tmp    = a[icalc][iarea][n][iex]*(\
                            np.exp(-1.0j*p[iex])*t*(coscm + 1.j*sincm) -\
                            np.exp(+1.0j*p[iex])*t*(coscp + 1.j*sincp) )
                        tmp    = [np.real(tmp),np.imag(tmp)]
                        for irc in range(nrc):
                            f[icalc][iarea][n][irc][i] += tmp[rc[irc]]
    return f

##------------------------------------------------------------------------------#
## Butterworth lowpass filter
##------------------------------------------------------------------------------#
#def butter_lowpass(cutoff, fs, order=5):
#    return 
#
#def butter_lowpass_filter(data, cutoff, fs, order=5):
#    b, a = butter(order, cutoff, fs=fs, btype='low', analog=False)
#    b, a = butter(order, cutoff, fs=fs, btype="low", analog=False)
#    y = lfilter(b, a, data)
#    return y
