#[1] A. Bruner, D. LaMaster, K. Lopata, "Accelerated Broadband Spectra Using i
#Transition Dipole Decomposition and Pade Approximants", JCTC 12, 3741-3750 (2016)

#Mathtools module for calculating for loops fast

#Import Libraries
import numpy as np
import numba as nb


#Routine for calculating the PADE-series
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
