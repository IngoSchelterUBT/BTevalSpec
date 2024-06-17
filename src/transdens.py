#==============================================================================#
# Transition-density decoupling
#==============================================================================#
import errorHandler as err
import numpy as np

#==============================================================================#
#in:
#   densft: List containing the Fourier-Transformed densities at the energies given in densex 
#   densen: energies at which the Fourier-Transformed densities are given; ideally these are the excitation energies
#   excit:  List of excitations (class Excitations)
#   T:      Free-propagation time
#   Ef:     Field strength
#   Ep:     Field polarization
#   Hw:     List of field-profile at given excitation energies
#   jcalc:  Calculation the densft belong to
#
# out:
#   transdens: List containing the Transition densities (as 1d arrays on the grid)
#==============================================================================#
def decouple(densft,densen,excit,T,Ef,Ep,Hw,jcalc=0,dbg=0,imagonly=False):

    # Check
    nex = len(excit.exlist)
    if any(nex!=nn for nn in [len(densft),len(densen)]): err.err(2,"Number of densities and energies must match number of excitations")

    # Compute matrix
    b = np.zeros((nex,nex))
    for iex, ex in enumerate(excit.exlist):
        c= -Ef*np.dot(Ep,ex.dipole)*np.abs(Hw[iex])#*1/hbar, which is one in Ry a.u.
        for ien, en in enumerate(densen):
            wm     = en-ex.energy
            wp     = en+ex.energy
            sincm  = np.sinc(wm*T/np.pi) #np.sinc is defined as sin(pi*x)/(pi*x)
            coscm  = (1.-np.cos(wm*T))/(wm*T) if abs(wm)>0. else 0.
            sincp  = np.sinc(wp*T/np.pi)
            coscp  = (1.-np.cos(wp*T))/(wp*T) if abs(wp)>0. else 0.
            b[ien,iex] = c*(\
                np.exp(-1.0j*ex.phase)*T*(coscp + 1.j*sincp) -\
                np.exp(+1.0j*ex.phase)*T*(coscm + 1.j*sincm))
    if imagonly: b=np.imag(b) #This makes B real-valued -> also use (real-valued) imag part of densft below

    #Verbose output
    if dbg>0:
        print("")
        print("Transition-density correlation matrix")
        print(b)

    # Invert matrix
    binv = np.linalg.inv(b)

    #Verbose output
    if dbg>0:
        print("")
        print("Inverted matrix")
        print(binv)

    # Apply inverted matrix
    transdens = []
    for iex, ex in enumerate(excit.exlist):
        transdens.append(np.empty(np.shape(densft[0]),dtype=np.cdouble))
        for ien, dens in enumerate(densft):
            if imagonly:
                transdens[iex] += binv[iex,ien]*np.imag(dens) #Check index order
            else:
                transdens[iex] += binv[iex,ien]*        dens
    return transdens

#==============================================================================#
# Compute the dipole of a given density
#==============================================================================#
def getDipole(dens,org,dxyz):
    nn  = np.shape(dens)
    dip = np.zeros((len(dens),3))
    for n in range(nn[0]):
        for ix in range(nn[1]):
            for iy in range(nn[2]):
                for iz in range(nn[3]):
                    dip[n,0] += dens[n,ix,iy,iz]*(org[0]+ix*dxyz[0])
                    dip[n,1] += dens[n,ix,iy,iz]*(org[1]+iy*dxyz[1])
                    dip[n,2] += dens[n,ix,iy,iz]*(org[2]+iz*dxyz[2])
                    norm += dens[n,ix,iy,iz]
    return dip

