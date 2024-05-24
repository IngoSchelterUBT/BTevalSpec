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
def decouple(densft,densen,excit,T,Ef,Ep,Hw,jcalc=0,dbg=0):

    # Check
    nex = len(excit.exlist)
    if any(nex!=nn for nn in [len(densft),len(densen)]): err.err(2,"Number of densities and energies must match number of excitations")

    # Compute matrix
    b = np.zeros((nex,nex))
    for iex, ex in enumerate(excit.exlist):
        c= -Ef*np.dot(Ep,ex.dipole)*np.abs(Hw[iex])#*1/hbar, which is one in Ry a.u.
        for ien, en in enumerate(densen):
            wm     = en-ex.energy
            wp     = en-ex.energy
            sincm  = np.sinc(wm*T/np.pi) #np.sinc is defined as sin(pi*x)/(pi*x)
            coscm  = (1.-np.cos(wm*T))/(wm*T) if abs(wm)>0. else 0.
            sincp  = np.sinc(wp*T/np.pi)
            coscp  = (1.-np.cos(wp*T))/(wp*T) if abs(wp)>0. else 0.
            b[ien,iex] = c*(\
                np.exp(-1.0j*ex.phase)*T*(coscp + 1.j*sincp) -\
                np.exp(+1.0j*ex.phase)*T*(coscm + 1.j*sincm))

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
        transdens.append(np.zeros(np.shape(densft[0])))
        for ien, dens in enumerate(densft):
            transdens[iex] += binv[iex,ien]*dens #Check index order

    return transdens
