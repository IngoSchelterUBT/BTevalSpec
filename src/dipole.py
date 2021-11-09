#This contains the Data class where all the information for the calculation is stored
import numpy as np

#own module
import errorHandler as err
import inout

class Dipole:

  def __init__(self,fileName,dipoleId):
    if isinstance(fileName, list) and len(fileName) == 3:
      #------------------------------------------------------------#
      ##############################################################
      # Read all three dipole files for building trace
      ##############################################################
      print(len(fileName))

    else:
      #------------------------------------------------------------#
      ##############################################################
      # Read Diople-File for creating dipole object
      ##############################################################
      dipoleFile = inout.readDipole(fileName)
      #ID for the dipole file, if there are multiple dipole files
      self.dipoleId = dipoleId
      
      #Data of dipole file
      #number of atoms
      self.natoms = int(dipoleFile.get('NATOMS',[0])[0])
      if self.natoms == 0:
        err.err(1,'There are no atoms in dipole file!')
      #number of atom types
      self.ntypes = int(dipoleFile.get('NTYPES',[0])[0])
      if self.ntypes == 0: 
        err.err(1,'There are no atom types in dipole file!')
      #array of number of atoms per type
      self.atomspertype = dipoleFile.get('ATOMSPERTYPE','empty')
      if 'empty' in str(self.atomspertype):
        err.err(1,'There are no atomtypes in dipole file!')
      #total number of electrons, electrons for spin-up and spin-down (1., 2. and 3. position)
      self.nelec = dipoleFile.get('NELEC','empty')
      if 'empty' in str(self.nelec):
        err.err(1, 'There are no electrons in spin-up and spin-down in dipole file!')
      totelec = sum(self.nelec)
      self.nelec = np.insert(self.nelec, 0, totelec, axis=0)
      #spacing of the used grid in bohr
      self.dxyz = dipoleFile.get('DX',[0.])[0]
      if self.dxyz == 0.:
        err.err(1, 'There is no grid spacing in dipole file!')
      #start time of propagation in rydberg units
      self.t_start = dipoleFile.get('T_START',[0.])[0]
      #length of time step in propagation
      self.dt = dipoleFile.get('DT',[0.])[0]
      if self.dt == 0.:
        err.errr(1, 'There is no time step in dipole file!')
      #boost energy of boost calculation in rydberg
      self.boostenergy = dipoleFile.get('BOOSTENERGY',[0.])[0]
      #k-vector of boost excitation 
      self.kvec = dipoleFile.get('BOOSTKVEC','empty')
      #dipole moment for dipole file
      self.dipData = np.loadtxt(fileName,comments='#')
