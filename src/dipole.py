#Import Libraries
import numpy as np

#own module
import errorHandler as err
import inout

#==============================================================================#
# Class Dipole
#==============================================================================#
# The Class consists of the following variables:
# - dipoleId: Id of dipole file (equals 3 for trace)
# - dipData: Matrix of dipole moment data (time | x | y | z) for trace
#            (time | trace), for guess (time | x+y+z)
# - kvec: k-vector of excitation
# - natoms: Number of atoms of molecule
# - ntypes: Number of types of atoms of molecule
# - atomspertype: Number of atoms per type
# - nelec: Number of valence electrons [total, spin-up, spin-down]
# - dxyz: Grid spacing
# - t_start: Starttime of propagation
# - dt: Time step of propagation
# - boostenergy: Energy of boost excitation in Ry
#------------------------------------------------------------------------------#

class Dipole:
  def __init__(self,fileName,dipoleId,calcFlag='no'):
    self.calcFlag = calcFlag
    if isinstance(fileName, list) and len(fileName) == 3 and calcFlag == 'trace':
      #------------------------------------------------------------#
      ##############################################################
      # Read all three dipole files for building trace
      ##############################################################
      dipoleHeader = []
      kvecMatrix = np.empty((0,3), float)
      dipole = []
      for i, file in enumerate(fileName):
        dipoleHeader.append(self.readDipoleHeader(fileName[i]))
        kvec = dipoleHeader[i].get('BOOSTKVEC','empty')
        notzero = np.nonzero(kvec)[0][0]
        dipole.append(np.loadtxt(file,comments='#')[:,notzero+1])

      trace = np.zeros(len(dipole[0]))
      for i in range(len(dipole[0])):
        for j in range(len(dipole)):
          trace[i] = trace[i] + dipole[j][i]
      #print(kvecMatrix)
      #print(dipole[0])
      trace = np.column_stack([np.loadtxt(fileName[0])[:,0],trace])

      #Set information in dipole object and overwrite kvec
      self.setHeadInformation(dipoleHeader[0])
      self.dipData = trace
      self.kvec = 0 #set to default value for trace-dipole moment object

    elif calcFlag == 'no':
      #------------------------------------------------------------#
      ##############################################################
      # Read Diople-File for creating dipole object
      ##############################################################
      #dipoleFile ist the head information of the dipole file as dictionary
      dipoleFile = self.readDipoleHeader(fileName)
      #ID for the dipole object, if there are multiple dipole objects
      self.dipoleId = dipoleId

      #Set data of dipole file in variables
      self.setHeadInformation(dipoleFile)

      #dipole moment for dipole file
      self.dipData = np.loadtxt(fileName,comments='#')

    elif calcFlag == 'guess':
      #------------------------------------------------------------#
      ##############################################################
      # Build the sum of one dipole moment in x-, y- and z-direction
      ##############################################################
      #dipoleFile ist the head information of the dipole file as dictionary
      dipoleFile = self.readDipoleHeader(fileName)

      #ID for the dipole object
      self.dipoleId = dipoleId

      #Set data of dipole file in variables
      self.setHeadInformation(dipoleFile)

      #load dipole File
      data = np.loadtxt(fileName,comments='#')

      #safe as dipData the sum of x-, y- and z-component
      self.dipData = np.column_stack([data[:,0],np.sum(data[:,1:],axis=1)])


  #Routine reads the head-information of the dipolefile and returns it as dictionary
  def readDipoleHeader(self,dipolefile):
    dipole = {}
    with open(dipolefile,'r') as f:
      for line in f:
        if '#!BT' in line:
          line_array = line.split()
          if len(line_array) == 3:
            (key, val) = (line_array[1],line_array[2])
            dipole[key] = [val]
          if len(line_array) == 4:
            (key, [val1, val2]) = (line_array[1],[line_array[3],line_array[2]])
            dipole[key] = [val1, val2]

    #Check units of Dipole File
    if dipole.get('DX')[-1] != 'a0':
      err.warning(1,'Unit of DX unknown!')
    elif dipole.get('T_START')[-1] != 't_ry':
      err.warning(1,'Unit of T_START unknown!')
    elif dipole.get('DT')[-1] != 't_ry':
      err.warning(1,'Unit of DT unknown!')
    elif dipole.get('BOOSTENERGY')[-1] != 'Ry':
      err.warning(1,'Unit of BOOSTENERGY unknown!')
    elif dipole.get('LASERFREQ')[-1] != 'Ry':
      err.warning(1,'Unit of LASERFREQ unknown!')
    elif dipole.get('LASERINT')[-1] != 'W/cm^2':
      err.warning(1,'Unit of LASERINT unknown!')
    elif dipole.get('LASEREFIELD')[-1] != 'E_ry':
      err.warning(1,'Unit of LASEREFIELD unknown!')
    elif dipole.get('LASERSTART')[-1] != 't_ry':
      err.warning(1,'Unit of LASERSTART unknown!')
    elif dipole.get('LASEREND')[-1] != 't_ry':
      err.warning(1,'Unit of LASEREND unknown!')
    elif dipole.get('DIP0')[-1] != 'dip_ry':
      err.warning(1,'unit of DIP0 unknown!')

    for k, v in dipole.items():
      if ',' in v[0]:
        dipole[k] = np.array([float(i) for i in v[0].split(',')]) #modifiziere v[0]: dort steht der value drinnen, d. h. wenn etwas Komma seperiertes ist dann aendere den value zu einer numpy list
      if self.isfloat(v[0]):
        v[0] = float(v[0])

    return dipole

  #Routine checks if it is a float in a String
  def isfloat(self,value):
    try:
      float(value)
      return True
    except ValueError:
      return False

  #This routine sets the information of the header of the dipole file into variables
  #dipoleFile is a dictionary
  def setHeadInformation(self,dipoleFile):
    #Data of dipole file
    #number of atoms
    self.natoms = int(dipoleFile.get('NATOMS',[0])[0])
    if self.natoms == 0:
      err.warning(1,'There are no atoms in dipole file!')
    #number of atom types
    self.ntypes = int(dipoleFile.get('NTYPES',[0])[0])
    if self.ntypes == 0:
      err.warning(1,'There are no atom types in dipole file!')
    #array of number of atoms per type
    self.atomspertype = dipoleFile.get('ATOMSPERTYPE','empty')
    if 'empty' in str(self.atomspertype):
      err.warning(1,'There are no atomtypes in dipole file!')
    #total number of electrons, electrons for spin-up and spin-down (1., 2. and 3. position)
    self.nelec = dipoleFile.get('NELEC','empty')
    if 'empty' in str(self.nelec):
      err.err(1, 'There are no electrons in spin-up and spin-down in dipole file!')
    totelec = sum(self.nelec)
    self.nelec = np.insert(self.nelec, 0, totelec, axis=0)
    #spacing of the used grid in bohr
    self.dxyz = dipoleFile.get('DX',[0.])[0]
    if self.dxyz == 0.:
      err.warning(1, 'There is no grid spacing in dipole file!')
    #start time of propagation in rydberg units
    self.t_start = dipoleFile.get('T_START',[0.])[0]
    #length of time step in propagation
    self.dt = dipoleFile.get('DT',[0.])[0]
    if self.dt == 0.:
      err.warning(1, 'There is no time step in dipole file!')
    #boost energy of boost calculation in rydberg
    self.boostenergy = dipoleFile.get('BOOSTENERGY',[0.])[0]
    #k-vector of boost excitation
    self.kvec = dipoleFile.get('BOOSTKVEC','empty')
    if 'empty' in str(self.kvec):
      err.err(1, 'There was no k-vector in dipole file!')
