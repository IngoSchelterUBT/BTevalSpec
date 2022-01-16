#==============================================================================#
# Class Excitations
#==============================================================================#
# The Class consists of the following variables:
# - list: List of dictionaries of excitation lines
#------------------------------------------------------------------------------#

class Excitations:
  def __init__(self,config):
    #List of dictionieries of excitation lines
    self.list = [None]*config.excitations
    default_excite = self.defaultExcitation()

    for i, value in enumerate(self.list):
      self.list[i] = default_excite | config.excitations[i]

#-----------------------------------------------------------------------------#
#   Methods of the class Excitations
#-----------------------------------------------------------------------------#
  #Routine for adding a new dictionary for a new excitation line to the list
  def addNewLine(self,name='S',fix=False,energy=0.,strength=0.,
               phase=0.,amplitude1=[0.,0.,0.],amplitude2=[0.,0.,0.],
               amplitude3=[0.,0.,0.]):
    exct = {}
    exct['name'] = name
    exct['fix'] = fix
    exct['energy'] = energy
    exct['strength'] = strength
    exct['phase'] = phase
    exct['amplitudes_file1'] = amplitude1
    exct['amplitudes_file2'] = amplitude2
    exct['amplitudes_file3'] = amplitude3

    self.list.append(exct)

  #Routine for sorting excitations after energies
  def sortLinesEnergy(self):
    self.list = sorted(self.list, key=lambda x: x['energy'])

  #Routine for deleting entry of list of dictionary
  def deleteKey(self,key):
    for i, exct in enumerate(self.list):
      del self.list[i][key]

  #Routine for getting energies as numpy array out of list of dict
  def list2energies():
    energies = [x['energy'] for x in self.list]

    return energies

  #Routine for default excitation
  def defaultExcitation():
    exct = {'name': 'S', 'fix': False, 'energy': 0.0, 'strength': None,
            'phase': 0.0, 'amplitudes_file1': [None]*3,
            'amplitudes_file2': [None]*3, 'amplitudes_file3': [None]*3}
    
    return exct

    
