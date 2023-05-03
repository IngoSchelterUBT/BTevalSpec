#==============================================================================#
# Utilities
#==============================================================================#
import errorHandler as err

#--------------------------------------------------------------------------#
# Convert index to direction
#--------------------------------------------------------------------------#
def getDir(i):
  if i == 0:
    return 'x'
  elif i == 1:
    return 'y'
  elif i == 2:
    return 'z'
  else:
    err.err(1,'There can not be more directions than x, y and z!')

#--------------------------------------------------------------------------#
# Routine checks if input can be converted into float
#--------------------------------------------------------------------------#
def isFloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False
