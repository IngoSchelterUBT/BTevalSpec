#Routines which are needed in a lot of modules

#import own modules
import errorHandler as err

#Get direction
def getDir(i):
  if i == 0:
    return 'x'
  elif i == 1:
    return 'y'
  elif i == 2:
    return 'z'
  else:
    err.err(1,'There can not be more directions than x, y and z!')
