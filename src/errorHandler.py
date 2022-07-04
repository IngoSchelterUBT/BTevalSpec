import sys

#Module for output of ERRORS and WARNINGS

def err(code,msg):
  print
  print(msg)
  sys.exit(code)

def warning(code,msg):
  print('WARNING: ' + msg)
