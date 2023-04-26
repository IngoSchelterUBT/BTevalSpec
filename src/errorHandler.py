import sys

#Module for output of ERRORS and WARNINGS

def err(code,msg):
    print
    print(msg)
    sys.exit(code)

def warn(msg):
    print('WARNING: ' + msg)
