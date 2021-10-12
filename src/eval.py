#!/usr/bin/env python
#------------------------------------------------------------------------------#
# author: Rian Richter
# Fit of spectrum
import numpy as np
import re
import sys
import os.path
import getopt
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#only for tests
from pprint import pprint


#Import Modules
import Data


#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():
  data = Data.Data()
  #data.createConfig()
  #data.readConfig()
  #pprint(vars(data))

#------------------------------------------------------------------------------#
# Call main
#------------------------------------------------------------------------------#
if __name__=="__main__":
    main()
