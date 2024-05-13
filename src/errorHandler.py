import sys

#Module for output of ERRORS and WARNINGS

#------------------------------------------------------------------------------#
# Error handler
# Input: code (Error code), msg (Error Message)
#------------------------------------------------------------------------------#
def err(code,msg):
    print("ERROR: " + msg)
    usage()
    sys.exit(code)

#------------------------------------------------------------------------------#
# Warning
# Input: msg (Error Message)
#------------------------------------------------------------------------------#
def warn(msg):
    print('WARNING: ' + msg)

#------------------------------------------------------------------------------#
# Print usage information
#------------------------------------------------------------------------------#
def usage():
    print("--------------------")
    print("")
    print("Usage:")
    print("")
    print("  ./eval.py [--debug=<dbg>] [<cmd-opt>] cmd arg")
    print("")
    print(" Create a new eval.yaml file")
    print("  ./eval.py new [name.yaml]")
    print("    Generates a new yaml file (default: eval.yaml) file.")
    print("    After creation, the user must add files containing the dipole moments and extern profile")
    print("")
    print(" Fourier transform")
    print("  ./eval.py [--minpw=<pw>] [--smooth=<smooth>] [--window=<window>] ft")
    print("    Fourier transforms the dipole moment files, writes the transformations, and updates eval.yaml. ")
    print("")
    print(" Pade approx")
    print("  ./eval.py [--wmax=<wmax>] [--dw=<dw>] [--smooth=<smooth>] [--thin=<thin>] pade")
    print("    Pade approximates the dipole moment files, writes the transformations, and updates eval.yaml.")
    print("")
    print(" Make a new guess based on the pade approximation")
    print("  ./eval.py [<pade opt>] [--thres=<thres>] guess <lb,rb>")
    print("    Creates a new guess, sets the plot range, and updates eval.yaml.")
    print("    Does a prior Pade approximation with the given <pade opt> options if not already done.")
    print("")
    print(" Remove excitation without fit")
    print("  ./eval.py rm <exlist>")
    print("    Remove <exlist> from the list of excitations, where exlist is a comma-separated list of excitation labels, e.g., S1,S2,S3")
    print("")
    print(" Add excitation without fit")
    print("  ./eval.py [--label=<label>] [--nsig=<nsig>] add [<energy>]")
    print("    Add new excitation with given label (optional) at given energy (optional).")
    print("    Do no fit but guess phase and dipoles moments.")
    print("    Add excitations that are larger than the mean value of maxima of ft-fit plus <nsig> times the standard deviation of maxima heights")
    print("")
    print(" Reset energy range to the standard pi/T interval around the current excitations' energy values")
    print("  ./eval.py reset")
    print("")
    print(" Fit")
    print("  ./eval.py [<ft opt>] [<guess opt>] [--skip] [--single] [--range=<lb,rb>] [--imag] fit [<nadd>] [<listOfExcitations>]")
    print("    Fit current excitations (if not --skip), add and fit <nadd> excitations one after the other, and update eval.yaml.")
    print("    --imag: only use imaginary part for fitting; automatically true for boost excitation.")
    print("    <nadd> defaults to +0")
    print("    <exs>    Comma-separated list of excitations to fit")
    print("    --single: Instead of fitting all excitations at once, fit them one after the other from high to low strength ")
    print("    Does a prior Fourier transform if not already done.")
    print("    Does a prior Guess if not already done; requires a range option.")
    print("")
    print(" Error")
    print("  ./eval.py error")
    print("    Compute error value")
    print("")
    print(" Significance")
    print("  ./eval.py significance [<listOfExcitations>]")
    print("    Compute significances of all lines or the excitations given in the comma-separated list of <listOfExcitations>")
    print("")
    print(" Plot")
    print("  ./eval.py [--exclude=<listOfExcitations>] plot <listOfMeasures>")
    print("    Plots <listOfMeasure> in {pade, ft, fit, err, spectrum} with excitations")
    print("    --exclude=<listOfExcitations>: Excitations that are excluded from the (fitted) data")
    print("")
    print(" Fix")
    print("  ./eval.py [--invert] fix {<listOfExcitations>|<energyRange>}")
    print("    Fix the given list of excitations or energy range")
    print("    --invert: Do the same but inverted")
    print("")
    print(" Release")
    print("  ./eval.py release")
    print("  Un-fix all excitations")
    print("")
    print("--------------------")
