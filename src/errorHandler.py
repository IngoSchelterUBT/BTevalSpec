import sys

#Module for output of ERRORS and WARNINGS

#------------------------------------------------------------------------------#
# Error handler
# Input: code (Error code), msg (Error Message)
#------------------------------------------------------------------------------#
def err(code,msg,cmd=""):
    print("ERROR: " + msg)
    usage(cmd=cmd)
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
def usage(cmd=""):
    print("--------------------")
    print("")
    print("Usage:")
    print("")
    print("  ./eval.py [-h |--help] [{-v |--verbose=}<dbg>] [{-f |--file=}<fname.yaml>] [<cmd-opt>] cmd arg")
    print("")
    print("  -h, -v, and -f are referenced as general options (<gen-opt>) below")
    if cmd in ["","new"]:
        print("")
        print(" Create a new configuration file")
        print("  ./eval.py [<gen-opt>] new")
        print("    Generates a new yaml file (default: eval.yaml) file.")
        print("    After creation, the user must add files containing the dipole moments and extern profile")
    if cmd in ["","ft"]:
        print("")
        print(" Fourier transform")
        print("  ./eval.py [<gen-opt>] [--minpw=<pw>] [--smooth=<smooth>] [--window=<window>] [--no-rmDC] ft")
        print("    Fourier transforms the dipole moment files, writes the transformations, and updates eval.yaml. ")
    if cmd in ["","pade"]:
        print("")
        print(" Pade approx")
        print("  ./eval.py [<gen-opt>] [--wmax=<wmax>] [--dw=<dw>] [--smooth=<smooth>] [--thin=<thin>] pade")
        print("    Pade approximates the dipole moment files, writes the transformations, and updates eval.yaml.")
    if cmd in ["","guess"]:
        print("")
        print(" Make a new guess based on the pade approximation")
        print("  ./eval.py [<gen-opt>] [--guess=pade [<pade opt>] [--thres=<thres>]| --guess=ft [<ft opt>] [--nsig=<nsig>]] [--range=<lb,rb>] guess")
        print("    Creates a new guess, sets the plot range, and updates eval.yaml.")
        print("    Does a prior Pade approximation with the given <pade opt> options if not already done.")
    if cmd in ["","rm"]:
        print("")
        print(" Remove excitation without fit")
        print("  ./eval.py [<gen-opt>] rm <exlist>")
        print("    Remove <exlist> from the list of excitations, where exlist is a comma-separated list of excitation labels, e.g., S1,S2,S3")
    if cmd in ["","add"]:
        print("")
        print(" Add excitation without fit")
        print("  ./eval.py [<gen-opt>] [--nsig=<nsig> | --nex=<nex> | --energy=<list-of-energies>] [--nofix] add")
        print("    Add new excitations either at given energies (optional) or automatically using a sigma-threshold (default: 2.0), i.e.,")
        print("    add excitations that are larger than the mean value of maxima of ft-fit plus <nsig> times the standard deviation of maxima heights")
        print("    Do no fit but guess phase and dipoles moments.")
        print("    Before adding new excitations, the existing ones are fixed. Switch that off with --nofix")
    if cmd in ["","reset"]:
        print("")
        print(" Reset energy range to the standard pi/T interval around the current excitations' energy values")
        print("  ./eval.py [<gen-opt>] reset")
    if cmd in ["","fit"]:
        print("")
        print(" Fit")
        print("  ./eval.py [<gen-opt>] [<ft opt>] [<guess opt>] [--skip] [--single] [--signif] [--range=<lb,rb>] [--crit=<error criterion> | --nadd=<nadd>] [--nsig=<nsig>] [--niter=<niter>] [--imag] [--fitphase] fit")
        print("    Fit current excitations (if not --skip), add and fit <nadd> excitations one after the other, and update eval.yaml.")
        print("    --imag: only use imaginary part for fitting; automatically true for boost excitation. (not supported yet)")
        print("    <nadd> defaults to 0")
        print("    --single: Instead of fitting all excitations at once, fit them one after the other from high to low strength ")
        print("    Does a prior Fourier transform if not already done.")
        print("    Does a prior Guess if not already done; requires a range option.")
        print("    --signif: Compute significances that require fitting")
        print("    --fitphase: Use phase as fit parameter")
        print("    --crit: Convergence criterion (only relevant if nadd==0 or not present)")
        print("    --niter: Maximum number of add-excitation iterations")
    if cmd in ["","plot"]:
        print("")
        print(" Plot")
        print("  ./eval.py [<gen-opt>] [--exclude=<listOfExcitations>] plot <listOfMeasures>")
        print("    Plots <listOfMeasure> in {pade, ft, fit, err, spectrum} with excitations")
        print("    --exclude=<listOfExcitations>: Excitations that are excluded from the (fitted) data")
    if cmd in ["","fix"]:
        print("")
        print(" Fix")
        print("  ./eval.py [<gen-opt>] [--invert] fix <listOfExcitations>")
        print("    Fix the given list of excitations or energy range")
        print("    --invert: Do the same but inverted")
    if cmd in ["","release"]:
        print("")
        print(" Release")
        print("  ./eval.py [<gen-opt>] [--invert] release <listOfExcitations>")
        print("  Release (un-fix) all excitations")
    if cmd in ["","decouple"]:
        print("")
        print(" Decouple")
        print("  ./eval.py [<gen-opt>] [--jcalc=<calc-idx>] [--ftype=<file type>] decouple")
        print("  Decouple given transition densities.")
        print("  Requires Fourier-transfor of the density at as many energies as there are excitations (ideally the excitation energies themselves) and a calculation identifier (calculation index)")
        print("  ftype can be 'compact' or 'cube' (default: 'cube')")
        print("  jcalc determines the calculations index from which the density stems (default: 0)")
    print("")
    print("--------------------")
