# Purpose

In general, Bayreuth spectrum evaluation script `BTevalSpec.py` takes the induced time-dependent dipole moment, e.g., from an electronic real-time TDDFT calculation, computes the corresponding spectrum, and evaluates the latter for excitation energies and oscillator strengths by fitting spectral lines with the proper shape.
If densities `n(r,w_k)` are given at energies `w_k` close to the actual excitation energies, BTevalSpec.py can also compute the true transition densities.

The script is able to evaluate dipole-moment data emerging from a so called boost excitation (a spectrally broad kick-like excitation) or a time-dependent electric dipole field.
The dipole-moment data may also be devided into different spatial regions that are, e.g., associated with different molecules to compute molecular contributions to supermolecular excitations.
Finally, the script accepts dipole-moment data from different calculations with different external-field polarizations to improve the evaluation.

# Reference

The underlying theory and algorithms are explained in

 - Ingo Schelter and Stephan Kümmel, "Accurate Evaluation of Real-Time Density Functional Theory Providing Access to Challenging Electron Dynamics", JCTC 14, 1910-1927 (2018), doi: 10.1021/acs.jctc.7b01013,
 - Ingo Schelter, Johannes M. Foerster, Rian Richter, Nils Schild, and Stephan Kümmel, unpublished (2025)

Please, cite the above references if you use BTevalSpec.py.

# Installation

## Requirements

BTevalSpec.py requires Python3 (tested for 3.10.5) with the following modules:
 - numpy
 - sys
 - os
 - scipy
 - concurrent.futures
 - matplotlib
 - getopt
 - re
 - ruamel.yaml
 - lmfit
 - numba
 - pyinstaller

## Directory

The following directories exist:
 - `src/` contains the source code (main script `eval.py`) and install script (`install.sh`), which can build a one-file executable using `pyinstaller`.
 - `testTemplates/` provides several test cases including input files. This directory is a template that can be copied for actual test directories.
 - `tools/` contains templates for auxiliary gnuplot and shell scripts, e.g., for visualizing data.

# Input file(s)

There are two kinds of input files: dipole moments and the laser profile (cf. the test directory for examples).

## Dipole moment

A dipole-moment file consist of a header and a data section.
The data section must provide four columns containing the equidistant time grid (col 1) and the x/y/z components of the induced(!) dipole moment (col 2-3), all in Rydberg atomic units.

```
#time                 |           x          |           y          |           z          |
 0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00  0.000000000000000E+00
 2.067068666828068E-01 -7.324170588265575E-17  8.349368755249521E-16  3.133584764152223E-10
 4.134137333656137E-01 -1.851177562732433E-16  1.160744093695872E-15  3.729813588188685E-09
 6.201206000484205E-01 -2.702214153295676E-16  4.022430641288225E-16  1.749995353858235E-08
 8.268274667312273E-01 -6.148812091645332E-16  2.360695496697065E-15  5.373312737625561E-08
 ...
```

The header contains meta information such as the type, strength, and polarization of the external perturbation.
Some quantities carry units, others don't.

```
# Column
# 1  2                    3      4                                          description
#
#Currently used keys
#!BT NELEC                       1.000000E+00,1.000000E+00                  Number of electrons (spin-dependent)    
#!BT BOOSTMASK                   no                                         Spatial-area in which a boost was applied ("no": no boost, "global": global boost, further masks are possible)
#!BT BOOSTENERGY          Ry     0.000000E+00                               Boost energy = Nelec*hbar^2/2m kboost^2
#!BT BOOSTKVEC                   0.000000E+00,0.000000E+00,0.000000E+00     Boost polarization (auto-normalized, boost strength is taken from BOOSTENERGY)
#!BT EXCITATION                  laser_mono                                 Type of excitation ("laser_mono" for an electric dipole field)
#!BT LASEREFIELD          E_ry   2.387186E-05                               Amplitude of the electric field in Ry a.u. so that the time-dependent electric field is LASERFIELD*LASERPOL(X,Y,Z)*laser_profile.dat
#!BT LASERPOLX                   0.000000E+00                               Laser polarization: x-component
#!BT LASERPOLY                   0.000000E+00                               .. y-component                
#!BT LASERPOLZ                   1.000000E+00                               .. z-component (auto-normalized)
#!BT LASEREND             t_ry   2.067069E+02                               End of the external field (in Ry time units) = begin of free propagation

#Further keys with potential usage in the future
#!BT VERSION                     v3.5.0                                     Version string (e.g., BTDFT)
#!BT DATE                        05-APR-2024_12:55:24_+0200                 Date 
#!BT SYSID                       [none]                                     System id string  
#!BT NATOMS                      2                                          Number of atoms
#!BT NTYPES                      1                                          Number of atom types
#!BT ATOMSPERTYPE                2                                          Number of atoms per atom type 
#!BT DX                   a0     7.000000E-01                               Grid spacing of the real-space grid
#!BT T_START              t_ry   0.000000E+00                               Global starting time
#!BT DT                   t_ry   2.067069E-01                               Time-step size
#!BT EXCITMASK                   global                                     Excitation mask if external field is applied with a spatial profile, e.g., in a certain area
#!BT NLASER                      1                                          Number of laser excitations (only one supported at the moment)
#!BT LASERFREQ            Ry     1.400013E-01                               Laser frequency
#!BT LASERINT             W/cm^2 1.000000E+07                               Laser intentity  
#!BT LASERSTART           t_ry   0.000000E+00                               Laser start time  
#!BT OBSERVABLE                  dipole                                     How the observable in this file is called 
#!BT STRIDE                      1                                          The time-step stride with which this observable was written from the program
#!BT MASK                        global1                                    The identifier of the spatial mask that defines the spatial area in which the data in this file were recorded
#!BT NAREAS                      1                                          The number of spatial areas
#!BT AREA                        1                                          The number of this area (starting at 1)
#!BT AREA_NAME                   global                                     The name of the spatial area 
#!BT SPIN                        tot                                        If this is a spin (up/dn) density of the total (tot) density
#!BT DIP0                 dip_ry 4.759665E-14,-6.683266E-14,1.188267E-12    The initial dipole moment (which was subtracted from the time-dependent dipole moment to get the induced dipole moment in the data section)
```

All header lines must begin with a "#!BT" (col 1) followed by a key (col 2).
If a key expects a dimensional value, column 3 contains a unit string.
Otherwise, it is empty.
Finally, column 4 contains the actual value.
All existing key-value pairs are described above (some are currently not and may never used).

## Laser profile

The laser-profile file together with the dipole-moment header keys "LASERFIELD" and "LASERPOL[XYZ]" defines the time-depend electric field and is only required for a electric-dipole-field excitation.
It must contain 1+NLASER columns with the time grid (col 1) and the time-profile of the electric dipole field in col (2:NLASER+1).
Currently, only a single laser is supportet.

# Run

`BTevalSpec.py` requires some custom modules that are included in the `src` directory.
Therefore, you can either call `BTevalSpec.py` directly from the src directory or create a symbolic link or build a one-file executable using `pyinstaller` and put into a location in your binary path.

```
./BTevalSpec.py -h
--------------------
Usage:

  ./BTevalSpec.py [-h |--help] [{-v |--verbose=}<dbg>] [{-f |--file=}<fname.yaml>] [<cmd-opt>] cmd arg

  General options
     -h, -v, and -f are referenced as general options (<gen-opt>) below
     -v (verbose) has several levels: 0 produces minimal output, 1 is usually useful, 2 is for debugging.
     -f <file>    specify a different input yaml file (default: eval.yaml)

  Help
    -h all   Prints all help messages
    -h gen   Prints general information and workflow
    -h yaml  Prints the eval.yaml-specific help message.
    -h cmd   Prints an one-line overview over all commands
    -h <cmd> Prints the <cmd>-specific help message.

  Reference:
    The algorithms are explained in
      I. Schelter, S. Kuemmel,         JCTC 14, 1910-1927 (2018)
      I. Schelter, S. Kuemmel, et al., unpub. (2025) (will be updated when the paper is published)

  Testing / Examples:
    Examples/Tests are provided in the 'testTemplates' directory. Copy this directory an run the included tests.
--------------------
```

## Overview

To run `BTevalSpec.py`
 1. Create a directory
 2. Put in or link the input files (dipole moments and laser profile if required)
 3. Either link `BTevalSpec.py` into the directory or call a pyinstalled-version of the script
 4. Create a new `eval.yaml` using the `new` command below (see below)
 5. Edit the newly created `eval.yaml` file (e.g., enter the correct file names and the fit range)
 6. [optional] Do a separate Fourier transform or Pade spectrum to get a non-evaluated spectrum
 7. Run the script repeatedly adding more excitations until you are satisfied ([optional] backup the current eval.yaml after each step so you can fall back to a previous stage)
 8. Evaluate the excitations found using the implemented significance measures
 9. [optional] Decouple Fourier-transformed densities into transition densities

```
./BTevalSpec.py -h gen

--------------------
  Workflow
    1) Copy/Link all required files into the evaluation directory (dipole moments, external-field profile, n(r,w_k), BTevalSpec.py)

    2) Create a new input-file template
       ./BTevalSpec.py new

    3) Edit the newly created eval.yaml
       Add dipole-moment files and external-field profile (if required)
       Adjust the fit range (if you don't know one yet, see below)

    Note: After each call to BTevalSpec.py in the following, the eval.yaml file is updated and automatically used for the next call to BTevalSpec.py
          This way, one can manually run the script step-by-step and intervene whenever it is necessary.
          It is also possible to manually adjust the eval.yaml file.

    4) Add & fit new excitations one after the other (repeat the following steps)
       0a) [optional] If you don't know how the spectrum looks like, first compute the Fourier- or Pade spectra initially and inspect them
          ./BTevalSpec.py {ft|pade}
       0b) [optional] If you want to make an initial guess without fitting, call
          ./BTevalSpec.py [--guess={ft|pade}] guess
       a) Run one iteration with the fit command:
          ./BTevalSpec.py --skip --niter=1 [--nsig=<val>] fit
             --skip skips the first fit of the input excitations
             --niter=1 makes one iteration consisting of (i) guess a set of new excitations (ii) fit new excitations (iii) fit all excitations together (iv) compute error and significance measures
             --nsig=<val> Adjust the threshold that determines if the script accepts a peak as an excitation. A larger value accepts less peaks.
          This automatically evaluates the current difference between computed spectrum and the fit, finds new excitations at peak positions, and fits the excitations.
          Note: This command automatically computes the Fourier and Pade spectra if not already done. Later calls will just read the latter.
       b) [optional] backup the updated eval.yaml so you can return to this stage later if necessary
       c) check significance measures (in eval.yaml or comprehensively in excit_1.dat) and decide to
          - stop the fit
          - adjust the excitations manually
          - do another iteration -> 4a)
          - plot spectra/fit/...
             ./BTevalSpec.py plot
          Note: A low significance can indicate that an excitation is erroneous. However, it can also show that there are still excitations missing in the vicinity.

    5) Evaluate transition densities (requires the DENSFT section in eval.yaml with as many n(r,w_k) as there are excitations)
       ./BTevalSpec.py decouple
--------------------
```

## eval.yaml

```
./BTevalSpec.py -h yaml
--------------------
 eval.yaml file format:

    Note: Most fields in the yaml file are filled automatically.
          After creating a new eval.yaml, you only need to specify/adjust
           - The file paths (dipole moments and laser profile if required)
           - The fit range (if you run a fit)
           - wmax (PADE) if you want a Pade spectrum

    Note: After each call to BTevalSpec.py in the following, the eval.yaml file is updated and automatically used for the next call to BTevalSpec.py
          This way, one can manually run the script step-by-step and intervene whenever it is necessary.
          It is also possible to manually adjust the eval.yaml file.

DESCRIPTION:                        # Description of evaluation
DIPOLE:                             # List of dipole moment files
                                    #   Outer list: dipole moments from different calculations (different external-field polarizations)
                                    #   Inner list: dipole moments from different spatial areas from the same calculation
- [dipole_calc1_area1.dat, dipole_calc1_area2.dat]
- [dipole_calc2_area1.dat, ...]
  #...
DENSFT:                             # Specify Fourier-transformed densities n(r,w_k) and w_k for transition-density calculation
  densft:                           #   List of n(r,w_k) (as many as excitations at w_k close to the excitation energies)
  - densft01.compact                #     Either as complex-valued BTcompact files
  - densft02.compact
  #OR
  - [densft01r.cube,densft01i.cube] #     or lists of [real,imag] Gaussian cube files
  - [densft02r.cube,densft02i.cube]
  #OR
  - densft01i.cube                  #     or imag-valued Gaussian cube files with the --imag to the decouple command option
  - densft02i.cube
  #...
  densen:                           #   the list of w_k (same order)
  - 0.131990
  - 0.134330
  #...
  jcalc: 0                          #   the calculation which created the n(r,w_k) (in case dipole-moment data came from different calculations). Starts with '0' (default)
OPT:                                # Options & values for
  FT:                               #   Fourier transform
    calc: false                     #     true: compute, false: read from files
    minpw: 17                       #     zero padding up to (at least) 2^<minpw>
    window: 0.0                     #     Kaiser-Bessel Window parameter (leave 0 for fit)
    smooth: 0.0                     #     Damping rate (leave 0 for fit)
    rmDC: false                     #     Subtract omega=0 component from the dipole moment (leave false)
  Pade:                             #   Pade approximation
    calc: false                     #     true: compute, false: read from files
    wmax: 0.5                       #     Upper energy bound [Ry] for Pade approximation (increase to approximate a larger part of the spectrum)
    dw: 1e-05                       #     Energy sampling (distance between omega sampling points)
    smooth: 0.0022278747967549423   #     Artificial decay rate used for Pade approximation (decrease to get sharper lines)
    thin: 0                         #     Use every 2^<thin>'s time step to evaluate the Pade spectrum (artificially increases the time-step size for comput. efficiency; leave zero)
  Fit:                              #   Fit
    range:                          #     Range of spectrum which should be fitted in Ry
    - 0.10
    - 0.28
    imagonly: false                 #     [in, optional] Set true if only the imaginary part shall be used (automatically true for a boost)
    fiterr: 0.01029415706259413     #     [out         ] Comprehensive fit error
    fitphase: false                 #     [in, optional] If true, the phase is used as a fit parameter (otherwise, it determined from the external-field profile)
EXT:                                # External-field
  profile: laser_profile.dat        #   [in, if laser excitation] Laser-profile file
  invertPhase: false                #   [in, optional           ] If true, multiplies the laser profile with '-1' (e.g., for BTDFT until v3.6.0)
SPEC:                               # Spectral information (fit)
- name: S1                          #   Excitation name (defaults to S<number>) but can be changed
  fix: false                        #   True, if the excitation is fixed (can be set manually for the next run)
  energy: 0.13198769111844885       #   Current value of the excitation's energy
  phase: 0.0                        #   Current value of the excitation's phase
  dipoles:                          #   Area-contributions to the transition dipole (consistent over calculations)
  - - 0.0014764204066536474
    - 0.3074363126526123
    - 0.3473751252687603
  ampl:                             #   Calculation- and Area-contributions to the vector amplitudes (different for different calculations)
  - - - -1.977874864015696e-06
      - -0.0004118546129144092
      - -0.000465358325824443
  dipole:                           #   Global transition dipole
  - 0.0014764204066536474
  - 0.3074363126526123
  - 0.3473751252687603
  strength: 0.004733710908363272    #   Oscillator strength
                                    #   Significance measures:
  signifFit: 0.0                    #   Only computed with --signif option to the 'fit' command
  signifErr: 0.9334766695559699
  signifAng: 0.9037784539834844
  signifExc: 0.0                    #   Only computed with --signif option to the 'fit' command
  signifRng: 0.9999997972734961
  signifPha: 1.0
    #...                            #   There are further quantities such as error estimates etc. for the single quantities
- name: S2                          # Further excitations
  #...
--------------------
```

## Commands

### Overview

```
./BTevalSpec.py -h cmd
--------------------
 Command overview:

  yaml     - [only with -h option] Print eval.yaml-specific help message
  cmd      - [only with -h option] Print eval.yaml-specific help message
  new      - create a new configuration file
  ft       - compute Fourier spectrum
  pade     - compute Pade    spectrum
  guess    - make an initial guess for the fit
  fit      - fit and add new excitations
  plot     - Plot certain spectra and fits
  rm       - remove an excitation
  add      - add    an excitation
  fix      - fix    an excitation (don't change its energy/dipole parameters during the fit)
  release  - opposite of fix
  reset    - reset the excitation's energy-interval that restricts the fit
  decouple - compute transition densities from Fourier-transformed densities (after the fit is complete)

    Note: After each call to BTevalSpec.py in the following, the eval.yaml file is updated and automatically used for the next call to BTevalSpec.py
          This way, one can manually run the script step-by-step and intervene whenever it is necessary.
          It is also possible to manually adjust the eval.yaml file.
```

### New

```
./BTevalSpec.py -h new
--------------------
 Create a new configuration file

  ./BTevalSpec.py [<gen-opt>] new

    Generates a new yaml file (default: eval.yaml).
    After creation, the user must add dipole-moment files and potentially a laser-profile file
    Moreover, the fit range [Ry] can be adjusted.

  Example (eval.yaml):

...
DIPOLE:                             # List of dipole moment files. Different calculations are only allowed to differ in the boost/laser polarization
  - [dipole_calc01_area01.dat, dipole_calc01_area02.dat]
  - [dipole_calc02_area01.dat, dipole_calc02_area02.dat]
EXT:
  profile: laser_profile.dat        # Profile of excitation
  invertPhase: false                # Before BTDFT v3.6.0, the laser profile missed a factor '-1'. InvertPhase==True compensates this error.
OPT:
  ...
  Fit:
    range:                          # Range of spectrum which should be fitted in Ry
    - 0.10                          # If you don't know a proper fit range yet, you separately compute the FT/Pade spectrum and inspect it visally.
    - 0.40
...
--------------------
```

### Fourier Transform

```
./BTevalSpec.py -h ft
--------------------
 Fourier transform

  ./BTevalSpec.py [<gen-opt>] [--minpw=<pw>] [--smooth=<smooth>] [--window=<window>] [--no-rmDC] ft

    Fourier transforms the dipole moment files, writes the transformations, and updates eval.yaml. 

  Note: If not disabled in the eval.yaml file, the Fourier transform is automatically computed if you use the 'fit' command.
        This command is useful, if you only want the Fourier spectrum in the first place without any fitting.

  Note: If the Fourier transform was computed previously, it is not computed again but read from the previously generated files
        which is automatically enforced by setting 'FT: {calc: false}' in the eval.yaml file.
        If, for some reason (e.g. you deleted the FT files), you want to compute the Fourier transform again, you have to set this option to 'true' again in eval.yaml.
--------------------
```

### Pade approximation

```
./BTevalSpec.py -h pade
--------------------
 Pade approx

  ./BTevalSpec.py [<gen-opt>] [--wmax=<wmax>] [--dw=<dw>] [--smooth=<smooth>] [--thin=<thin>] pade

    Pade approximates the dipole moment files, writes the transformations, and updates eval.yaml.

  Note: If not disabled in the eval.yaml file, the Pade spectrum is automatically computed if you use the 'fit' command.
        This command is useful, if you only want the Pade spectrum in the first place without any fitting.

  Note: If the Pade spectrum was computed previously, it is not computed again but read from the previously generated files
        which is automatically enforced by setting 'Pade: {calc: false}' in the eval.yaml file.
        If, for some reason (e.g. you deleted the FT files), you want to compute the Fourier transform again, you have to set this option to 'true' again in eval.yaml.
--------------------
```

### New guess

```
./BTevalSpec.py -h guess
--------------------
 Make a new guess based on the Fourier (default) or Pade spectrum

  ./BTevalSpec.py [<gen-opt>] [--guess=pade [<pade opt>] [--thres=<thres>]| --guess=ft [<ft opt>] [--nsig=<nsig>]] [--range=<lb,rb>] guess

    Creates a new guess, sets the plot range, and updates eval.yaml.
    The guess is based on the Pade (--guess=pade) or Fourier (--guess=ft, default) spectrum
    with the given <pade opt> or <ft opt> (if the latter were not computed previously).
--------------------
```

### Remove excitation

```
./BTevalSpec.py -h rm
--------------------
 Remove excitation without fit

  ./BTevalSpec.py [<gen-opt>] rm <exlist>

    Remove <exlist> from the list of excitations, where exlist is a comma-separated list of excitation labels, e.g., S1,S2,S3
    Instead, you can also remove teh excitation manually from eval.yaml.
--------------------
```

### Add excitation

```
./BTevalSpec.py -h add
--------------------
 Add excitation without fit

  ./BTevalSpec.py [<gen-opt>] [--nsig=<nsig> | --nex=<nex> | --energy=<list-of-energies>] [--nofix] add

    Add new excitations either at given energies (optional) or automatically using a sigma-threshold (default: 2.0), i.e.,
    add excitations that are larger than the mean value of maxima of ft-fit plus <nsig> times the standard deviation of maxima heights
    Do no fit but guess phase and dipoles moments.
    Before adding new excitations, the existing ones are fixed. Switch that off with the --nofix option
--------------------
```

### Reset excitation energy interval

```
./BTevalSpec.py -h reset
--------------------
 Reset energy range to the standard pi/T interval around the current excitations' energy values

  ./BTevalSpec.py [<gen-opt>] reset

  Note: When an excitation is added by the code, its energy is restricted to a symmetric interval around the initial energy value to stabilize the fit.
        If you note that a rng significance value of an excitation is not close to 1, it may be useful to reset the energy ranges and, thus, give more freedom to the energy fit
--------------------
```

### Fit

```
./BTevalSpec.py -h fit
--------------------
 Fit

  ./BTevalSpec.py [<gen-opt>] [<ft opt>] [<guess opt>] [--reset] [--skip] [--single] [--signif] [--range=<lb,rb>] [--crit=<error criterion> | --nadd=<nadd>] [--nsig=<nsig>] [--niter=<niter>] [--imag] [--fitphase] fit

    Fit current excitations (if not --skip), add and fit <nadd> excitations one after the other, and update eval.yaml.
    --skip: Skip the first collective fit of existing excitations (useful if eval.yaml was changed manually)
    --imag: only use imaginary part for fitting; automatically true for boost excitation.
    <nadd> defaults to 0
    --single: Instead of fitting all excitations at once, fit them one after the other from high to low strength 
    Does a prior Fourier transform if not already done.
    Does a prior Guess if not already done; requires a range option.
    --signif: Compute the significance-values that require fitting (are costly)
    --fitphase: Use phase as fit parameter
    --crit: Convergence criterion (only relevant if nadd==0 or not present): Add excitations until the fit error drops below this criterion.
    --niter: Maximum number of add-excitation iterations (in each iteration, the number of added excitations is determined by the <nsig> threshold)
    --reset: Reset energy range

  Note: The simplest approach to fit a spectrum is to do
        1) one iteration at a time using
           `./BTevalSpec.py --skip [--nsig=<nsig>] --niter=1 fit`
        2) then evaluate the added excitations (e.g. by inspecting the significance measures)
        3) backup the resulting eval.yaml
        4) go to 1) (use --skip to suppress the initial fitting of the input excitations, use the --nsig option to change the threshold for adding new excitation)

 Note: If any significance measure is not close to 1, this can (!) indicate that an excitation is erroneous. Try to remove the the excitation using
           `./BTevalSpec.py rm <excitation identifier>`
       and call
           `./BTevalSpec.py fit`
       to just fit the remaining excitations again
       If the rng significance of an excitation is not close to 1, resetting the energy range using --reset gives the fitting more freedom, which helps sometimes.

 Note: At the end of this command, the script prints what it would do next (i.e., which excitations it would add next).
       Don't be confused, it looks as if the script had actually done another fit iteration, however, it is just an information.

 Note: Several output files are generated that, e.g., contain the spectra, the fit, the fit-objective function (=kind of f-resolved fit error), a table of excitations, etc.
       Also, during the fit, the scripts opens a window with the current fit-objective and nsig threshold before fitting in each iteration to show what is actually done.
--------------------
```

### Plotting
       
```
./BTevalSpec.py -h plot
--------------------
 Plot

  ./BTevalSpec.py [<gen-opt>] [--exclude=<listOfExcitations>] plot <listOfMeasures>

    Plots <listOfMeasure> in {pade, ft, fit, err, spectrum} with excitations
    --exclude=<listOfExcitations>: Excitations that are excluded from the (fitted) data
--------------------
```

### Fix excitations

```
./BTevalSpec.py -h fix
--------------------
 Fix

```
  ./BTevalSpec.py [<gen-opt>] [--invert] fix <listOfExcitations>

    Fix the given list of excitations or energy range (those are excluded from the fit)
    --invert: Do the same but inverted
    You can do this manually in eval.yaml as well
--------------------
```

### Release excitations

```
./BTevalSpec.py -h release
--------------------
 Release

  ./BTevalSpec.py [<gen-opt>] [--invert] release <listOfExcitations>

  Release (un-fix) all excitations, cf. 'fix' command
--------------------
```

### Compute transition densities

```
./BTevalSpec.py -h decouple
--------------------
 Decouple

  ./BTevalSpec.py [<gen-opt>] [--jcalc=<calc-idx>] decouple

  Decouple given Fourier-transformed densities n(r,omega) at given omega to get the proper transition densities.
  Requires Fourier-transform of the density at as many energies as there are excitations (ideally the excitation energies themselves) and a calculation identifier (calculation index).
  jcalc determines the calculations index from which the density stems (default: 0).

DENSFT:
  densft:                           #   List of n(r,w_k) (as many as excitations at w_k close to the excitation energies)
  - densft01.compact                #     Either as complex-valued BTcompact files (requires BTDFT-modules)
  - densft02.compact
  #OR
  - [densft01r.cube,densft01i.cube] #     or lists of [real,imag] Gaussian cube files
  - [densft02r.cube,densft02i.cube]
  #OR
  - densft01i.cube                  #     or imag-valued Gaussian cube files with the --imag to the decouple command option
  - densft02i.cube
  #...
  densen:                           #Actual omega [Ry] at which the FT densities above are given
  - 0.131990
  - 0.134330
    ...
  jcalc: 0
--------------------
```

## Examples/Testing

Copy the `testTemplates/` directory into a directory, e.g., `test~/` (the tilde makes git ignore this directory)
```
cp -rP testTemplates test~
```
for actual testing.
The `-P` option leaves symbolic links as such.

Here is an overview of the single test cases:

[Todo]

Confer the `readme.md` files in the single test directories for instructions.

