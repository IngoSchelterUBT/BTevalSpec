# BTevalSpec.py - Documentation

## Authors & Correspondence

Authors: Ingo W. Schelter (<ingo.schelter@uni-bayreuth.de>) and Rian R. Richter

## Cite

The underlying theory and algorithms are explained in

- [Schelter2018] Ingo Schelter and Stephan Kümmel, "Accurate Evaluation of Real-Time Density Functional Theory Providing Access to Challenging Electron Dynamics", JCTC 14, 1910-1927 (2018), doi: 10.1021/acs.jctc.7b01013,
- [Schelter2025] Ingo Schelter, Johannes M. Foerster, Rian Richter, Nils Schild, and Stephan Kümmel, unpublished (2025)

Please, cite the latest reference if you use BTevalSpec.py.

## Purpose

In general, Bayreuth spectrum evaluation script `BTevalSpec.py` takes the induced time-dependent dipole moment, e.g., from an electronic real-time TDDFT calculation, computes the corresponding spectrum, and evaluates the latter for excitation energies and oscillator strengths by fitting spectral lines with the proper shape.
If densities `n(r,w_k)` are given at energies `w_k` close to the actual excitation energies, BTevalSpec.py can also compute the true transition densities.

The script is able to evaluate dipole-moment data emerging from a so called boost excitation (a spectrally broad kick-like excitation) or a time-dependent electric dipole field.
The dipole-moment data may also be devided into different spatial regions that are, e.g., associated with different molecules to compute molecular contributions to supermolecular excitations.
Finally, the script accepts dipole-moment data from different calculations with different external-field polarizations to improve the evaluation.

## Installation

### Requirements

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
- [optional] pyinstaller
- [optional] mpi4py (for BTDFT support)

### Directory

The following directories exist:

- `src/` contains the source code (main script `eval.py`) and install script (`install.sh`), which can build a one-file executable using `pyinstaller`.
- `testTemplates/` provides several test cases including input files. This directory is a template that can be copied for actual test directories.
- `tools/` contains templates for auxiliary gnuplot and shell scripts, e.g., for visualizing data.

## Input file(s)

There are two kinds of input files: dipole moments and the laser profile (cf. the test directory for examples).

### Dipole moment

A dipole-moment file consist of a header and a data section.
The data section must provide four columns containing the equidistant time grid (col 1) and the x/y/z components of the induced(!) dipole moment (col 2-3), all in Rydberg atomic units.

```text
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

```text
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

All header lines must begin with a `#!BT` (col `1`) followed by a key (col `2`).
If a key expects a dimensional value, column `3` contains a unit string.
Otherwise, it is empty.
Finally, column `4` contains the actual value.
All existing key-value pairs are described above (some are currently not and may never used).

### Laser profile

The laser-profile file together with the dipole-moment header keys `LASERFIELD` and `LASERPOL[XYZ]` defines the time-depend electric field and is only required for a electric-dipole-field excitation.
It must contain `1+NLASER` columns with the time grid (col `1`) and the time-profile of the electric dipole field in col `2:NLASER+1`.
Currently, only a single laser is supportet.

## Add-Line Objective & Error Suppression

### Description

The script frequently plots a kind of spectrum which is called the "add-line objective" in the following and the examples/tests.
This function (cf. [Schelter2025]) is a comprehensive, energy-resolved measure for how much the current fit deviates from the actual data (i.e., the Fourier spectra from all provided dipole-moment files from different calculations, spatial areas, Cartesian components, and complex components).

The add-line objective is used to find new excitations and can be inspected during the fit procedure by the user to evaluate the quality of the fit.

### Line-Shape error & Error Compensation

Often, an excitation does not match its analytical line shape perfectly, which then leads to the "line-shape error".
This line-shape error remains visible in the add-line objective after fitting this excitation.
If there is a small excitation next to the large, already fitted excitation, it is superimposed by the latter's line-shape error.
To allow finding this small excitation "behind" the line-shape errors of larger excitations, `BTevalSpec.py` provides an error-suppression mechanism (cf. [Schelter2025]).
If the error suppression is active, there are two add-line objectives plotted: One with error suppression and one without as a reference.

Here is an example from the `Na2DA` test including the new-line guess (crosses) based on the error-suppressed add-line objective (dark solid line):

![Na2DA_addLineObj1_wref](testTemplates/Na2DA/readmeFiles/addLineObj1_wref.png   "Add-line objective (for fit iteration 2, wref=1) from the Na2DA test")

Remember that the add-line objective without error suppression (light gray line) shows the remaining part of the spectrum that is missing in the fit.
The error suppression in this case removes a large part of the add-line objective at abount `0.15 Ry`.
The difference between the line is the line-shape error, that is removed by the error suppression, i.e., the error from fitting previous lines at about `0.15 Ry`.
The previously fitted lines at that energy were more than 100 times larger than the remaining ones that are visible in the error-suppressed add-line objective.

## Run

`BTevalSpec.py` requires some custom modules that are included in the `src` directory.
Therefore, you can either call `BTevalSpec.py` directly from the src directory or create a symbolic link or build a one-file executable using `pyinstaller` and put into a location in your binary path.

```text
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

### Workflow

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

```text
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

### eval.yaml

```text
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

### Commands

#### Command Overview

```text
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

#### New

```text
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

#### Fourier transform

```text
./BTevalSpec.py -h ft

--------------------
Fourier transform

  ./BTevalSpec.py [<gen-opt>] [--minpw=<pw>] [--smooth=<smooth>] [--window=<win>] [--rmDC] ft

    Fourier-transforms the induced dipole moment, writes them into files, and updates eval.yaml. 

    <gen-opt>         General options
    --minpw=<pw>      Pad the induced dipole moment with zeros up to 2^<pw> to increase the frequency-space sampling rate (without changing the FT itself) (default: 17)
    --smooth=<smooth> Multiply the induced dipole moment with a decaying exponential function e^(-<smooth>*t) for smoother lines (default: 0)
    --window=<win>    Applies a Kaiser-Bessel window with parameter <win> to the induced dipole moment (default: 0)
    --rmDC            Remove the DC (omega=0) component of the induced dipole moment d(t) by subtracting int d(t) dt. This is usually not recommended.

  Attention: A smoothed or windowed signal cannot be used for subsequent fitting but is only intended for visual inspection.

  Note: The Fourier transform is automatically computed if you use the 'fit' command.
        The 'ft' command is useful, if you only want the Fourier spectrum in the first place without any fitting.

  Note: If the Fourier transform was computed previously, it is not computed again but read from the previously generated files
        This is automatically enforced by setting 'FT: {calc: false}' in the eval.yaml file after computing the Fourier transform once.
        If, for some reason (e.g. you deleted the FT files), you want to compute the Fourier transform again, you have to set this option to 'true' again in eval.yaml of use the 'ft' command manually.

--------------------
```

#### Pade approximation

```text
./BTevalSpec.py -h pade

--------------------
 Pade approx

  ./BTevalSpec.py [<gen-opt>] [--wmax=<wmax>] [--dw=<dw>] [--smooth=<smooth>] [--thin=<thin>] pade

    Computes the Pade spectrum of the induced dipole moment, writes them into files, and updates eval.yaml.

    <gen-opt>         General options
    --wmax=<wmax>     Upper omega boundary (in Ry units) (default: 0.5)
    --smooth=<smooth> Multiply the induced dipole moment with a decaying exponential function e^(-<smooth>*t) for smoother lines. This is alway applied for computing the Pade spectrum.
                      Smaller <smooth> values (=decay rates) result in sharper lines (but still smooth in contrast to the Fourier spectrum). Too small values, however, can produce artifacts. (default: ln(100)/<propagation time>)
    --dw=<dw>         Frequency-space sampling step size (default: 1E-5)
    --thin=<thin>     Thins out the dipole-moment samples using only every <thin>'s value to lower the computation burden and file size (default: 0)

  Note: The Pade spectrum is automatically computed if you use the 'fit' command.
        The 'pade' command is useful, if you only want the Pade spectrum in the first place without any fitting.

  Note: If the Pade spectrum was computed previously, it is not computed again but read from the previously generated files
        This is automatically enforced by setting 'Pade: {calc: false}' in the eval.yaml file after computing the Pade spectrum once.
        If, for some reason (e.g. you deleted the Pade files), you want to compute the Pade spectrum again, you have to set this option to 'true' again in eval.yaml of use the 'pade' command manually.

--------------------
```

#### New guess

```text
./BTevalSpec.py -h guess

--------------------
 Make a new guess based on the Fourier (default) or Pade spectrum)

  ./BTevalSpec.py [<gen-opt>] [--guess=pade [<pade opt>] [--thres=<thres>]| --guess=ft [<ft opt>] [--nsig=<nsig>]] [--range=<lb,rb>] [--imag] guess

    Creates an initial guess for the first batch of spectral lines without fitting and updates eval.yaml.
    The guess is based on the Pade (--guess=pade) or Fourier (--guess=ft, default) spectrum
    with the given <pade opt> or <ft opt> (if the latter were not computed previously).

    <gen-opt>         General options
    --guess=<guess>   Guess is based on the FT (<guess>=ft, default) or Pade (<guess>=pade) spectrum. If you use either, you can use the corresponding options of the 'ft' or 'pade' command, respectively.
    --thres=<thres>   If <guess>=pade: Real number in ]0,1] that determines how large a line in the Pade spectrum needs to be relative to the largest lines in order to be considered as for the initial guess. (default: 0.05, i.e., lines with a height of at least 5% of the largest one are considered)
                      Smaller <smooth> values (=decay rates) result in sharper lines (but still smooth in contrast to the Fourier spectrum). Too small values, however, can produce artifacts. (default: ln(100)/<propagation time>)
    --nsig=<nsig>     For the current iteration, determines the threshold for adding new lines in the following way: (default: 2.0)
                      The code takes the residue spectrum, which is essentially the Fourier spectrum minus the fit spectrum, which shows several peaks from the sine-cardinal shaped spectral lines.
                      In order to distinguish main peaks from side peaks, the code finds all peak in the given spectral range and computes the mean peak height (hbar) and the peak-height standard deviation (hsig).
                      Peaks that exceed hbar+<nsig>*hsig are accepted as new lines.
    --range=<lb,rb>   Determines the considered frequency range (in Ry) (default: 0.0,0.4 )
                      Attention: This option is written into eval.yaml and remembered for future runs.
    --imag            Only consider imaginary part. (default: true for boost, false else)
                      Attention: This option is written into eval.yaml and remembered for future runs. To undo it, set the 'imagonly' entry in eval.yaml to 'false' again.

  Note: The guess is automatically computed if you use the 'fit' command.
        The 'guess' command is useful, if you only want a guess without fitting.
--------------------
```

#### Remove excitation

```text
./BTevalSpec.py -h rm

--------------------
 Remove excitation without fit

  ./BTevalSpec.py [<gen-opt>] rm <exlist>

    Remove <exlist> from the list of excitations, where exlist is a comma-separated list of excitation labels, e.g., S1,S2,S3
    Alternatively, you can also remove the excitation manually from eval.yaml.
--------------------
```

#### Add excitation

```text
./BTevalSpec.py -h add

--------------------
 Add excitation without fitting

  ./BTevalSpec.py [<gen-opt>] [--nsig=<nsig> | --nadd=<nex> | --energy=<en>] [--nofix] add

    Add new excitations either at given energies (optional) or automatically using a sigma-threshold without fitting
    Before adding new excitations, the existing ones are fixed. Switch that off with the --nofix option.
    By deault, one batch of excitations is added with default parameters (as does the 'fit' command with '--niter=1' but without fitting)

    <gen-opt>         General options
    --nsig=<nsig>     For the current iteration, determines the threshold for adding new lines in the following way: (default: 2.0)
                      The code takes the residue spectrum, which is essentially the Fourier spectrum minus the fit spectrum, which shows several peaks from the sine-cardinal shaped spectral lines.
                      In order to distinguish main peaks from side peaks, the code finds all peak in the given spectral range and computes the mean peak height (hbar) and the peak-height standard deviation (hsig).
                      Peaks that exceed hbar+<nsig>*hsig are accepted as new lines.
    --nadd=<nex>      Add <nex> of the next largest excitations (default: 0, i.e., disabled)
    --energy=<en>     Add excitations at the given energies (as comma-separated list, e.g., <en>=0.160,0.171). Initial values for line heights are extracted from the spectrum (default: , i.e., disabled)
    --nofix           Old lines are automatically fixed (i.e., they gain a 'fix=True' in eval.yaml and are not fitted until released). The --nofix flag suppresses this behaviour.
    --range=<lb,rb>   Determines the considered frequency range (in Ry) (default: 0.0,0.4 )
                      Attention: This option is written into eval.yaml and remembered for future runs.
    --imag            Only consider imaginary part (default: true for boost, false else)
                      Attention: This option is written into eval.yaml and remembered for future runs. To undo it, set the 'imagonly' entry in eval.yaml to 'false' again.
    --wref=<wref>     Reference-omega for error suppression - suppresses the line-shape error of already fitted lines (that produce spurious lines in the residue spectrum) when searching for new lines.
                      Usually: small error suppression for <wref>=0.01, large suppression for <wref>=1.0. (default: 0.0)
--------------------
```

#### Reset excitation energy interval

```text
./BTevalSpec.py -h reset

--------------------
 Reset energy range to the standard pi/T interval around the current excitations' energy values

  ./BTevalSpec.py [<gen-opt>] reset

    <gen-opt>         General options
    When an excitation is added by the code, its energy parameter during the fit is restricted to a symmetric interval around the initial energy value to stabilize the fit.
    If you note that the rng significance of an excitation is not close to 1, you can try to reset the energy ranges once and, thus, give more freedom to the energy fit.
--------------------
```

#### Fit

```text
./BTevalSpec.py -h fit

--------------------
Fit

  ./BTevalSpec.py [<gen-opt>] [<ft opt>] [<guess opt>] [--reset] [--skip] [--single] [--signif] [--range=<lb,rb>] [--crit=<convcrit> | --nadd=<nadd>] [--nsig=<nsig>] [--niter=<niter>] [--imag] [--wref=<wref>] [--fitphase] fit

    Fit excitations following this scheme:
      1) Fit the present (from eval.yaml) excitations (use --skip to skip this step)
      2) Fix the present (old) excitations
      3) Add one batch of new excitations. How many excitations are added in one iteration is controlled by the --nsig (and/or the --nadd) parameter)
      4) Fit new excitations together (or singly if --single option is given)
      5) Release the old excitations
      6) Fit all excitations together
      7) Go to 2) until --niter iterations are done or --nadd excitations were added or the comprehensive error falls below the convergence criterion <convcrit>
      8) Compute significances (if --signif is given, also the computationally demanding ones are computed)

    <gen-opt>         General options
    <ft-opt>          'ft' command options, if ft is was not done before
    <guess-opt>       'guess' command options, if no lines are present yet
    --skip            Skip the first collective fit of existing excitations (useful if eval.yaml was changed manually)
    --nadd=<nadd>     Add up to <nadd> excitations (default: 0)
    --single          Instead of fitting all excitations at once, fit them one after the other from high to low strength 
                      Does a prior Fourier transform if not already done.
                      Does a prior Guess if not already done; requires a range option.
    --signif          Compute the significance-values that require fitting (costly)
    --fitphase        Use phase as fit parameter
                      Attention: This option is written into eval.yaml and remembered for future runs. To undo it, set the 'fitphase' entry in eval.yaml to 'false' again.
    --crit=<convcrit> Convergence criterion (only relevant if nadd==0 or not present): Add excitations until the fit error drops below this criterion. (default: 0., i.e., disabled)
    --niter           Maximum number of add-excitation iterations (in each iteration, the number of added excitations is determined by the <nsig> threshold)
    --reset           Reset energy range (cf. 'reset' command)
    --nsig=<nsig>     For the current iteration, determines the threshold for adding new lines in the following way: (default: 2.0)
                      The code takes the residue spectrum, which is essentially the Fourier spectrum minus the fit spectrum, which shows several peaks from the sine-cardinal shaped spectral lines.
                      In order to distinguish main peaks from side peaks, the code finds all peak in the given spectral range and computes the mean peak height (hbar) and the peak-height standard deviation (hsig).
                      Peaks that exceed hbar+<nsig>*hsig are accepted as new lines.
    --range=<lb,rb>   Determines the considered frequency range (in Ry) (default: 0.0,0.4 )
                      Attention: This option is written into eval.yaml and remembered for future runs.
    --imag            Only consider imaginary part (default: true for boost, false else)
                      Attention: This option is written into eval.yaml and remembered for future runs. To undo it, set the 'imagonly' entry in eval.yaml to 'false' again.
    --wref=<wref>     Reference-omega for error suppression - suppresses the line-shape error of already fitted lines (that produce spurious lines in the residue spectrum) when searching for new lines.
                      Usually: small error suppression for <wref>=0.01, large suppression for <wref>=1.0. (default: 0.0)

  Note: The simplest approach to fit a spectrum is to do
        1) one iteration at a time using
           ./BTevalSpec.py --skip [--nsig=<nsig>] --niter=1 fit
        2) then evaluate the added excitations (e.g. by inspecting the significance measures)
        3) backup the resulting eval.yaml
        4) go to 1) (use --skip to suppress the initial fitting of the input excitations, use the --nsig option to change the threshold for adding new excitation)

 Note: If any significance measure is not close to 1, this can (!) indicate that an excitation is erroneous. Try to remove the the excitation using
           ./BTevalSpec.py rm <excitation identifier>
       and call
           ./BTevalSpec.py fit
       to just fit the remaining excitations again
       If the rng significance of an excitation is not close to 1, resetting the energy range using --reset gives the fitting more freedom, which helps sometimes.

 Note: At the end of this command, the script prints what it would do next (i.e., which excitations it would add next).
       Don't be confused, it looks as if the script had actually done another fit iteration, however, it is just an information.

 Note: Several output files are generated that, e.g., contain the spectra, the fit, the fit-objective function (=kind of f-resolved fit error), a table of excitations, etc.
       Also, during the fit, the scripts opens a window with the current fit-objective and nsig threshold before fitting in each iteration to show what is actually done.
--------------------
```

#### Plotting

```text
./BTevalSpec.py -h plot

--------------------
 Plot

  ./BTevalSpec.py [<gen-opt>] [--range=<lb,rb>] [--imag] [--wref=<wref>] [--exclude=<excit>] plot <listOfMeasures>

    Plots <listOfMeasure> in {pade, ft, fit, err, spectrum} with excitations.

    <gen-opt>         General options
    <listOfMeasures>  [not supported yet] Measures to plot
    --exclude=<excit> [not supported yet] Comma-separated list of excitations that are excluded from the (fitted) data (default: '' )
    --range=<lb,rb>   Determines the considered frequency range (in Ry) (default: 0.0,0.4 )
                      Attention: This option is written into eval.yaml and remembered for future runs.
    --imag            Only consider imaginary part (default: true for boost, false else)
                      Attention: This option is written into eval.yaml and remembered for future runs. To undo it, set the 'imagonly' entry in eval.yaml to 'false' again.
    --wref=<wref>     Reference-omega for error suppression - suppresses the line-shape error of already fitted lines (that produce spurious lines in the residue spectrum) when searching for new lines.
                      Usually: small error suppression for <wref>=0.01, large suppression for <wref>=1.0. (default: 0.0)
--------------------
```

#### Fix excitations

```text
./BTevalSpec.py -h fix

--------------------
 Fix

  ./BTevalSpec.py [<gen-opt>] [--invert] fix {<listOfExcitations>|all}

    Fix the given list of excitations (by name). Those remain but are excluded from the fit (i.e., parameters are fixed).
    If you call with 'all' instead of a list of excitations, all excitations are fixed

    <gen-opt>         General options
    --invert: Do the same but inverted
    You can do this manually in eval.yaml as well by setting the 'fix=True' for the respective excitations.
--------------------
```

#### Release excitations

```text
./BTevalSpec.py -h release

--------------------
 Release

  ./BTevalSpec.py [<gen-opt>] [--invert] release {<listOfExcitations>|all}

    Release (un-fix) the given list of excitations (by name), cf. 'fix' command
    If you call with 'all' instead of a list of excitations, all excitations are released
    You can do this manually in eval.yaml as well by setting the 'fix=False' for the respective excitations.
--------------------
```

#### Compute transition densities

```text
./BTevalSpec.py -h decouple

--------------------
 Decouple

  ./BTevalSpec.py [<gen-opt>] [--jcalc=<calc-idx>] decouple

  Decouple given Fourier-transformed densities n(r,omega) at given omega to get the proper transition densities n_j.
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

  After the decoupling, you can the following output with one lines per excitation:

Transition density | abs. norm | real norm | imag norm | sin^2(dAng) | dAbs
S1                      0.0259     99.92 %      0.08 %      0.00        0.90 %
...

- Col 1) Excitation name
- Col 2) int |n_j| d^3r
- Col 3) Real-valued fraction (must be close to 100%)
- Col 4) Imag-valued fraction (must be close to   0%)
- Col 5) Compares the angle between the transition dipole from the fit and the one evaluated from the transition density (must be close to 0)
- Col 6) Compares the modulus of the transition dipole from the fit and the one evaluated from the transition density (must be close to 0)

If any of the columns 4-6 shows unrealistic values (significantly different from 0), something went wrong.
--------------------
```

### Examples/Testing

#### Overview

- Na2                 - Simple test using a global dipole moment from a laser excitation on a Na2 cluster. Also shows the error-supression.
- Na4_laser_transdens - Simple test using a global dipole moment from a laser excitation on a Na4 cluster with subsequent transition-density evaluation [cf. Schelter et al., JCTC 14, 1910 (2018)].
- Na4_boost           - Test using a global dipole moment from a boost excitation on a Na4 cluster.
- Na2DA               - Imperfect H-aggregate consisting of a Na2-Na2 donor-acceptor system with donor/acceptor specific dipole moments from a laser excitation with a gaussian profile to test fit to different spectral regions. Significances are considered there.
- B302                - Global dipole moment from a boost calculation on a bacteriochlorophyll. This test tries to fit a broad excitation band with many close-lying excitations.
- B302_2calc          - This test is comparable with B302 but uses dipole moments from two different calculations with different external-field polarizations.

#### Instructions

Copy the `testTemplates/` directory into a directory, e.g., `test~/` (the tilde makes git ignore this directory)

```bash
cp -rP testTemplates test~
```

for actual testing.
The `-P` option leaves symbolic links as such.

Confer the `readme.md` files in the single test directories for instructions.
