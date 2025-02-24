# Test Na2

Simple test using a global dipole moment from a laser excitation on a Na2 cluster.
 - BTevalSpec.py         - Link to the script
 - dipole_global1_01.dat - dipole-moment data (confer header for meta information on the calculation)
 - eval.yaml.0           - initial eval.yaml: (i) call `./BTevalSpec.py new` and (ii) adding the dipole- and laser-profile files, set invertPhase to true, and adjust the fit range
 - laser_profile.dat     - time-dependent laser profile (cf. meta data in the dipole file for information about  polarization etc.)
 - plotFit.plt           - A simple gnuplot script for plotting the fit called by plotFit.sh
 - plotFit.sh            - A script that calls plotFit.plt to plot different contributions to the spectrum

Call
 -  or
 - call `./BTevalSpec.py [-v 1] -f eval.yaml.guess guess` to create a new eval.yaml file, add the dipole- and laser-profile files, adjust the fit range and call `./BTevalSpec.py fit` to fit.

 - start from scratch: call `./BTevalSpec.py new` to create a new eval.yaml file and add the dipole- and laser-profile files, adjust the fit range,  and call `./BTevalSpec.py fit` to fit or

