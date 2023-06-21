reset

FILE='eval.yaml.excit.dat'

plot \
	FILE u 1:3:2:4 with errorbars ti 'Strength',\
	FILE u 1:5 ti 'signifFit',\
	FILE u 1:6 ti 'signifAng',\
	FILE u 1:7 ti 'signifExc'
