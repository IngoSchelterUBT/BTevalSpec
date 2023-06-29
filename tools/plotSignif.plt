reset

FILE='eval.yaml.1.res.excit.dat'

set xlabel "Sj"

plot \
	FILE u 0:5 ti 'signifFit',\
	FILE u 0:6 ti 'signifAng',\
	FILE u 0:7 ti 'signifExc',\
	FILE u 0:8 ti 'signifErr',\
	FILE u 0:9 ti 'signifRng'
