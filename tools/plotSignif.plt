reset

FILE='eval.yaml.excit.dat'

set xlabel "Sj"

plot \
	FILE u ($0-1):5 ti 'signifFit',\
	FILE u ($0-1):6 ti 'signifAng',\
	FILE u ($0-1):7 ti 'signifExc',\
	FILE u ($0-1):8 ti 'signifErr'
