reset

FILE1='eval.yaml.2.res.excit.dat'
FILE2='eval.yaml.3.res.excit.dat'

set xlabel "Sj"

plot \
	FILE1 u ($0-1):5 ti 'signifFit',\
	FILE1 u ($0-1):6 ti 'signifAng',\
	FILE1 u ($0-1):7 ti 'signifExc',\
	FILE1 u ($0-1):8 ti 'signifErr',\
	FILE2 u ($0-1):5 ti 'signifFit',\
	FILE2 u ($0-1):6 ti 'signifAng',\
	FILE2 u ($0-1):7 ti 'signifExc',\
	FILE2 u ($0-1):8 ti 'signifErr'
