reset

FILE='eval.yaml.6.res.excit.dat'

plot \
	FILE u 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars ti 'Strength'
