reset

FILE1='eval.yaml.1.res.excit.dat'
FILE2='../B302/eval.yaml.1.res.excit.dat'

plot \
	FILE1 u 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars ti 'y-boost' ,\
	FILE2 u 1:2:($1-$3):($1+$3):($2-$4):($2+$4) with xyerrorbars ti 'xz-boost'
