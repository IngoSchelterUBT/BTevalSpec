reset

plot \
	'B302/eval.yaml.0.res.excit.dat'       u 1:2:3:4 with errorbars ti 'B302 - 1 xz'     ,\
	'B302_2calc/eval.yaml.0.res.excit.dat' u 1:2:3:4 with errorbars ti 'B302 - 2 xz & y'
