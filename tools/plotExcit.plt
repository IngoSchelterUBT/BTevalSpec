reset

FILE='eval.yaml.excit.dat'

plot \
	FILE u 1:3:2:4 with errorbars ti 'Strength'
