reset

plot \
	'< paste dipole_global1_01_xz_ft_x.dat dipole_global1_01_xz_fit_x.dat' u 4:($3-$6) w l ,\
	'< paste dipole_global1_01_xz_ft_y.dat dipole_global1_01_xz_fit_y.dat' u 4:($3-$6) w l ,\
	'< paste dipole_global1_01_xz_ft_z.dat dipole_global1_01_xz_fit_z.dat' u 4:($3-$6) w l ,\
	'< paste dipole_global1_01_y_ft_x.dat  dipole_global1_01_y_fit_x.dat'  u 4:($3-$6) w l ,\
	'< paste dipole_global1_01_y_ft_y.dat  dipole_global1_01_y_fit_y.dat'  u 4:($3-$6) w l ,\
	'< paste dipole_global1_01_y_ft_z.dat  dipole_global1_01_y_fit_z.dat'  u 4:($3-$6) w l 
