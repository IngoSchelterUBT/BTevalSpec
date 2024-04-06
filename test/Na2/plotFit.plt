reset

#n=1
#d="x"
#r=2

set xrange [0.12:0.18]
set yrange [-100:100]
set y2tics

DAT = "dipole_global1_0".n."_ft_".d.".dat
FIT = "dipole_global1_0".n."_fit_".d.".dat
LAS = "laser_profile_ft.dat"

set title "Area ".n." | Direction ".d." | RC ".r
plot DAT u 1:r w l, FIT u 1:r w l, LAS u 1:(atan2($3,$2)) axis x1y2 w l, LAS u 1:(sqrt($2**2+$3**2)) w l
