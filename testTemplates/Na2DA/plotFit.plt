reset

#n=1
#d="x"
#r=2

set xrange [0.12:0.24]
set yrange [-5:5]

DAT = "dipole_atomgroups1_0".n."_ft_".d.".dat
FIT = "dipole_atomgroups1_0".n."_fit_".d.".dat

set title "Area ".n." | Direction ".d." | RC ".r
plot DAT u 1:r w l, FIT u 1:r w l
