reset

#n=1
#d="x"
#r=2

set xrange [0.1:0.28]
set yrange [-330:330]
set y2tics

DAT = "dipole_global1_0".n."_ft_".d.".dat
FIT = "dipole_global1_0".n."_fit_".d.".dat

set title "Area ".n." | Direction ".d." | RC ".r
plot DAT u 1:r w l, FIT u 1:r w l
