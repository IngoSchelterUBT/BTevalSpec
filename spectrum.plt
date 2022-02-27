reset
set sample 1000
set xrange [3.0:103.0]
set xlabel 'Energy (eV)'
set ylabel 'Oscillator strength'

###############################################################################

ev = 13.605693122994
n = 0.025

#Excitations
w_1 = 0.2259951836201931
f_1 = 0.016679419524174198
w_2 = 0.2432196317128414
f_2 = 0.27974417998403334
w_3 = 0.24825646068787896
f_3 = 0.0034937473623454078
w_4 = 0.28812853484268325
f_4 = 0.11158985873526993
w_5 = 0.3553641750763505
f_5 = 0.01867908074225572
w_6 = 4.715422186843155
f_6 = -11187.745432423866
w_7 = 7.539260463066363
f_7 = 28727.49684155919

#Single lines
s_1(x) = f_1*exp(-((x/ev-w_1)/(n/ev))**2)
s_2(x) = f_2*exp(-((x/ev-w_2)/(n/ev))**2)
s_3(x) = f_3*exp(-((x/ev-w_3)/(n/ev))**2)
s_4(x) = f_4*exp(-((x/ev-w_4)/(n/ev))**2)
s_5(x) = f_5*exp(-((x/ev-w_5)/(n/ev))**2)
s_6(x) = f_6*exp(-((x/ev-w_6)/(n/ev))**2)
s_7(x) = f_7*exp(-((x/ev-w_7)/(n/ev))**2)

#Spectrum
sp(x) = \
        s_1(x) + s_2(x) + s_3(x) + s_4(x) + s_5(x) + s_6(x) + s_7(x)

unset arrow

#Arrows
set arrow from w_1*ev,0 to w_1*ev,f_1 lw 2 head filled
set label 'S' at w_1*ev,f_1 offset 0,0.5
set arrow from w_2*ev,0 to w_2*ev,f_2 lw 2 head filled
set label 'S' at w_2*ev,f_2 offset 0,0.5
set arrow from w_3*ev,0 to w_3*ev,f_3 lw 2 head filled
set label 'S' at w_3*ev,f_3 offset 0,0.5
set arrow from w_4*ev,0 to w_4*ev,f_4 lw 2 head filled
set label 'S' at w_4*ev,f_4 offset 0,0.5
set arrow from w_5*ev,0 to w_5*ev,f_5 lw 2 head filled
set label 'S' at w_5*ev,f_5 offset 0,0.5
set arrow from w_6*ev,0 to w_6*ev,f_6 lw 2 head filled
set label 'S' at w_6*ev,f_6 offset 0,0.5
set arrow from w_7*ev,0 to w_7*ev,f_7 lw 2 head filled
set label 'S' at w_7*ev,f_7 offset 0,0.5

#Plot
plot \
      sp(x) lc rgb 'black' dt 1 lw 2 ti 'Spectrum',\
      s_1(x) lc 4 dt 2 lw 2 ti 'Excitations',\
      s_2(x) lc 4 dt 2 lw 2 noti,\
      s_3(x) lc 4 dt 2 lw 2 noti,\
      s_4(x) lc 4 dt 2 lw 2 noti,\
      s_5(x) lc 4 dt 2 lw 2 noti,\
      s_6(x) lc 4 dt 2 lw 2 noti,\
      s_7(x) lc 4 dt 2 lw 2 noti,\
      0 lc rgb 'black' lw 2 noti


set output