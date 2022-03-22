reset
set sample 1000
set xrange [3.0:3277.0]
set xlabel 'Energy (eV)'
set ylabel 'Oscillator strength'

###############################################################################

ev = 13.605693122994
n = 0.025

#Excitations
w_1 = 0.22602672335797466
f_1 = 0.016714530515883907
w_2 = 0.24321840384066531
f_2 = 0.279347385704451
w_3 = 0.2481873454609313
f_3 = 0.0038954563210598482
w_4 = 0.2881213411332013
f_4 = 0.11190774956181075
w_5 = 0.34894863899743406
f_5 = 0.0025686380580015344
w_6 = 149.92181052654146
f_6 = 3650.0056227432747
w_7 = 240.8076405080726
f_7 = 49707.240954214845

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