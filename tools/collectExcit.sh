#!/bin/bash

grep 'energy:'      $1 | awk '{print $2}' > energy.dat
grep 'energyErr:'   $1 | awk '{print $2}' > energyErr.dat
grep 'strength:'    $1 | awk '{print $2}' > strength.dat
grep 'strengthErr:' $1 | awk '{print $2}' > strengthErr.dat
grep 'signifFit:'   $1 | awk '{print $2}' > signifFit.dat
grep 'signifAng:'   $1 | awk '{print $2}' > signifAng.dat
grep 'signifExc:'   $1 | awk '{print $2}' > signifExc.dat
grep 'signifErr:'   $1 | awk '{print $2}' > signifErr.dat
grep 'signifRng:'   $1 | awk '{print $2}' > signifRng.dat

paste energy.dat strength.dat energyErr.dat strengthErr.dat signifFit.dat signifAng.dat signifExc.dat signifErr.dat signifRng.dat > $1.tmp
echo "# energy | strength | energyErr | strengthErr | signifFit | signifAng | signifExc | signifErr | signifRng" > $1.excit.dat
cat $1.tmp >> $1.excit.dat
rm $1.tmp

#gnuplot --persist -e "FILE='$1.excit.dat'" plotExcit.plt
