#!/bin/bash

FIN=$1
FOUT="excit.dat"
TMP1="tmp1.dat"
TMP2="tmp2.dat"
TMP3="tmp3.dat"
TMP4="tmp4.dat"

grep 'energy:'   $FIN | awk '{print $2}' > $TMP1
grep 'strength:' $FIN | awk '{print $2}' > $TMP2
grep 'signif:'   $FIN | awk '{print $2}' > $TMP3
paste $TMP1 $TMP2 $TMP3 > $TMP4
echo "# energy | strength | signif" > $FOUT
cat $TMP4 >> $FOUT
rm $TMP1 $TMP2 $TMP3 $TMP4
