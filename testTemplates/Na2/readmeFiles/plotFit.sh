#!/bin/bash

#n=$1
#d="$2"
#r=$3

gnuplot --persist -e "n='1'" -e "d='z'" -e "r=2" plotFit.plt
gnuplot --persist -e "n='1'" -e "d='z'" -e "r=3" plotFit.plt
