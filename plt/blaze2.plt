#!/usr/bin/gnuplot

set terminal x11

set title "blaze function"
set xlabel "lambda - lambda-c [A]"
set ylabel "flux"

set grid
set palette rgbformulae 33, 13, 10


plot for [order= 1:62] "<awk -v order=". order. " '$1 == order && $NF == 1' blaze2.tmp" u 2:3:1 w l lc palette z lw 3 t order


pause -1
