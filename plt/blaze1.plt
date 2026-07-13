#!/usr/bin/gnuplot

set terminal postscript eps enhanced
set terminal x11

set xlabel "{/Symbol l} - {/Symbol l}_c [{\305}]" font ",12"
set ylabel "I [ADU]" font ",12"

set yr [10000:]

set xr [-40:80]

set grid
set palette rgbformulae 33,13,10

plot for [order= 1:62] "<awk -v order=". order. " '$1 == order && $NF == 1' blaze1.tmp" u 2:3:1 w l lc palette z lw 1.2 t order

pause -1

#set term png small font "{,12}"
#set out "blaze1.png"
#replot

q
