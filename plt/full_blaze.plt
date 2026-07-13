#!/usr/bin/gnuplot

set terminal postscript eps enhanced
set terminal x11
set colors classic 

#set title "blaze function"
set xlabel "{/Symbol l} - {/Symbol l}_c [{\305}]" font ",12"
set ylabel "I [ADU]" font ",12"

set grid
set palette rgbformulae 33,13,10

plot "<awk '($5 == 1) {print $2*10**10, $(NF)}' Spectra.dat" u 1:2 w l t 'blaze value',\
     "../../input_for_rs1+2/ktc00018.txt" u 1:2 t 'ktc00018' w l
pause -1

#set term png small font "{,12}"
#set out "blaze1.png"
#replot

q
