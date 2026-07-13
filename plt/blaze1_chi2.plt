#!/usr/bin/gnuplot
set terminal postscript eps enhanced
set terminal x11

# =============================================
# run addchi2toblaze.sh to generate input files
# =============================================

set xlabel "{/Symbol l} - {/Symbol l}_c [{\305}]" font ",12"
set ylabel "I_{/Symbol l} [ADU] " font ",12"
set zlabel "chi2" font ",20"

set grid
set palette rgbformulae 33,13,10

set logscale cb 10
set cbtics 5 logscale
set cbr [1000:50000]

plot for [order= 1:62] "<awk -v order=". order. " '$1 == order && $NF == 7' blaze1_chi2.tmp" u 2:3:4 w l lc palette z lw 1.2 t order

pa -1

set term png small
set out "blaze_m_const.png"
replot

q
