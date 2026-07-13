#!/usr/bin/gnuplot

set colors classic
set term x11

set xl "iter"
set yl "chi^2"

set logscale y
set ytics 10
set mytics 10
set key outside samplen 1.5


plot "<awk '{ i++; print i, $(NF-1); }' chi2_funcB.tmp" u 1:2 t "chi2" w lp lt 1 ps 0.5,\

pa -1

replot

q
