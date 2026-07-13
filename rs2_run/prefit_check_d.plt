#!/usr/bin/gnuplot

set colors classic
set term x11

d1 = 1.65632377936872 
d2 = 0.218792776045715 
d3 = 0.0604985543548616 
d4 = -0.00365884507002047 
d5 = 4.45283016821683e-05

f(x) = d1 + d2*x + d3*x**2 + d4*x**3 + d5*x**4
fit f(x) 'prefit_check_d.tmp' u 1:2 via d1,d2,d3,d4,d5

plot "prefit_check_d.tmp" u 1:2 t 'free delta' w p lc 'red',\
	"prefit_check_d.tmp" u 1:3 t 'resimplex prefit' w l lc 'blue',\
	"prefit_check_d.tmp" u 1:(f($1)) t 'gnuplot prefit' w l lc 'black' ,\

pa -1


