#!/usr/bin/gnuplot

set colors classic
set term x11

e1 = -3.55653322391977 
e2 = -6.66832855743865 
e3 = 0.644119212700238 
e4 = -0.0181860057337128 
e5 = 0.000157231465904982

f(x) = e1 + e2*x + e3*x**2 + e4*x**3 +e5*x**4
fit f(x) 'prefit_check_e.tmp' u 1:2 via e1,e2,e3,e4,e5

plot "prefit_check_e.tmp" u 1:2 t 'free epsilon' w p lc 'red',\
	"prefit_check_e.tmp" u 1:3 t 'resimplex prefit' w l lc 'blue',\
	"prefit_check_e.tmp" u 1:(f($1)) t 'gnuplot prefit' w l lc 'black' ,\

pa -1


