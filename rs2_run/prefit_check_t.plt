#!/usr/bin/gnuplot

set colors classic
set term x11

t1 = 84.3290946424168 
t2 = -46.872236222112 
t3 = 3.73592405100025 
t4 = -0.0973022041341739 
t5 = 0.000794060040210317

f(x) = t1+t2*x+t3*x**2+t4*x**3+t5*x**4
fit f(x) 'prefit_check_t.tmp' u 1:2 via t1,t2,t3,t4,t5

plot "prefit_check_t.tmp" u 1:2 t 'free theta' w p lc 'red',\
	"prefit_check_t.tmp" u 1:3 t 'resimplex prefit' w l lc 'blue',\
	"prefit_check_t.tmp" u 1:(f($1)) t 'gnuplot prefit' w l lc 'black' ,\

pa -1


