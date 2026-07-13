#!/usr/bin/gnuplot

set colors classic
set term x11

g1 = -47.6314157257169 
g2 = 6.16085169487821 
g3 = -0.39279075085156 
g4 = 0.0129306784606281 
g5 = -0.000130973540092256

f(x) = g1+g2*x+g3*x**2+g4*x**3+g5*x**4
fit f(x) 'prefit_check_g.tmp' u 1:2 via g1,g2,g3,g4,g5

plot "prefit_check_g.tmp" u 1:2 t 'free gamma' w p lc 'red',\
	"prefit_check_g.tmp" u 1:3 t 'resimplex prefit' w l lc 'blue',\
	"prefit_check_g.tmp" u 1:(f($1)) t 'gnuplot prefit' w l lc 'black',\

pa -1


