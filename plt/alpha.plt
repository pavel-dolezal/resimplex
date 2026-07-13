#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
set yr [0.6:1.6]
#set xr [0:48]
#set xr [29:48]
#set xr [17:28]
#set xr [:55]
#set xr [45:]

set palette rgbformulae 33,13,10

a1=1
a2=10
a3=10
a4=100
a5=1
a6=1
a7=1
a8=1
a9=1
a10=1
a11=1
a12=1
a13=1
f(x) = a1+a2*x+a3*x**2+a4*x**3+a5*x**4
#+a6*x**5
#+a7*x**6+a8*x**7
#+a9*x**8
#+a10*x**9+a11*x**10+a12*x**11+a13*x**12
fit f(x) 'alpha.tmp' u 1:2 via a1,a2,a3,a4,a5
#,a6
#,a7,a8
#,a9
#,a10,a11,a12,a13
save fit "alpha_fit.asc"
plot	"alpha.tmp" u 1:2 t "alpha skip" w lp ,\
	"alpha.tmp" u 1:(f($1)) w l t "fit" lc 'blue' ,\
#	"alpha2.asc" u 1:2 t "a no skip" w lp,\

pa -1

#set term png small
#set out "alpha.png"
#replot

q
