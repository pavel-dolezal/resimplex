#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
set yr [-300:200]
#set xr [29:48]
#set xr [17:28]
#set xr [:17]
#set xr [0:48]
set xr [46:]

set palette rgbformulae 33,13,10

d1=1
d2=1
d3=1
d4=1
d5=1
d6=1
d7=1
d8=1
d9=1
d10=1
d11=1
d12=1
d13=0.00001
d14=0.00001
f(x) = d1+d2*x+d3*x**2+d4*x**3+d5*x**4
#+d6*x**5+d7*x**6
#+d8*x**7+d9*x**8+d10*x**9+d11*x**10+d12*x**11+d13*x**12+d14*x**13
fit f(x) 'delta.asc' u 1:2 via d1,d2,d3,d4,d5
#,d6,d7
#,d8,d9,d10,d11,d12,d13,d14

save fit "delta_fit.asc"

plot	"delta.asc" u 1:2 t "delta" w lp ,\
	"delta.asc" u 1:(f($1)) w l t "fit" lc 'blue' ,\
#	"delta2.asc" u 1:2 t "old" w lp,\

pa -1

#set term png small
#set out "delta.png"
#replot

q
