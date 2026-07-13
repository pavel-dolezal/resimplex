#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
set yr [-200:200]
#set xr [29:48]
#set xr [17:28]
#set xr [:17]
#set xr [0:48]
#set xr [46:]

set palette rgbformulae 33,13,10

e1=1
e2=1
e3=1
e4=1
e5=1
e6=1
e7=1
e8=1
e9=1
e10=1
e11=1
e12=1
e13=1
f(x) = e1+e2*x+e3*x**2+e4*x**3+e5*x**4
#+e6*x**5+e7*x**6
#+e8*x**7+e9*x**8
#+e10*x**9+e11*x**10+e12*x**11+e13*x**12
fit f(x) 'epsilon.asc' u 1:2 via e1,e2,e3,e4,e5
#,e6,e7
#,e8,e9
#,e10,e11,e12,e13

save fit "epsilon_fit.asc"

#,a7
#,a8,a9,a10,a11,a12,a13
plot	"epsilon.asc" u 1:2 t "epsilon" w lp ,\
	"epsilon.asc" u 1:(f($1)) w l t "fit" lc 'blue' ,\
#	"epsilon2.asc" u 1:2 t "old" w lp ,\

pa -1

#set term png small
#set out "epsilon.png"
#replot

q
