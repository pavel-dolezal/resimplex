#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
#set yr [-500:400]
#set xr [29:48]
#set xr [17:28]
#set xr [:17]
#set xr [0:48]
#set xr [46:]

set palette rgbformulae 33,13,10

t1=1
t2=1
t3=1
t4=1
t5=1
t6=1
t7=1
t8=1
t9=1
t10=1
t11=1
t12=1
t13=1
f(x) = t1+t2*x+t3*x**2+t4*x**3+t5*x**4
#+t6*x**5+t7*x**6
#+t8*x**7+t9*x**8+t10*x**9+t11*x**10+t12*x**11+t13*x**12
fit f(x) 'theta.asc' u 1:2 via t1,t2,t3,t4,t5
#,t6,t7
#,t8,t9,t10,t11,t12,t13

save fit "theta_fit.asc"

#,a5,a6
#,a7
#,a8,a9,a10,a11,a12,a13
plot	"theta.asc" u 1:($2) t "theta" w lp ,\
	"theta.asc" u 1:(f($1)) w l t "fit" lc 'blue',\
#	"theta2.asc" u 1:2 t "old" w lp ,\

pa -1

#set term png small
#set out "theta.png"
#replot

q
