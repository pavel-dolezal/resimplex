#!/usr/bin/gnuplot
set colors classic
set terminal postscript eps enhanced
set term x11
set key bottom
set xl "r"
set yl "{/Symbol a}"
#set yr [-8:10]
#set xr [0:55]
set palette rgbformulae 33,13,10
l1 = 52500
l2  = 1400
l3  = -74
l4  = -27
l5  = 4
l6  = 0.02
l7  = 0.001 
l8  = 0.0001
l9  = 0.00001
l10 = 0.000001
l11 = 0.0000001
l12 = 0.00000001
l13 = 0.000000001
l14 = 0.0000000001
l15 = 0.00000000001
l16 = 0.000000000001
l17 = 0.0000000000001
l18 = 0.00000000000001
l19 = 0.000000000000001
f(x) = l1+l2*x+l3*x**2+l4*x**3+l5*x**4+l6*x**5+l7*x**6+l8*x**7+l9*x**8+l10*x**9+l11*x**10+l12*x**11+l13*x**12+l14*x**13+l15*x**14+l16*x**15+l17*x**16+l18*x**17+l19*x**18
#+l20*x**19+l21*x**20+l22*x**21

fit f(x) 'max_blaze.asc' u 1:2 via l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19
#,l20,l21,l22

save fit "max_blaze_fit.asc"

#plot	"max_blaze.asc" u 1:(f($1)-$2) t "residua" w lp ,\

plot "max_blaze.asc" u 1:(f($1)) w l t "gnuplot fit" ,\
     "max_blaze.asc" u 1:2 w l t "max_blaze" ,\
#     "max_blaze.asc" u 1:(g($1)) w l t "new-poly",\

pa -1

#set term png small
#set out "max_blaze.png"
#replot

q
