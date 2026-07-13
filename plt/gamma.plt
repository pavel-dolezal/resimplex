#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
#set yr [-800:400]
#set xr [29:48]
#set xr [17:28]
#set xr [:17]
#set xr [0:48]
set xr [46:]

set palette rgbformulae 33,13,10

g1=1
g2=1
g3=1
g4=1
g5=1
g6=1
g7=1
g8=1
g9=1
g10=1
g11=1
g12=1
g13=1
f(x) = g1+g2*x+g3*x**2+g4*x**3+g5*x**4
#+g6*x**5+g7*x**6
#+g8*x**7+g9*x**8+g10*x**9+g11*x**10+g12*x**11+g13*x**12
fit f(x) 'gamma.asc' u 1:2 via g1,g2,g3,g4,g5
#,g6,g7
#,g8,g9,g10,g11,g12,g13
save fit "gamma_fit.asc"

a1=-114.77216133424587       
a2=0.55598787365005808       
a3=0.63273090675306776       
a4=-1.6433140520886902E-002   
a5=1.1646659455917733E-004
g(x)=a1+a2*x+a3*x**2+a4*x**3+a5*x**4

b1=20.230849014315602        
b2=3.7350902308198712        
b3=1.9280172056563294E-002  
b4=-1.6118247576501158E-003   
b5=6.8777471475951756E-006
h(x)=b1+b2*x+b3*x**2+b4*x**3+b5*x**4

plot	"gamma.asc" u 1:2 t "gamma" w lp ,\
	"gamma.asc" u 1:(f($1)) w l t "fit" lc 'blue' ,\
	#"gamma.asc" u 1:(g($1)) t "gamma4" w lp,\
	#"gamma.asc" u 1:(h($1)) t "g4 from zero" w lp,\

pa -1

#set term png small
#set out "gamma.png"
#replot

q
