#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
#set yr [-2:7]

#set xr [:17]
#set xr [15:29]
#set xr [27:48]
set xr [27:43]
#set xr [41:]
#set xr [46:]
#set xr [:48]

set palette rgbformulae 33,13,10

l1 = 2380.8435548261323        
l2 = 3.9218744126266972E-003  
l3 = -1.0692775259323601        
l4 = 151.38734288778659       
l5 = -4.7455440029226832       
l6 = -1.3860422798997006E-002   
l7 = 2.4593702208279225E-004  
l8 = -3.0377329791284924E-006
l9 = 0.00000001
l10 = 0.00000000001
l11 = 0.000000000000001
l12 = 0.0000001
l13 = 0.0000000001
l14 = 0.00000000001
l15 = 0.00000000000001
l16 = 0.0000000000000001
l17 = 0.0000000000000001 
l18 = 0.000000000000000001
l19 = 0.000000000000000000001
l20 = 0.0000000000000000000001
l21 = 0.0000000000000000000000001

g1 = 2380.8435548261323        
g2 = 3.9218744126266972E-003  
g3 = -1.0692775259323601        
g4 = 151.38734288778659       
g5 = -4.7455440029226832       
g6 = -1.3860422798997006E-002   
g7 = 2.4593702208279225E-004  
g8 = -3.0377329791284924E-006
s1 = 0.15
s2 = 0.3
s3 = 0.0001
s4 = 0.00001

#f(x) = l1*tan(l2*x-l3)+l4+l5*x+l6*x**2+l7*x**3+l8*x**4+l9*x**5+s1*sin(s2*x-s3)
f(x) = l1 + l2*x + l3*x**2 + l4*x**3 + l5*x**4 + l6*x**5 + l7*x**6 + l8*x**7 + l9*x**8

#+l10*x**6+l11*x**7+l12*x**8+l13*x**9+l14*x**10+l15*x**11+l16*x**12+l17*x**13+l18*x**14+l19*x**15+l20*x**16+l21*x**17+
#+l10*x**6
#+l11*x**7
#g(x) = g1*tan(g2*x-g3)+g4+g5*x+g6*x**2+g7*x**3+g8*x**4
#+g9*x**5+g10*x**6+g11*x**7

#fit f(x) 'lambda_c.asc' u 1:2 via l1,l2,l3,l4,l5,l6,l7,l8,l9,s1,s2,s3
fit f(x) 'lambda_c.asc' u 1:2 via l1,l2,l3,l4,l5,l6,l7,l8,l9
#,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20,l21
#,l7,l8,l9
#,l10,l11
#,l9
#,s1,s2,s3
#,l10
#,l11

save fit "lambda_c_fit.asc"

#fit g(x) 'lambda_c.asc' u 1:2 via g1,g2,g3,g4,g5,g6,g7,g8

#save fit "lambda_c_fit_exp.asc"

plot "lambda_c.asc" u 1:(f($1)-$2) w l t "gnuplot-fit" lc 'red' ,\
#     "lambda_c.asc" u 1:(h($1)) w l t "sinus" lc 'blue' ,\
#     "lambda_c.asc" u 1:(g($1)-$2) w l t "resimplex3" lc 'blue' ,\

pa -1

#set term png small
#set out "lambda_c.png"
#replot

q
