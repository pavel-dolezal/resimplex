#!/usr/bin/gnuplot

set colors classic
set terminal postscript eps enhanced
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"
#set yr [-8:10]
set xr [15:29]

set palette rgbformulae 33,13,10


l1 =3028.6755672893028        
l2 = 4.6064662786312006E-003 
l3 = -0.98134683012185042       
l4 = -34.897859240575073       
l5 = -9.4798152066424670       
l6 = -3.6927177217153627E-002 
l7 =  1.4152525286865268E-004  
l8 = -6.5166189984160432E-006

f(x) = l1*tan(l2*x-l3)+l4+l5*x+l6*x**2+l7*x**3+l8*x**4

fit f(x) 'lambda_c.asc' u 1:6 via l1,l2,l3,l4,l5,l6,l7,l8

plot	"lambda_c.asc" u 1:6 t "lambda-c" w lp ,\
	"lambda_c.asc" u 1:(f($1)) w l t "fit" lc 'blue' ,\
#	"lambda_c2.asc" u 1:6 t "old" w lp ,\
#plot "lambda_c.asc" u 1:(f($1)-$2) w l t "residua" lc 'red'

pa -1

#set term png small
#set out "lambda_c.png"
#replot

q
