#!/usr/bin/gnuplot

set colors classic
set term x11

ang = 1.e-10  # m
nm = 1.e-9  # m

set xl "lambda [ang]"
set yl "I_lambda [] (shifted by dataset number)"

#set yr [0.9:1.25]
#set xr [5000:5300]

#set xr [5000:6000]
#set yr [0.67:1.5]

fac = 0.25
set ytics fac
set mytics 1
set grid ytics mytics
set zeroaxis
set bar 0.5

folder_name = system("pwd")
set title "folder:" . folder_name

set palette model RGB defined (1 "blue", 2 "red")

bind "k" "unset xrange; unset yrange; replot"
bind "1" "set xr [4500:4800]; replot"
bind "2" "set xr [4800:4930]; replot"
bind "3" "set xr [5300:6100]; replot"
bind "4" "set xr [6450:6680]; replot"
bind "5" "set xr [5000:5200]; replot"
bind "6" "set xr [4490:8900]; replot"
bind "7" "set xr [6150:6525]; replot"


set label "split2" at 5812.0,graph 1.01 center; set arrow from 5812.0,graph 0 rto 0,graph 1 nohead lt 0
set label "split1" at 5178.0,graph 1.01 center; set arrow from 5178.0,graph 0 rto 0,graph 1 nohead lt 0

p \
  "<awk '{print $1, $2, $(NF)}' ktc00018.asc" u 1:2:3 not w l lc palette z,\

pa -1

#set term png small
#set out "Spectra.png"
#rep

q

