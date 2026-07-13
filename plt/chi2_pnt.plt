#!/usr/bin/gnuplot

set colors classic
#set term qt font "Arial,8"
set term x11

ang = 1.e-10  # m
nm = 1.e-9  # m

set xl "lambda [ang]"
set yl "I_lambda [] (shifted by dataset number)"

set yr [0.0:1.1]
set xr [4500:4800]

fac = 0.25
set ytics fac
set mytics 1
set grid ytics mytics
set zeroaxis
set bar 0.5

bind "k" "unset xrange; unset yrange; replot"


set palette model RGB defined (1 "blue", 2 "red")


p \
  "<awk '{print $1, $2}' ktc00018.asc" u 1:2 t 'spectrum' w l lc 'blue',\
  "<awk '{print $1, $6/20}' ktc00018.asc" u 1:2 t 'chi2' w l pt 7 ps 0.8 lc 'red' ,\
  "<awk '($4==17) {print $2*10**10, $3}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '{print $1, $2+0.1, $6}' ktc00058.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.2, $6}' ktc00087.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.3, $6}' ktc00123.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.4, $6}' ktc00143.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.5, $6}' ktc00179.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.6, $6}' ktc00248.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.7, $6}' ktc00282.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.8, $6}' ktc00303.asc" u 1:2:3 not w l lc palette z,\
#  "<awk '{print $1, $2+0.9, $6}' ktc00337.asc" u 1:2:3 not w l lc palette z,\
 # "<awk '{print $1, $2}' ktc00337-10spectra.asc" u 1:2 t '10-spectra' w l lc 'red',\

 # "<awk '($4==57) {print $2*10**10, $3+0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "Spectra.dat" u ($2*10**10):7 w l lc 'black',\
#  "<awk '($5!=i){ print s; } { print; i=$5; }' ../../resimplex2-7-no-eps-degree2/output/Spectra.dat" u ($2/ang):($3+fac*($5-1)) t 'with lines' w l lt 3 lc 'red',\


pa -1

set term png small
set out "Spectra.png"
rep

q

  "Spectra.dat" u ($2/nm):($3+1*$5):4 w err lt 3 pt 1 ps 0

