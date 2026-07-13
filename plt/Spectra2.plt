#!/usr/bin/gnuplot

set colors classic
set term qt

ang = 1.e-10  # m
nm = 1.e-9  # m

set xl "lambda [ang]"
set yl "I_lambda [] (shifted by dataset number)"

#set yr [3.15:3.3]
#set yr [2.25:3.5]
#set xr [5650:5915]
set xr [4490:8900]

folder_name = system("pwd")
#folder_name = system("basename `pwd`") 
set title "folder:" . folder_name

#set xr [5881:6585]
#set yr [2.57:2.81]

fac = 0.25
set ytics fac
set mytics 1
set grid ytics mytics
set zeroaxis
set bar 0.5

#call "line.plt" "split2" 5812.0
#call "line.plt" "split1" 5178.0

set label "47" at 7127.0,graph 1.01 center; set arrow from 7127.0,graph 0 rto 0,graph 1 nohead lt 0
set label "48" at 7219.0,graph 1.01 center; set arrow from 7212.0,graph 0 rto 0,graph 1 nohead lt 0
set label "49" at 7312.0,graph 1.01 center; set arrow from 7310.0,graph 0 rto 0,graph 1 nohead lt 0
set label "50" at 7408.0,graph 1.01 center; set arrow from 7400.0,graph 0 rto 0,graph 1 nohead lt 0
set label "51" at 7507.0,graph 1.01 center; set arrow from 7500.0,graph 0 rto 0,graph 1 nohead lt 0
set label "52" at 7608.0,graph 1.01 center; set arrow from 7600.0,graph 0 rto 0,graph 1 nohead lt 0	
set label "split2" at 5812.0,graph 1.01 center; set arrow from 5812.0,graph 0 rto 0,graph 1 nohead lt 0
set label "split1" at 5178.0,graph 1.01 center; set arrow from 5178.0,graph 0 rto 0,graph 1 nohead lt 0

bind "k" "unset xrange; unset yrange; replot"
bind "1" "set xr [4500:4800]; replot"
bind "2" "set xr [4800:4930]; replot"
bind "3" "set xr [5300:6100]; replot"
bind "4" "set xr [6450:6680]; replot"
bind "5" "set xr [5000:5200]; replot"
bind "6" "set xr [4490:8900]; replot"
bind "7" "set xr [5750:6500]; set yr [0.75:1.80]; replot"

m=1
n=4
indices  = "1  2  3   4   5   6   7  8    9   10
datasets = "17 57 86 122 142 178 197 214 235 269"

p \
 "<awk -v m=" .m. " -v n=" .n. " '($5 >= m && $5 <= n) && ($5!=i) { print s; } ($5 >= m && $5 <= n) { print; i=$5; }' Spectra.dat" u ($2/ang):($3+fac*($5-1)) t 'results' w l lt 3 ,\
 "<awk '($4==235) {print $2*10**10, $3+2*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
 "<awk '($4==197) {print $2*10**10, $3+0*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
 "<awk '($4==214) {print $2*10**10, $3+1*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\

# "<awk '($4==269) {print $2*10**10, $3+9*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
#  "<awk '($4==57) {print $2*10**10, $3+0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '($4==86) {print $2*10**10, $3+2*0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '($4==122) {print $2*10**10, $3+3*0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '($4==142) {print $2*10**10, $3+4*0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '($4==178) {print $2*10**10, $3+5*0.25}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\
#  "<awk '($4==17) {print $2*10**10, $3}' ../synthetic.dat" u 1:2 t 'synthetic' w l lc 'orange',\

pa -1

set term png small
set out "Spectra.png"
rep

q

  "Spectra.dat" u ($2/nm):($3+1*$5):4 w err lt 3 pt 1 ps 0
