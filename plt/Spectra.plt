#!/usr/bin/gnuplot

set colors classic
set term qt

ang = 1.e-10  # m
nm = 1.e-9  # m

set xl "lambda [ang]"
set yl "I_lambda [] (shifted by dataset number)"

#set xr [:8000]
#set yr [0.25:3.5]
#set yr [1.75:3.5]
#set yr [2.25:3.5]
#set xr [5650:5915]
#set xr [4490:8900]
#set xr [5000:6000]
#set xr [5881:6585]
#set yr [2.57:2.81]

folder_name = system("pwd")
set title "folder:" . folder_name

fac = 0.25
set ytics fac
set mytics 1
set grid ytics mytics
set zeroaxis
set bar 0.5

set label "48" at 7212.0,graph 1.01 center; set arrow from 7212.0,graph 0 rto 0,graph 1 nohead lt 0
set label "49" at 7310.0,graph 1.01 center; set arrow from 7310.0,graph 0 rto 0,graph 1 nohead lt 0
set label "50" at 7400.0,graph 1.01 center; set arrow from 7400.0,graph 0 rto 0,graph 1 nohead lt 0
set label "51" at 7500.0,graph 1.01 center; set arrow from 7500.0,graph 0 rto 0,graph 1 nohead lt 0
set label "52" at 7600.0,graph 1.01 center; set arrow from 7600.0,graph 0 rto 0,graph 1 nohead lt 0	
set label "split2" at 5812.0,graph 1.01 center; set arrow from 5812.0,graph 0 rto 0,graph 1 nohead lt 0
set label "split1" at 5178.0,graph 1.01 center; set arrow from 5178.0,graph 0 rto 0,graph 1 nohead lt 0

bind "k" "unset xrange; unset yrange; replot"
bind "1" "set xr [4500:4800]; replot"
bind "2" "set xr [4800:4930]; replot"
bind "3" "set xr [5300:6100]; replot"
bind "4" "set xr [6450:6680]; replot"
bind "5" "set xr [5000:5200]; replot"
bind "6" "set xr [4490:8900]; replot"
bind "7" "set xr [6150:6525]; replot"

m=1
n=30

do for [j=m:n] {
	set label sprintf("%d",j) at graph 0.9, first 1.07+(j-1)*fac;
}

p \
 "<awk -v m=" .m. " -v n=" .n. " '($5 >= m && $5 <= n) && ($5!=i) { print s; } ($5 >= m && $5 <= n) { print; i=$5; }' Spectra.dat" u ($2/ang):($3+fac*($5-1)) t 'results' w l lt 3 ,\
 "<awk '($4==18) {print $2, $3+0*0.25}' ../../common_input/synthetic.dat" u ($1/ang):2 t 'synthetic-example' w l lc 'orange',\
 "<awk '($5==1) {print $2, $(NF-1)*0.5}' Spectra.dat" u ($1/ang):2 w l t 'interpolation' lc 'red' ,\
# "<awk '($4==18) {print $2*10**10, $3+0*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==8) {print $2*10**10, $3+1*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==7) {print $2*10**10, $3+2*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==9) {print $2*10**10, $3+3*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==235) {print $2*10**10, $3+4*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==197) {print $2*10**10, $3+5*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\
# "<awk '($4==214) {print $2*10**10, $3+6*0.25}' ../synthetic.dat" u 1:2 not w l lc 'orange',\

pa -1

#set term png small
#set out "Spectra.png"
#rep

q
