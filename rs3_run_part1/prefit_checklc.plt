#!/usr/bin/gnuplot

set colors classic
set term x11

#a1 = 0.951467966477249 
#a2 = -0.00310883933775447 
#a3 = 0.000289276989305378 
#a4 = -7.4296819634584e-06 
#a5 = 6.75457672210317e-08
#f(x) = a1+a2*x+a3*x**2+a4*x**3+a5*x**4
#
#fit f(x) 'prefit_checklc.tmp' u 1:2 via a1,a2,a3,a4,a5

plot \
    for [dataset=2:10] sprintf("< awk -v dataset=%d '$4==dataset' prefit_checklc.tmp", dataset) \
        using 1:2 not w lp lc 'red', \
    "< awk '$4==1' prefit_checklc.tmp" \
        using 1:2 title 'free max-blaze' with lp lc 'red', \
    for [dataset=2:10] sprintf("< awk -v dataset=%d '$4==dataset' prefit_checklc.tmp", dataset) \
        u 1:3 not w lp lc 'blue' ,\
	"<awk '$4==1' prefit_checklc.tmp" u 1:2 t 'resimplex prefit' w lp lc 'blue',\

#	"prefit_checklc.tmp" u 1:3 t 'resimplex-prefit' w l lc 'blue',\
#	"prefit_checklc.tmp" u 1:(f($1)) t 'gnuplot prefit' w l lc 'black' ,\

pa -1


