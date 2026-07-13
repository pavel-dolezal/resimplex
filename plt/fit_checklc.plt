#!/usr/bin/gnuplot

set colors classic
set term x11

plot for [x=2:10] "<awk -v x=". x. " '$NF == x' fit_checklc.asc" u 1:($2-$3) w l not lc 'red' ,\
	"<awk '$NF == 1' fit_checklc.asc" u 1:($2-$3) w l t 'residua' lc 'red',\

#plot for [x=2:10] "<awk -v x=". x. " '$NF == x' fit_checkA.asc" u 1:2 w p not lc 'red',\

#	"fit_checklc.asc" u 1:3 w lp lc 'blue',\

pa -1


