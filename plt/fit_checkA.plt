#!/usr/bin/gnuplot

set colors classic
set term x11

set key bottom
set xl "r"
set yl "{/Symbol a}"


plot for [x=2:10] "<awk -v x=". x. " '$NF == x' fit_checkA.asc" u 1:2 w p not lc 'red',\
     for [y=2:10] "<awk -v y=". y. " '$NF == y' fit_checkA.asc" u 1:3 w l not lc 'blue',\
     "<awk '$NF==1' fit_checkA.asc" u 1:2 w p t 'presimplex free' lc 'red' ,\
     "<awk '$NF==1' fit_checkA.asc" u 1:3 w l t 'resimplex fit' lc 'blue' ,\

pa -1


