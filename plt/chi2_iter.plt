#!/usr/bin/gnuplot

chi2 = `awk 'BEGIN{ min=1e38; }{ chi=$(NF-1); if (chi<min){ min=chi; } }END{ print min; }' chi2_func.tmp`
nfree = 47

set colors classic
set term x11

set xl "iter"
set yl "chi^2"

#set yr [2e5:1e7]
#set yr [820000:836000]
set logscale y
set ytics 10
set mytics 10
set key outside samplen 1.5

set arrow from nfree,graph 1 to nfree,graph 0 lt 0 nohead
set label sprintf(" %.1f ", chi2) at graph 1,first chi2*1.5 right

plot "<awk '{ i++; print i, $(NF-1); }' chi2_func.tmp"      u 1:2 t "chi2" w lp lt 1 ps 0.5,\
	chi2 w l lt 0 not
#     "<awk '{ i++; print i, sqrt(($1)^2); }' chi2_func.tmp" u 1:2 t "a(1)" w l ,\
#     "<awk '{ i++; print i, $2; }' chi2_func.tmp" u 1:2 t "a(2)" w l ,\
#     "<awk '{ i++; print i, sqrt(($3)^2); }' chi2_func.tmp" u 1:2 t "a(3)" w l ,\
#     "<awk '{ i++; print i, sqrt(($4)^2); }' chi2_func.tmp" u 1:2 t "a(4)" w l ,\
#     "<awk '{ i++; print i, sqrt(($5)^2); }' chi2_func.tmp" u 1:2 t "a(5)" w l ,\
#     "<awk '{ i++; print i, sqrt(($6)^2); }' chi2_func.tmp" u 1:2 t "a(6)" w l ,\
#     "<awk '{ i++; print i, $8; }' chi2_func.tmp" u 1:2 t "lambda_c" w l,\
# "<awk '{ i++; print i, $6; }' chi2_func.tmp" u 1:2 t "theta" w l,\



pa -1

replot

q
