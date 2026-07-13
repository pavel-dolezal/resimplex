#!/usr/bin/gnuplot

#chi2 = `awk 'BEGIN{ min=1e38; }{ chi=$NF; if (chi<min){ min=chi; } }END{ print min; }' chi2_func.tmp`
nfree = 47

set colors classic
set term x11

set xl "iter"
set yl "chi^2"

#set yr [638000:644000]

set logscale y
set ytics 10
set mytics 10
set key outside samplen 1.5
#set yr [60000:71000]
#set arrow from nfree,graph 1 to nfree,graph 0 lt 0 nohead
#set label sprintf(" %.1f ", chi2) at graph 1,first chi2*1.5 right

plot "<awk '{ i++; print i, $(NF-1); }' chi2_func_presx.tmp" u 1:2 t "chi2" w lp lt 1 ps 0.5,\
#     "<awk '{ i++; print i, $1; }' chi2_func2.tmp" u 1:2 t "alpha" w l ,\
#     "<awk '{ i++; print i, $2; }' chi2_func2.tmp" u 1:2 t "delta" w l ,\
#     "<awk '{ i++; print i, $3; }' chi2_func2.tmp" u 1:2 t "theta" w l ,\
#     "<awk '{ i++; print i, $4; }' chi2_func2.tmp" u 1:2 t "max_blaze" w l ,\
#     "<awk '{ i++; print i, $5; }' chi2_func2.tmp" u 1:2 t "lambda-c" w l ,\
#     "<awk '{ i++; print i, $3; }' chi2_func2.tmp" u 1:2 t "epsilon" w l ,\
#     "<awk '{ i++; print i, $8; }' chi2_func.tmp" u 1:2 t "lambda_c" w l,\
# "<awk '{ i++; print i, $6; }' chi2_func.tmp" u 1:2 t "theta" w l,\

#  chi2 w l lt 0 not

pa -1

replot

q
