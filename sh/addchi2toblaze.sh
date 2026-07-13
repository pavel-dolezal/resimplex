#!/bin/bash

awk 'FNR==NR {a[NR]=$2; next} {print $1,FS, $3,FS, $4,FS, a[$1], FS, $NF}' chi2_order.tmp blaze1.tmp  >blaze1_chi2.tmp


