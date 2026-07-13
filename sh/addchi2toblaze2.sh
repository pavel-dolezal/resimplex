#!/bin/bash

awk 'FNR==NR {a[NR]=$2; next} {print $1,FS, $3,FS, $4,FS, a[$1], FS, $5}' chi2_order.tmp blaze2.tmp  >blaze2_chi2.tmp


