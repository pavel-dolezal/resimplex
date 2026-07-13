#!/bin/sh

./Spectra.awk list ktc*.asc > Spectra.tmp
sort -g < Spectra.tmp > Spectra.dat
rm Spectra.tmp

