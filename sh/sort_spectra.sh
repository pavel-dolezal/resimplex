#!/bin/bash

for file in ktc*.asc; do 
	sort -g -o "sorted_spectra/$file" "$file"
		done
