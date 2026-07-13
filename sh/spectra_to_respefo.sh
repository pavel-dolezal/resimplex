#!/bin/bash

for file in ktc*.asc; do
	awk '/^[[:space:]]*#/ {print $0; next} { print $1, $2}' "$file" > spectra_to_respefo/"$file"
done
