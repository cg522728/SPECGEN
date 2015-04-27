#!/bin/bash
let v=100
let i=6
mkdir tmp
cd tmp
../../OUTPUT/tspec 1 $v. $i. 'CsI' 4 0. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../nofilter.dat
../../OUTPUT/tspec 1 $v. $i. 'CsI' 42 250. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../mo250.dat
../../OUTPUT/tspec 1 $v. $i. 'CsI' 40 125. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../zr125.dat
../../OUTPUT/tspec 1 $v. $i. 'CsI' 29 100. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../cu100.dat
../../OUTPUT/tspec 1 $v. $i. 'CsI' 13 500. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../al500.dat
../../OUTPUT/tspec 1 $v. $i. 'CsI' 13 100. 'Fe' 1. 0.1 0.
cat OUTPUT_T.DAT > ../al100.dat
cd ..
rm -r tmp
./filter.p
echo $v
