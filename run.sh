#!/bin/bash
rm result.dat
touch result.dat
for i in `seq 10 80`;
do
clear
./SPECGEN 1 100. 6. 'CsI' 4 0. $i 50. 0.01 1.
cat OUTPUT.DAT >> result.dat
done
