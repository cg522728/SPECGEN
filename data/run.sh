#!/bin/bash
rm OUTPUT_MM.DAT
rm result*.dat
rm sort_*.dat
touch result_csi.dat
for i in `seq 12 79`;
do
clear
echo '#####################################################'
echo $i '/74'
echo 'CsI'
./mmsens 1 100. 6. 'CsI' 4 0. $i 50. 0.1 0.001
cat OUTPUT_MM.DAT >> result_csi.dat
#./plot.p
#clear
#echo '#####################################################'
#echo $i '/74'
#echo 'Ge'
#./mmsens 1 100. 6. 'Ge' 4 0. $i 50. 0.01 0.001
#cat OUTPUT_MM.DAT >> result_ge.dat
#cat result_ge.dat | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 ' > sort_ge.dat
#./plot.p
#clear
#echo '#####################################################'
#echo $i '/74'
#echo 'KBr'
#./mmsens 1 100. 6. 'KBr' 4 0. $i 50. 0.01 0.001
#cat OUTPUT_MM.DAT >> result_kbr.dat
#cat result_kbr.dat | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 ' > sort_kbr.dat
#./plot.p
#clear
#echo '#####################################################'
#echo $i '/74'
#echo 'Mo'
#./mmsens 1 100. 6. 'Mo' 4 0. $i 50. 0.01 0.001
#cat OUTPUT_MM.DAT >> result_mo.dat
#cat result_mo.dat | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 ' > sort_mo.dat
#./plot.p
#clear
#echo '#####################################################'
#echo $i '/74'
#echo 'W'
#./mmsens 1 100. 6. 'W' 4 0. $i 50. 0.01 0.001
#cat OUTPUT_MM.DAT >> result_w.dat
#cat result_w.dat | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 ' > sort_w.dat
#./plot.p
done





