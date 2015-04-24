#!/usr/bin/gnuplot
set terminal wxt 1 size 1680,1050 persist
set multiplot layout 5,1
set xrange [1:92]
set lmargin at screen 0.05
plot '../OUTPUT/OUTPUT_MM.DAT' using 1:2 with lines title 'MM1'
plot '../OUTPUT/OUTPUT_MM.DAT' using 1:3 with lines title 'MM2'
plot '../OUTPUT/OUTPUT_MM.DAT' using 1:4 with lines title 'MM3'
plot '../OUTPUT/OUTPUT_MM.DAT' using 1:($2+$3+$4) with lines title 'ISAM'
stats "< cat ../OUTPUT/OUTPUT_MM.DAT | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 '" using 1:(($2+$3+$4)/50)
m1Max = STATS_mean_y
stats 'dataK.dat' using 1:($4/$3)
eMax = STATS_mean_y
plot 'dataK.dat' using 1:($4/$3/eMax) with points title 'Exp.', '../OUTPUT/OUTPUT_MM.DAT' using 1:(($2+$3+$4)/50/m1Max) with lines title 'Model'
unset multiplot
