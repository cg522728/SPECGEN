#!/usr/bin/gnuplot
set terminal wxt 1 persist size 1024,768
stats "< cat ../OUTPUT/OUTPUT_MM.DAT | grep -E '12 |14 |17 |19 |22 |23 |24 |25 |26 |28 |29 |32 |34 |45 |46 |47 |50 |79 '" using 1:($2+$3+$4)
m1Max = STATS_mean_y
stats 'dataK.dat' using 1:($4/$3)
eMax = STATS_mean_y
plot 'dataK.dat' using 1:($4/$3/eMax) with points
replot '../OUTPUT/OUTPUT_MM.DAT' using 1:(($2+$3+$4)/m1Max) with linespoints


