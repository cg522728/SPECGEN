#!/usr/bin/gnuplot
set datafile missing 'NaN'
set terminal dumb
plot 'data.dat' using 1:($5/$4)
ymax2 = GPVAL_DATA_Y_MAX
plot 'result.dat' using 1:($5/$4)
ymax3 = GPVAL_DATA_Y_MAX
set terminal wxt persist
plot 'data.dat' using 1:(($5/$4)/ymax2) with points title 'Experimenteel'
replot 'result.dat' using 1:(($5/$4)/ymax3) with linespoints title 'Model'
