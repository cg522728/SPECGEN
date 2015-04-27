#!/usr/bin/gnuplot
set terminal pdfcairo size 28cm,18cm
set output 'filter.pdf'
load 'palette.pal'
plot 'nofilter.dat'  using 1:2 with lines title 'No filter',\
	'mo250.dat' using 1:2 with lines title 'Mo, 250um',\
	'zr125.dat' using 1:2 with lines title 'Zr, 125um',\
	'cu100.dat' using 1:2 with lines title 'Cu, 100um',\
	'al500.dat' using 1:2 with lines title 'Al, 500um',\
	'al100.dat' using 1:2 with lines title 'Al, 100um'
unset output
set terminal wxt
replot
