#!/usr/bin/gnuplot
set datafile missing '***'
stats 'dataK.dat' using 1:($4/$3)
eCsI = STATS_mean_y
stats 'result_csi.dat' using 1:($5/$4)
mCsI = STATS_mean_y
stats 'dataK.dat' using 1:($6/$3)
eGe = STATS_mean_y
stats 'result_ge.dat' using 1:($5/$4)
mGe = STATS_mean_y
stats 'dataK.dat' using 1:($8/$3)
eKBr = STATS_mean_y
stats 'result_kbr.dat' using 1:($5/$4)
mKBr = STATS_mean_y
stats 'dataK.dat' using 1:($10/$3)
eMo = STATS_mean_y
stats 'result_mo.dat' using 1:($5/$4)
mMo = STATS_mean_y
stats 'dataK.dat' using 1:($12/$3)
eW = STATS_mean_y
stats 'result_w.dat' using 1:($5/$4)
mW = STATS_mean_y
set terminal pdfcairo size 21cm, 20cm
#set terminal wxt size 1026,768 persist
set xrange [10:80]
set xlabel 'E(keV)'
set ylabel 'I(cps/10keV)'
set output 'plot.pdf'
set title 'CsI'
plot 'dataK.dat' using 1:(($4/$3)/(eCsI)) with points title 'Experimenteel', 'result_csi.dat' using 1:(($5/$4)/(mCsI)) with lines title 'Model'
set title 'Ge'
plot 'dataK.dat' using 1:(($6/$3)/(eGe)) with points title 'Experimenteel' , 'result_ge.dat' using 1:(($5/$4)/(mGe)) with lines title 'Model'
set title 'KBr'
plot 'dataK.dat' using 1:(($8/$3)/(eKBr)) with points title 'Experimenteel', 'result_kbr.dat' using 1:(($5/$4)/(mKBr)) with lines title 'Model'
set title 'Mo'
plot 'dataK.dat' using 1:(($10/$3)/(eMo)) with points title 'Experimenteel', 'result_mo.dat' using 1:(($5/$4)/(mMo)) with lines title 'Model'
set title 'W'
plot 'dataK.dat' using 1:(($12/$3)/(eW)) with points title 'Experimenteel', 'result_w.dat' using 1:(($5/$4)/(mW)) with lines title 'Model'
set terminal wxt
replot
