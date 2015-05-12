echo "set grid" > gpfile
echo "set xrange [0:100]" > gpfile
plot="plot \"result.dat\" u 1:2 w l"
echo $plot >> gpfile
tail -f  gpfile | gnuplot  '-'
