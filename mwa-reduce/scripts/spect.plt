set terminal postscript enhanced color
#set logscale y
#set xrange [0.001:]
#set yrange [-8:2]
set output "spect.ps"
set key bottom left
set xlabel "Frequency (MHz)"
set ylabel "Flux density (Jy)"
plot \
"spectrum.txt" using (column(1)/1000000):2 with lines title '' lw 2.0
