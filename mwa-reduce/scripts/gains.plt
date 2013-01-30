set terminal postscript enhanced color
set logscale y
#set xrange [0.1:]
set yrange [0.1:]
set output "gains.ps"
set key bottom left
set xlabel "Channel (index)"
set ylabel "Gain (Jy)"
plot \
"gains.txt" every ::2 using 1:2 with lines title 'A1' lw 2.0, \
"gains.txt" every ::2 using 1:3 with lines title 'A2' lw 2.0, \
"gains.txt" every ::2 using 1:4 with lines title 'A3' lw 2.0, \
"gains.txt" every ::2 using 1:5 with lines title 'A4' lw 2.0, \
"gains.txt" every ::2 using 1:6 with lines title 'A5' lw 2.0, \
"gains.txt" every ::2 using 1:7 with lines title 'A6' lw 2.0, \
"gains.txt" every ::2 using 1:8 with lines title 'A7' lw 2.0, \
"gains.txt" every ::2 using 1:9 with lines title 'A8' lw 2.0, \
"gains.txt" every ::2 using 1:10 with lines title 'A9' lw 2.0, \
"gains.txt" every ::2 using 1:11 with lines title 'A10' lw 2.0, \
"gains-ref.txt" every ::2 using 1:2 with lines title '' lw 2.0 lt 1, \
"gains-ref.txt" every ::2 using 1:3 with lines title '' lw 2.0 lt 2, \
"gains-ref.txt" every ::2 using 1:4 with lines title '' lw 2.0 lt 3, \
"gains-ref.txt" every ::2 using 1:5 with lines title '' lw 2.0 lt 4, \
"gains-ref.txt" every ::2 using 1:6 with lines title '' lw 2.0 lt 5, \
"gains-ref.txt" every ::2 using 1:7 with lines title '' lw 2.0 lt 6, \
"gains-ref.txt" every ::2 using 1:8 with lines title '' lw 2.0 lt 7, \
"gains-ref.txt" every ::2 using 1:9 with lines title '' lw 2.0 lt 8, \
"gains-ref.txt" every ::2 using 1:10 with lines title '' lw 2.0 lt 9, \
"gains-ref.txt" every ::2 using 1:11 with lines title '' lw 2.0 lt 10
