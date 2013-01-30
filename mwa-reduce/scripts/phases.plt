set terminal postscript enhanced color
#set logscale y
#set xrange [0.001:]
#set yrange [-8:2]
set output "phases.ps"
set key bottom left
set xlabel "Channel (index)"
set ylabel "Phase offset (deg)"
plot \
"phases.txt" every ::1 using 1:(column(2)*180/3.1415) with lines title 'A1' lw 2.0, \
"phases.txt" every ::1 using 1:(column(3)*180/3.1415) with lines title 'A2' lw 2.0, \
"phases.txt" every ::1 using 1:(column(4)*180/3.1415) with lines title 'A3' lw 2.0, \
"phases.txt" every ::1 using 1:(column(5)*180/3.1415) with lines title 'A4' lw 2.0, \
"phases.txt" every ::1 using 1:(column(6)*180/3.1415) with lines title 'A5' lw 2.0, \
"phases.txt" every ::1 using 1:(column(7)*180/3.1415) with lines title 'A6' lw 2.0, \
"phases.txt" every ::1 using 1:(column(8)*180/3.1415) with lines title 'A7' lw 2.0, \
"phases.txt" every ::1 using 1:(column(9)*180/3.1415) with lines title 'A8' lw 2.0, \
"phases.txt" every ::1 using 1:(column(10)*180/3.1415) with lines title 'A9' lw 2.0, \
"phases.txt" every ::1 using 1:(column(11)*180/3.1415) with lines title 'A10' lw 2.0, \
"phases.txt" every ::1 using 1:((column(2)*180/3.1415 + column(3)*180/3.1415 + column(4)*180/3.1415 + column(5)*180/3.1415 + column(6)*180/3.1415 + column(7)*180/3.1415 + column(8)*180/3.1415 + column(9)*180/3.1415 + column(10)*180/3.1415 + column(11)*180/3.1415 + column(12)*180/3.1415 + column(13)*180/3.1415 + column(14)*180/3.1415 + column(15)*180/3.1415 + column(16)*180/3.1415 + column(17)*180/3.1415 + column(18)*180/3.1415 + column(19)*180/3.1415 + column(20)*180/3.1415 + column(21)*180/3.1415 + column(22)*180/3.1415 + column(23)*180/3.1415 +column(24)*180/3.1415 + column(25)*180/3.1415 + column(26)*180/3.1415 + column(27)*180/3.1415 + column(28)*180/3.1415 + column(29)*180/3.1415 + column(30)*180/3.1415 + column(31)*180/3.1415 + column(32)*180/3.1415 + column(33)*180/3.1415 )/32) with lines title 'Avg' lw 4.0, \
"phases-ref.txt" every ::1 using 1:(column(2)*180/3.1415) with lines title '' lw 2.0 lt 1, \
"phases-ref.txt" every ::1 using 1:(column(3)*180/3.1415) with lines title '' lw 2.0 lt 2, \
"phases-ref.txt" every ::1 using 1:(column(4)*180/3.1415) with lines title '' lw 2.0 lt 3, \
"phases-ref.txt" every ::1 using 1:(column(5)*180/3.1415) with lines title '' lw 2.0 lt 4, \
"phases-ref.txt" every ::1 using 1:(column(6)*180/3.1415) with lines title '' lw 2.0 lt 5, \
"phases-ref.txt" every ::1 using 1:(column(7)*180/3.1415) with lines title '' lw 2.0 lt 6, \
"phases-ref.txt" every ::1 using 1:(column(8)*180/3.1415) with lines title '' lw 2.0 lt 7, \
"phases-ref.txt" every ::1 using 1:(column(9)*180/3.1415) with lines title '' lw 2.0 lt 8, \
"phases-ref.txt" every ::1 using 1:(column(10)*180/3.1415) with lines title '' lw 2.0 lt 9, \
"phases-ref.txt" every ::1 using 1:(column(11)*180/3.1415) with lines title '' lw 2.0 lt 10
