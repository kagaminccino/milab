set terminal x11 0

reset

set grid ytics xtics
set xrange [0:1]
set yrange [0:0.4]
set xlabel "x"
set ylabel "y"

set size ratio -1

set key right
plot "endpoint.dat" using 2:3 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "position.eps"
replot
set output

set terminal x11 1

reset

set grid ytics xtics
set xrange [0:2.5]
set yrange [-0.5:0.1]
set xlabel "t"
set ylabel "vx"

set key right
plot "endpoint.dat" using 1:4 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "vx.eps"
replot
set output

set terminal x11 2

reset

set grid ytics xtics
set xrange [0:2.5]
set yrange [0:0.5]
set xlabel "t"
set ylabel "vy"

set key right
plot "endpoint.dat" using 1:5 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "vy.eps"
replot
set output

set terminal x11 3

reset

set grid ytics xtics
set xrange [0:2.5]
set yrange [-1.5:1]
set xlabel "t"
set ylabel "ax"

set key right
plot "endpoint.dat" using 1:6 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "ax.eps"
replot
set output

set terminal x11 4

reset

set grid ytics xtics
set xrange [0:2.5]
set yrange [5:-1]
set xlabel "t"
set ylabel "ay"

set key right
plot "endpoint.dat" using 1:7 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "ay.eps"
replot
set output
