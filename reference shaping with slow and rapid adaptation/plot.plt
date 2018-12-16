set terminal x11 0

reset

set grid ytics xtics
# set xrange [0:1]
# set yrange [0:0.4]
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
# set xrange [0:2.5]
# set yrange [-0.5:0.1]
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
# set xrange [0:2.5]
# set yrange [0:0.5]
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
# set xrange [0:2.5]
# set yrange [-1.5:1]
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
# set xrange [0:2.5]
# set yrange [5:-1]
set xlabel "t"
set ylabel "ay"

set key right
plot "endpoint.dat" using 1:7 with line linewidth 5

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "ay.eps"
replot
set output

set terminal x11 5

reset

set grid ytics xtics
set xlabel "t"
set ylabel "x, x_o, x_d, x_g"

set key right
plot "endpoint.dat" using 1:2 with line linewidth 5 title "x"
replot "endpoint.dat" using 1:8 with line linewidth 5 title "x_o"
replot "endpoint.dat" using 1:10 with line linewidth 5 title "x_d"
replot "endpoint.dat" using 1:12 with line linewidth 5 title "x_g"

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "x, x_o, x_d, x_g.eps"
replot
set output

set terminal x11 6

reset

set grid ytics xtics
set xlabel "t"
set ylabel "y, y_o, y_d, y_g"

set key right
plot "endpoint.dat" using 1:3 with line linewidth 5 title "y"
replot "endpoint.dat" using 1:9 with line linewidth 5 title "y_o"
replot "endpoint.dat" using 1:11 with line linewidth 5 title "y_d"
replot "endpoint.dat" using 1:13 with line linewidth 5 title "y_g"

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "y, y_o, y_d, y_g.eps"
replot
set output

set terminal x11 7

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_x, F_x"

set key right
plot "force.dat" using 1:2 with line linewidth 5 title "Fe_x"
replot "force.dat" using 1:4 with line linewidth 5 title "F_x"

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "Fe_x, F_x.eps"
replot
set output

set terminal x11 8

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_y, F_y"

set key right
plot "force.dat" using 1:3 with line linewidth 5 title "Fe_y"
replot "force.dat" using 1:5 with line linewidth 5 title "F_y"

set terminal postscript eps enhanced color "GothicBBB-Medium-EUC-H"
set output "Fe_y, F_y.eps"
replot
set output
