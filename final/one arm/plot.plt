set terminal x11 0

reset

set grid ytics xtics
set xlabel "x"
set ylabel "y"

set size ratio -1

set key right
plot "endpoint.dat" using 2:3 with line linewidth 5

set terminal svg
set output "position.svg"
replot
set output

set terminal x11 1

reset

set grid ytics xtics
set xlabel "t"
set ylabel "vx"

set key right
plot "endpoint.dat" using 1:4 with line linewidth 5

set terminal svg
set output "vx.svg"
replot
set output

set terminal x11 2

reset

set grid ytics xtics
set xlabel "t"
set ylabel "vy"

set key right
plot "endpoint.dat" using 1:5 with line linewidth 5

set terminal svg
set output "vy.svg"
replot
set output

set terminal x11 3

reset

set grid ytics xtics
set xlabel "t"
set ylabel "ax"

set key right
plot "endpoint.dat" using 1:6 with line linewidth 5

set terminal svg
set output "ax.svg"
replot
set output

set terminal x11 4

reset

set grid ytics xtics
set xlabel "t"
set ylabel "ay"

set key right
plot "endpoint.dat" using 1:7 with line linewidth 5

set terminal svg
set output "ay.svg"
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

set terminal svg
set output "x.svg"
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

set terminal svg
set output "y.svg"
replot
set output

set terminal x11 6

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_x, F_x"

set key right
plot "force.dat" using 1:2 with line linewidth 5 title "Fe_x"
replot "force.dat" using 1:4 with line linewidth 5 title "F_x"

set terminal svg
set output "Fx.svg"
replot
set output

set terminal x11 7

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_y, F_y"

set key right
plot "force.dat" using 1:3 with line linewidth 5 title "Fe_y"
replot "force.dat" using 1:5 with line linewidth 5 title "F_y"

set terminal svg
set output "Fy.svg"
replot
set output
