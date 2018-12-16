set terminal x11 0

reset

set grid ytics xtics
set xlabel "x"
set ylabel "y"

set size ratio -1

set key right
plot "endpoint_exp.dat" using 2:3 with line linewidth 5

set terminal svg
set output "position_exp.svg"
replot
set output

set terminal x11 1

reset

set grid ytics xtics
set xlabel "t"
set ylabel "vx"

set key right
plot "endpoint_exp.dat" using 1:4 with line linewidth 5

set terminal svg
set output "vx_exp.svg"
replot
set output

set terminal x11 2

reset

set grid ytics xtics
set xlabel "t"
set ylabel "vy"

set key right
plot "endpoint_exp.dat" using 1:5 with line linewidth 5

set terminal svg
set output "vy_exp.svg"
replot
set output

set terminal x11 3

reset

set grid ytics xtics
set xlabel "t"
set ylabel "ax"

set key right
plot "endpoint_exp.dat" using 1:6 with line linewidth 5

set terminal svg
set output "ax_exp.svg"
replot
set output

set terminal x11 4

reset

set grid ytics xtics
set xlabel "t"
set ylabel "ay"

set key right
plot "endpoint_exp.dat" using 1:7 with line linewidth 5

set terminal svg
set output "ay_exp.svg"
replot
set output

set terminal x11 5

reset

set grid ytics xtics
set xlabel "t"
set ylabel "x, x_o, x_d, x_g of exp"

set key right
plot "endpoint_exp.dat" using 1:2 with line linewidth 5 title "x"
replot "endpoint_exp.dat" using 1:8 with line linewidth 5 title "x_o"
replot "endpoint_exp.dat" using 1:10 with line linewidth 5 title "x_d"
replot "endpoint_exp.dat" using 1:12 with line linewidth 5 title "x_g"

set terminal svg
set output "x_exp.svg"
replot
set output

set terminal x11 6

reset

set grid ytics xtics
set xlabel "t"
set ylabel "y, y_o, y_d, y_g of exp"

set key right
plot "endpoint_exp.dat" using 1:3 with line linewidth 5 title "y"
replot "endpoint_exp.dat" using 1:9 with line linewidth 5 title "y_o"
replot "endpoint_exp.dat" using 1:11 with line linewidth 5 title "y_d"
replot "endpoint_exp.dat" using 1:13 with line linewidth 5 title "y_g"

set terminal svg
set output "y_exp.svg"
replot
set output

set terminal x11 6

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_x, F_x of exp"

set key right
plot "force_exp.dat" using 1:2 with line linewidth 5 title "Fe_x"
replot "force_exp.dat" using 1:4 with line linewidth 5 title "F_x"

set terminal svg
set output "Fx_exp.svg"
replot
set output

set terminal x11 7

reset

set grid ytics xtics
set xlabel "t"
set ylabel "Fe_y, F_y of exp"

set key right
plot "force_exp.dat" using 1:3 with line linewidth 5 title "Fe_y"
replot "force_exp.dat" using 1:5 with line linewidth 5 title "F_y"

set terminal svg
set output "Fy_exp.svg"
replot
set output
