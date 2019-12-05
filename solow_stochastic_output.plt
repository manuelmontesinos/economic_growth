#set multiplot layout 3,2
m="./solow_stochastic_output.txt"
#set tmargin 2

set terminal qt 0
set key
set xlabel "Time"
set title "Capital stock"
plot m using 1:2 with lines title "exact" lt 7, m using 1:3 with lines title "approx" lt -1 dt 2

set terminal qt 1
set key
set xlabel "Time"
set title "Consumption"
plot m using 1:4 with lines title "exact" lt 7, m using 1:5 with lines title "approx" lt -1 dt 2

set terminal qt 2
set key
set xlabel "Time"
set title "Output"
plot m using 1:6 with lines title "exact" lt 7, m using 1:7 with lines title "approx" lt -1 dt 2

set terminal qt 3
set key
set xlabel "Time"
set title "Growth rate of output"
plot m using 1:8 with lines title "exact" lt 7, m using 1:9 with lines title "approx" lt -1 dt 2

set terminal qt 4
set key
set xlabel "Time"
set title "Logged output"
plot m using 1:10 with lines title "exact" lt 7, m using 1:11 with lines title "approx" lt -1 dt 2

set terminal qt 5
set key
set xlabel "Time"
set title "Marginal productivity"
plot m using 1:12 with lines title "exact" lt 7, m using 1:13 with lines title "approx" lt -1 dt 2

#unset multiplot