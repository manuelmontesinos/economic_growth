cd ''
m="./change_savings_output.txt"

set terminal qt 0
set key right bottom
set xlabel "Time"
set title "Effect of a change in the savings rate on the capital stock"
plot m using 1:2 with lines title "exact solution" lt 7 lc rgb 'blue', m using 1:3 with lines title "linear approx" lt -1 dt 2 lc rgb 'red', m using 1:4 with lines title "quadratic approx" lt -1 dt 2 lc rgb 'gold', m using 1:5 with lines title "exact continuous" lt -1 dt 2 lc rgb 'violet'

set terminal qt 1
set key right bottom
set xlabel "Time"
set title "Effect of a change in the savings rate on consumption"
plot m using 1:6 with lines title "exact solution" lt 7 lc rgb 'blue', m using 1:7 with lines title "linear approx" lt -1 dt 2 lc rgb 'red', m using 1:8 with lines title "quadratic approx" lt -1 dt 2 lc rgb 'gold', m using 1:9 with lines title "exact continuous" lt -1 dt 2 lc rgb 'violet'

set terminal qt 2
set key right bottom
set xlabel "Time"
set title "Effect of a change in the savings rate on output"
plot m using 1:10 with lines title "exact solution" lt 7 lc rgb 'blue', m using 1:11 with lines title "linear approx" lt -1 dt 2 lc rgb 'red', m using 1:12 with lines title "quadratic approx" lt -1 dt 2 lc rgb 'gold', m using 1:13 with lines title "exact continuous" lt -1 dt 2 lc rgb 'violet'

set terminal qt 3
set key right bottom
set xlabel "Time"
set title "Effect of a change in the savings rate on the growth rate of output"
plot m using 1:14 with lines title "exact solution" lt 7 lc rgb 'blue', m using 1:15 with lines title "linear approx" lt -1 dt 2 lc rgb 'red', m using 1:16 with lines title "quadratic approx" lt -1 dt 2 lc rgb 'gold', m using 1:17 with lines title "exact continuous" lt -1 dt 2 lc rgb 'violet'

set terminal qt 4
set key right bottom
set xlabel "Time"
set title "Effect of a change in the savings rate on the log of output"
plot m using 1:18 with lines title "exact solution" lt 7 lc rgb 'blue', m using 1:19 with lines title "linear approx" lt -1 dt 2 lc rgb 'red', m using 1:20 with lines title "quadratic approx" lt -1 dt 2 lc rgb 'gold', m using 1:21 with lines title "exact continuous" lt -1 dt 2 lc rgb 'violet'
