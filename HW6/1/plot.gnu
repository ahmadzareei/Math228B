 set terminal postscript eps enhanced color solid fontscale '0.5'
set out 'particle_epsilon_0.0.eps'

set xlabel 'x'
set ylabel 'y'
set yrange [-2.5:1.5]
set xrange [-0.5:1.5]
set title 'Boundary at different times'
plot '1_a/initial.dat' w l t 't=0', '1_a/middle.dat' w l t 't=0.05','1_a/middle2.dat' w l t 't=0.1', '1_a/final.dat' w l t 't=0.15'