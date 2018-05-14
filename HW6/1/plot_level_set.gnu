set terminal postscript eps enhanced color solid fontscale '0.5'
set out 'level_set_0.0.eps'

#set xlabel 'x'
#set ylabel 'y'
set yrange [-2.5:1.5]
set xrange [-0.5:1.5]

set contour
unset surface
set view 0,0,1

set cntrparam levels discrete 0
unset key
splot '1_b/initial.dat' w l t 't=0', '1_b/middle.dat' w l t 't=0.05','1_b/middle2.dat' w l t 't=0.1','1_b/middle3.dat' w l t 't=0.15', '1_b/final.dat' w l t 't=0.2'