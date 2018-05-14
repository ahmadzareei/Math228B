
All cpp extension files in here, contain the written code in c++

all codes are compiled and runned under gnu compiler (g++)
"gnuplot" is used for making the plots 

for compiling, under linux machine, simply write in the terminal:
>>g++ filename.cpp -o output.o

for running
>>./output.o

________________________________________________________

Gnuplot Scripts that are used you can find in here:

* For using these scripts, in your current directory make 3 folders with names "N512", "N1024" and "N2048" and change N in these codes and run the code in these folders, then open the terminal in your current folder and run the fillowing scripts:

//****************************************************
#Finding energy (plot of energy vs time)

set terminal postscript eps enhanced color solid lw 2
set out 'energy.eps'
set xlabel 't'
set ylabel 'Energy'
set yrange [0:3]
plot 'N2048/energy_time.dat' w l t 'N=2048'

//************************************************
#Finding propagating shape of wave over time
set terminal postscript eps enhanced color solid lw 2
set out 'shape.eps'

set xlabel 'x'
set ylabel 'u'
set xtics rotate
set multiplot layout 3,3

plot 'N2048/initial_shape.dat' w l t ''
plot 'N2048/data_1.dat' w l t ''
plot 'N2048/data_2.dat' w l t ''
plot 'N2048/data_3.dat' w l t ''
plot 'N2048/data_4.dat' w l t ''
plot 'N2048/data_5.dat' w l t ''
plot 'N2048/data_6.dat' w l t ''
plot 'N2048/data_7.dat' w l t ''
plot 'N2048/data_8.dat' w l t ''
//**************************************************
# convergence test for the scheme... 
set terminal postscript eps enhanced color solid lw 2
set out 'shape_time8.eps'

set xlabel 'x'
set ylabel 'u'

plot 'N2048/data_79.dat' w l t 'N=2048', 'N1024/data_79.dat' w l t 'N=1024' , 'N512/data_79.dat' w l t 'N=512'
 
//******************************************************
