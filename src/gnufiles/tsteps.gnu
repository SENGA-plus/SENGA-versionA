#
# gnuplot command file for plotting timesteps
#
#
#
#set xrange [0:1000]
#set yrange [0:1E-6]
set data style lines
plot "fort.3"
pause 10
#set term postscript eps
#set output "ranfst11.eps"
#plot "ranfst11.res"
#
exit
