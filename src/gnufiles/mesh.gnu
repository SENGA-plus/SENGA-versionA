#
# gnuplot command file for plotting flames
#
#
#
#set xrange [0:1.0E-2]
#set yrange [-1.0:1.0]
set style data lines
#
plot "meshpoints.res"
pause 10
#
plot "dhatdg.res"
pause 10
#
plot "strpoints.res"
pause 10
#
plot "dgdhat.res"
pause 10
#
plot "d2gdh2.res"
pause 10
#
exit
